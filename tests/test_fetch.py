"""Tests for downloads over several connections: segments of one file, resuming
them (and a .part left by a single-connection download), falling back to one
connection, servers that say they are busy, and several files at once."""

import hashlib
import http.server
import json
import re
import threading
import time

import pytest

from data_steward import cli, fetch as fetch_mod, sources
from data_steward.db import Catalog
from data_steward.fetch import FetchError, fetch, fetch_many, pending_bytes, retry_after

SIZE = 300_000


def blob(seed: int, size: int = SIZE) -> bytes:
    out, i = bytearray(), 0
    while len(out) < size:
        out += hashlib.sha256(f"{seed}-{i}".encode()).digest()
        i += 1
    return bytes(out[:size])


class RangeServer:
    """A local HTTP server for in-memory files that answers Range requests
    (unless ranges=False), and can be told to say it is busy, cut responses
    short, or answer slowly. It records the requests and how many it served
    at once."""

    def __init__(self, files: dict, ranges: bool = True):
        self.files, self.ranges = files, ranges
        self.log = []          # (path, Range header or None, status)
        self.busy = 0          # answer this many ranged requests (not the 1-byte probe) with 429
        self.cut_after = None  # stop each response after this many bytes
        self.delay = 0.0       # pause between 4 KiB pieces
        self.active = self.max_active = 0
        self.lock = threading.Lock()
        outer = self

        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def do_HEAD(self):
                self.answer(head=True)

            def do_GET(self):
                self.answer(head=False)

            def answer(self, head):
                path = self.path.split("?")[0].lstrip("/")
                data = outer.files.get(path)
                rng = self.headers.get("Range")
                if data is None:
                    return self.status(path, rng, 404)
                start, end, code = 0, len(data) - 1, 200
                m = re.fullmatch(r"bytes=(\d+)-(\d*)", rng or "")
                if m and outer.ranges:
                    start = int(m[1])
                    end = min(int(m[2]) if m[2] else len(data) - 1, len(data) - 1)
                    if start >= len(data):
                        return self.status(path, rng, 416)
                    with outer.lock:
                        busy = outer.busy > 0 and not (start == 0 and end == 0)
                        outer.busy -= busy
                    if busy:
                        return self.status(path, rng, 429, {"Retry-After": "0"})
                    code = 206
                with outer.lock:
                    outer.log.append((path, rng, code))
                self.send_response(code)
                self.send_header("Content-Length", str(end - start + 1))
                self.send_header("Accept-Ranges", "bytes" if outer.ranges else "none")
                if code == 206:
                    self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
                self.end_headers()
                if head:
                    return
                body = data[start:end + 1]
                if outer.cut_after is not None:
                    body = body[:outer.cut_after]
                with outer.lock:
                    outer.active += 1
                    outer.max_active = max(outer.max_active, outer.active)
                try:
                    for i in range(0, len(body), 4096):
                        self.wfile.write(body[i:i + 4096])
                        if outer.delay:
                            time.sleep(outer.delay)
                except (BrokenPipeError, ConnectionResetError):
                    pass
                finally:
                    with outer.lock:
                        outer.active -= 1
                    self.close_connection = True

            def status(self, path, rng, code, headers=None):
                with outer.lock:
                    outer.log.append((path, rng, code))
                self.send_response(code)
                for k, v in (headers or {}).items():
                    self.send_header(k, v)
                self.send_header("Content-Length", "0")
                self.end_headers()

        self.httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.httpd.daemon_threads = True
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"

    def remote(self, name: str, **kw) -> sources.RemoteFile:
        data = self.files[name]
        return sources.RemoteFile(name, f"{self.url}/{name}", len(data), hashlib.md5(data).hexdigest(), **kw)

    def served_ranges(self, name: str) -> list[tuple[int, int]]:
        out = []
        for path, rng, code in self.log:
            if path == name and code == 206 and rng != "bytes=0-0":
                a, b = rng.removeprefix("bytes=").split("-")
                out.append((int(a), int(b)))
        return out


@pytest.fixture(autouse=True)
def small_segments(monkeypatch):
    """Segment small test files the way large real ones are."""
    monkeypatch.setattr(fetch_mod, "SEGMENT_MIN_FILE", 10_000)
    monkeypatch.setattr(fetch_mod, "SPLIT_MIN", 8_000)
    monkeypatch.setattr(fetch_mod, "CHUNK", 4096)
    monkeypatch.setattr(fetch_mod, "SAVE_EVERY", 0.05)


@pytest.fixture
def server():
    s = RangeServer({"big.bin": blob(1), "other.bin": blob(2), **{f"f{i}.bin": blob(10 + i) for i in range(5)}})
    yield s
    s.httpd.shutdown()


def test_segmented_download_is_correct(server, tmp_path):
    f = server.remote("big.bin")
    out = fetch(f, tmp_path, connections=4)
    assert out.read_bytes() == server.files["big.bin"]
    assert not (tmp_path / "big.bin.part").exists() and not (tmp_path / "big.bin.part.json").exists()
    probe = [rng for _, rng, _ in server.log if rng == "bytes=0-0"]
    assert len(probe) == 1
    assert len(server.served_ranges("big.bin")) >= 4  # four segments, maybe more taken over at the end


def test_one_connection_is_the_old_single_stream(server, tmp_path):
    out = fetch(server.remote("big.bin"), tmp_path, connections=1)
    assert out.read_bytes() == server.files["big.bin"]
    assert [rng for _, rng, _ in server.log] == [None]


def test_interrupted_segmented_download_resumes_each_segment(server, tmp_path, monkeypatch):
    data = server.files["big.bin"]
    f = server.remote("big.bin")
    server.cut_after = 20_000  # every connection drops after 20 kB
    monkeypatch.setattr(fetch_mod, "RETRIES", 1)
    monkeypatch.setattr(fetch_mod, "SPLIT_MIN", 50_000)  # no taking over: 4 segments of 75 kB
    with pytest.raises(FetchError, match="re-run to resume"):
        fetch(f, tmp_path, connections=4)
    state = json.loads((tmp_path / "big.bin.part.json").read_text())
    missing = sum(b - a for a, b in state["missing"])
    assert state["size"] == SIZE and 0 < missing < SIZE
    assert pending_bytes(tmp_path) <= missing  # `stew downloads` doesn't count the gaps as done
    # Re-run: only the missing ranges are fetched, then the md5 checks the result.
    server.cut_after = None
    server.log.clear()
    assert fetch(f, tmp_path, connections=4).read_bytes() == data
    fetched = sum(b - a + 1 for a, b in server.served_ranges("big.bin"))
    assert fetched == missing
    assert not any(rng == "bytes=0-0" for _, rng, _ in server.log)  # resumed from its state, no probe
    assert not (tmp_path / "big.bin.part.json").exists()


def test_resumed_segments_can_continue_over_one_connection(server, tmp_path, monkeypatch):
    f = server.remote("big.bin")
    server.cut_after = 20_000
    monkeypatch.setattr(fetch_mod, "RETRIES", 1)
    with pytest.raises(FetchError):
        fetch(f, tmp_path, connections=4)
    server.cut_after = None
    assert fetch(f, tmp_path, connections=1).read_bytes() == server.files["big.bin"]


def test_part_from_a_single_stream_download_is_resumed_not_restarted(server, tmp_path):
    data = server.files["big.bin"]
    have = 123_457
    (tmp_path / "big.bin.part").write_bytes(data[:have])  # left by stew 0.6 (no .part.json)
    assert fetch(server.remote("big.bin"), tmp_path, connections=4).read_bytes() == data
    got = server.served_ranges("big.bin")
    assert min(a for a, _ in got) == have  # nothing before it was fetched again
    assert max(b for _, b in got) == SIZE - 1
    assert len(got) >= 4


def test_falls_back_to_one_connection_without_ranges(tmp_path):
    s = RangeServer({"big.bin": blob(3)}, ranges=False)
    try:
        f = s.remote("big.bin")
        assert fetch(f, tmp_path, connections=8).read_bytes() == s.files["big.bin"]
        assert [code for _, _, code in s.log] == [200, 200]  # the probe, then one plain download
        assert not (tmp_path / "big.bin.part.json").exists()
        # A .part from before is started again, as the server can't send the rest.
        (tmp_path / "x").mkdir()
        (tmp_path / "x/big.bin.part").write_bytes(s.files["big.bin"][:1000])
        assert fetch(f, tmp_path / "x", connections=8).read_bytes() == s.files["big.bin"]
    finally:
        s.httpd.shutdown()


def test_signed_urls_and_small_or_unsized_files_use_one_connection(server, tmp_path):
    f = server.remote("big.bin", resolve=lambda: f"{server.url}/big.bin")  # e.g. HCA: a fresh URL each time
    assert fetch(f, tmp_path / "a", connections=8).read_bytes() == server.files["big.bin"]
    unsized = sources.RemoteFile("big.bin", f"{server.url}/big.bin")
    assert fetch(unsized, tmp_path / "b", connections=8).read_bytes() == server.files["big.bin"]
    assert all(rng is None for _, rng, _ in server.log)


def test_busy_server_slows_down_instead_of_failing(server, tmp_path, capsys):
    server.busy = 3
    f = server.remote("big.bin")
    assert fetch(f, tmp_path, connections=4).read_bytes() == server.files["big.bin"]
    err = capsys.readouterr().err
    assert "server busy (429" in err and "1 connection(s)" in err


def test_busy_server_that_never_recovers_fails_resumably(server, tmp_path, monkeypatch):
    monkeypatch.setattr(fetch_mod, "THROTTLE_RETRIES", 2)
    server.busy = 1000
    with pytest.raises(FetchError, match="stayed busy"):
        fetch(server.remote("big.bin"), tmp_path, connections=2)
    assert (tmp_path / "big.bin.part.json").exists()


def test_retry_after_seconds_or_date():
    import email.utils
    import urllib.error

    def err(value):
        return urllib.error.HTTPError("u", 429, "Too Many Requests", {"Retry-After": value} if value else {}, None)

    assert retry_after(err("7")) == 7
    assert retry_after(err(None), default=3) == 3
    assert retry_after(err("100000")) == fetch_mod.MAX_WAIT
    assert 40 <= retry_after(err(email.utils.formatdate(time.time() + 45, usegmt=True))) <= 47


def test_several_files_at_once_within_the_connection_cap(server, tmp_path, capsys, monkeypatch):
    server.delay = 0.002
    open_now, most = [0], [0]
    lock = threading.Lock()
    real = fetch_mod.request

    class Counted:  # connections stew holds open at once, counted on its side
        def __init__(self, *a, **kw):
            self.r = real(*a, **kw)

        def __enter__(self):
            with lock:
                open_now[0] += 1
                most[0] = max(most[0], open_now[0])
            return self.r.__enter__()

        def __exit__(self, *exc):
            with lock:
                open_now[0] -= 1
            return self.r.__exit__(*exc)

    monkeypatch.setattr(fetch_mod, "request", Counted)
    files = [server.remote(f"f{i}.bin") for i in range(5)] + [sources.RemoteFile("tiny.txt", f"{server.url}/f0.bin")]
    server.files["tiny.txt"] = b"hello"
    files[-1] = server.remote("tiny.txt")
    fetch_many(files, tmp_path, jobs=3, connections=8, max_connections=6)
    for f in files:
        assert (tmp_path / f.name).read_bytes() == server.files[f.name]
    assert 3 <= most[0] <= 6  # files in parallel, never more connections than the cap
    out = capsys.readouterr().out
    assert all(f"[{i}/6]" in out for i in range(1, 7))


def test_a_failing_file_does_not_stop_the_others(server, tmp_path):
    files = [server.remote("f0.bin"), sources.RemoteFile("gone.bin", f"{server.url}/gone.bin", SIZE),
             server.remote("f1.bin")]
    with pytest.raises(FetchError, match="1 of 3 file"):
        fetch_many(files, tmp_path, jobs=2)
    assert (tmp_path / "f0.bin").exists() and (tmp_path / "f1.bin").exists()


def test_download_command_options_and_config(server, tmp_path, capsys):
    catalog = tmp_path / "data/.data-steward/catalog.db"
    Catalog.create(catalog)
    (catalog.parent / "config.toml").write_text("[download]\nconnections = 3\njobs = 2\n")
    url = f"{server.url}/big.bin"
    assert cli.main(["--catalog", str(catalog), "config"]) == 0
    assert "connections 3   jobs 2   max_connections 16" in capsys.readouterr().out
    assert cli.main(["--catalog", str(catalog), "download", url, "--yes", "--connections", "2"]) == 0
    d = Catalog.open(catalog).find("big")
    assert open(d["path"], "rb").read() == server.files["big.bin"]
    assert len(server.served_ranges("big.bin")) >= 2
    with pytest.raises(SystemExit, match="at least 1"):
        cli.main(["--catalog", str(catalog), "download", "zenodo:1", "--connections", "0"])
    (catalog.parent / "config.toml").write_text("[download]\nconections = 3\n")
    cli._SETTINGS.clear()
    assert cli.main(["--catalog", str(catalog), "download", url + "?x", "--yes"]) == 1
    assert "unknown [download] key conections" in capsys.readouterr().err
