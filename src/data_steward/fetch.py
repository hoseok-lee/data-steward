"""Downloads files, resuming partial downloads and checking size and (where the
source publishes one) the md5 or sha256 checksum.

A large file on a server that answers ranged requests is fetched over several
connections at once: its missing byte ranges are split into segments, each
written at its own offset in `<name>.part`. What is still missing is kept in
`<name>.part.json`, so an interrupted download resumes every segment where it
stopped. A `.part` without that file (left by a single-connection download) is
resumed too: its bytes count as done and only the rest is split. Connections
that finish early take over half of the largest segment still running.

Several files of a dataset can be fetched at once (`fetch_many`), with a cap on
the connections in use altogether."""

from __future__ import annotations

import email.utils
import hashlib
import http.client
import json
import os
import shutil
import sys
import threading
import time
import urllib.error
from dataclasses import dataclass
from pathlib import Path

from .sources import RemoteFile, SourceError, request
from .util import human_size

CHUNK = 1 << 20
RETRIES = 4                  # failures in a row of one connection (or segment) before giving up
THROTTLE_RETRIES = 12        # "busy" answers (429/503) in a row before giving up
SEGMENT_MIN_FILE = 64 << 20  # smaller files are fetched over one connection
SPLIT_MIN = 8 << 20          # segments are never split into pieces smaller than this
SAVE_EVERY = 5.0             # seconds between saves of the segment state
BUSY_WAIT = 30               # seconds to wait after a 429/503 without Retry-After
MAX_WAIT = 600               # longest Retry-After honoured

CONNECTIONS = 8              # per file
JOBS = 3                     # files at once
MAX_CONNECTIONS = 16         # altogether

NETWORK_ERRORS = (urllib.error.URLError, TimeoutError, OSError, SourceError, http.client.HTTPException)


class FetchError(Exception):
    pass


class _Retire(Exception):
    """A connection closes because the server asked us to slow down."""


class Stopped(Exception):
    """Another part of the download was interrupted; this file's state was saved."""


# --- progress ----------------------------------------------------------------------------

class Track:
    def __init__(self, board: Board, name: str, total: int | None, have: int, connections: int):
        self.board, self.name, self.total, self.connections = board, name, total, connections
        self.done = self.start = have
        self.t0 = time.monotonic()

    def reset(self, have: int):
        with self.board.lock:
            self.done = self.start = have
            self.t0 = time.monotonic()

    def line(self) -> str:
        rate = (self.done - self.start) / max(time.monotonic() - self.t0, 1e-6)
        pct = f"{100 * self.done / self.total:5.1f}% " if self.total else ""
        of = f" / {human_size(self.total)}" if self.total else ""
        conns = f"  ({self.connections} connections)" if self.connections > 1 else ""
        return f"  {self.name}: {pct}{human_size(self.done)}{of}  {human_size(int(rate))}/s{conns}"

    def rate(self) -> float:
        return (self.done - self.start) / max(time.monotonic() - self.t0, 1e-6)


class Board:
    """Progress of the files being downloaded. On a terminal: one line, redrawn
    twice a second (all running files on it). In a log: a line per running file
    every 30 seconds. A last line for each file when it is finished."""

    def __init__(self):
        self.lock = threading.RLock()
        self.tracks: list[Track] = []
        self.tty = sys.stderr.isatty()
        self.last = time.monotonic()
        self.drawn = False

    def start(self, name: str, total: int | None, have: int, connections: int = 1) -> Track:
        t = Track(self, name, total, have, connections)
        with self.lock:
            self.tracks.append(t)
        return t

    def update(self, t: Track, n: int):
        with self.lock:
            t.done += n
            now = time.monotonic()
            if now - self.last < (0.5 if self.tty else 30):
                return
            self.last = now
            if self.tty:
                self._draw()
            else:
                for x in self.tracks:
                    sys.stderr.write(x.line() + "\n")
                sys.stderr.flush()

    def _draw(self):
        width = max(40, shutil.get_terminal_size((100, 20)).columns - 1)
        if len(self.tracks) == 1:
            msg = self.tracks[0].line()
        else:
            parts = [f"{x.name.rsplit('/', 1)[-1][:24]} {100 * x.done / x.total:.0f}%" if x.total
                     else f"{x.name.rsplit('/', 1)[-1][:24]} {human_size(x.done)}" for x in self.tracks]
            rate = sum(x.rate() for x in self.tracks)
            msg = f"  {len(self.tracks)} files: " + " | ".join(parts) + f"  {human_size(int(rate))}/s"
        sys.stderr.write("\r" + msg[:width].ljust(width))
        sys.stderr.flush()
        self.drawn = True

    def _clear(self):
        if self.tty and self.drawn:
            width = max(40, shutil.get_terminal_size((100, 20)).columns - 1)
            sys.stderr.write("\r" + " " * width + "\r")
            self.drawn = False

    def say(self, msg: str, out=None):
        """Prints a message without mangling the progress line."""
        with self.lock:
            self._clear()
            stream = out or sys.stderr
            stream.write(msg + "\n")
            stream.flush()

    def finish(self, t: Track):
        with self.lock:
            if t in self.tracks:
                self.tracks.remove(t)
            self.say(t.line())


# --- one file ----------------------------------------------------------------------------

def digest(p: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def retry_after(e: urllib.error.HTTPError, default: int | None = None) -> int:
    """Seconds a 429/503 answer asks us to wait (Retry-After: seconds or a date)."""
    value = (e.headers.get("Retry-After") if e.headers else None) or ""
    value = value.strip()
    wait = BUSY_WAIT if default is None else default
    if value.isdigit():
        wait = int(value)
    elif value:
        try:
            wait = int(email.utils.parsedate_to_datetime(value).timestamp() - time.time()) + 1
        except (TypeError, ValueError, OverflowError):
            pass
    return max(0, min(wait, MAX_WAIT))


def _http_error(e: urllib.error.HTTPError, url: str) -> None:
    if e.code in (401, 403, 404, 410):
        raise FetchError(f"{url} answered {e.code} {e.reason}")


def state_file(part: Path) -> Path:
    return part.with_name(part.name + ".json")


def fetch(f: RemoteFile, dest_dir: Path, connections: int = CONNECTIONS, board: Board | None = None,
          stop: threading.Event | None = None) -> Path:
    """Downloads `f` into `dest_dir`, over up to `connections` connections. A
    finished file is left alone; a `<name>.part` left by an interrupted run is
    resumed."""
    dest = (dest_dir / f.name).resolve()
    if not dest.is_relative_to(dest_dir.resolve()):
        raise FetchError(f"refusing to write {f.name!r} outside {dest_dir}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    state = state_file(part)
    board = board or Board()
    stop = stop or threading.Event()
    if dest.exists() and (f.size is None or dest.stat().st_size == f.size):
        board.say(f"  {f.name}: already here")
        return dest
    ranges = _load_state(f, part, state, board)
    if ranges is None and connections > 1 and _splittable(f) and _accepts_ranges(f.url, f.size):
        have = part.stat().st_size if part.exists() else 0
        if have > f.size:
            have = 0
        ranges = [[have, f.size]] if have < f.size else []
        if have:
            board.say(f"  {f.name}: resuming at {human_size(have)}; fetching the rest over {connections} connections")
        with part.open("r+b" if have else "wb") as out:
            out.truncate(have)
        _save_state(state, f.size, ranges)
    if ranges is not None:
        Segmented(f, part, state, ranges, connections, board, stop).run()
        state.unlink(missing_ok=True)
    else:
        _fetch_single(f, part, board, stop)
    size = part.stat().st_size
    if f.size is not None and size != f.size:
        raise FetchError(f"{f.name}: got {size} bytes, expected {f.size} (re-run to resume)")
    for algorithm, want in (("md5", f.md5), ("sha256", f.sha256)):
        if want:
            got = digest(part, algorithm)
            if got != want.lower():
                part.unlink()
                raise FetchError(f"{f.name}: {algorithm} {got} does not match the published {want}; deleted it")
    part.rename(dest)
    return dest


def _splittable(f: RemoteFile) -> bool:
    # Signed URLs handed out per request (HCA) are fetched over one connection.
    return (f.size is not None and f.size >= SEGMENT_MIN_FILE and f.resolve is None
            and f.url.lower().startswith(("http://", "https://")))


def _accepts_ranges(url: str, size: int) -> bool:
    """Whether the server answers a ranged request with that range (206) of a file this size."""
    try:
        with request(url, headers={"Range": "bytes=0-0"}) as r:
            ok = r.status == 206 and r.headers.get("Content-Range", "").strip().endswith(f"/{size}")
            r.read(1)
            return ok
    except (urllib.error.HTTPError, *NETWORK_ERRORS):
        return False


def _load_state(f: RemoteFile, part: Path, state: Path, board: Board) -> list[list[int]] | None:
    """The missing ranges recorded by an interrupted segmented download, or None."""
    if not state.exists():
        return None
    try:
        d = json.loads(state.read_text())
        ranges = [[int(a), int(b)] for a, b in d["missing"]]
        size = int(d["size"])
    except (OSError, ValueError, KeyError, TypeError):
        d = None
    if d is None or not part.exists() or (f.size is not None and size != f.size):
        why = "its size changed" if d and f.size is not None and size != f.size else "the partial file is unusable"
        board.say(f"  {f.name}: {why}; starting it again")
        state.unlink(missing_ok=True)
        part.unlink(missing_ok=True)
        return None
    if ranges:
        left = sum(b - a for a, b in ranges)
        board.say(f"  {f.name}: resuming, {human_size(size - left)} of {human_size(size)} already here")
    return ranges


def _save_state(state: Path, size: int, ranges: list[list[int]]) -> None:
    tmp = state.with_name(state.name + ".tmp")
    tmp.write_text(json.dumps({"size": size, "missing": ranges}) + "\n")
    os.replace(tmp, state)


def pending_bytes(folder: Path) -> int:
    """Bytes of the `.part` files under `folder` that are counted in their
    apparent size but not downloaded yet (segments are written ahead)."""
    pending = 0
    try:
        states = list(folder.rglob("*.part.json")) if folder.is_dir() else []
    except OSError:
        return 0
    for st in states:
        try:
            d = json.loads(st.read_text())
            part = st.with_name(st.name.removesuffix(".json"))
            size = part.stat().st_size
            pending += sum(max(0, min(b, size) - a) for a, b in d["missing"] if a < size)
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return pending


def _fetch_single(f: RemoteFile, part: Path, board: Board, stop: threading.Event) -> None:
    """One connection, appending to the .part (resumed with a Range request)."""
    track = None
    failures = busy = 0
    try:
        while not stop.is_set():
            have = part.stat().st_size if part.exists() else 0
            if f.size is not None and have == f.size:
                break
            headers = {"Range": f"bytes={have}-"} if have and f.url.startswith("http") else {}
            wait = None
            try:
                url = f.resolve() if f.resolve else f.url
                with request(url, headers=headers) as r:
                    if have and r.status != 206:
                        have = 0  # the server ignored the range; start over
                    total = f.size
                    if total is None and r.headers.get("Content-Length", "").isdigit():
                        total = have + int(r.headers["Content-Length"])
                    if track is None:
                        track = board.start(f.name, total, have)
                    else:
                        track.total = total
                        track.reset(have)
                    with part.open("ab" if have else "wb") as out:
                        while chunk := r.read(CHUNK):
                            out.write(chunk)
                            board.update(track, len(chunk))
                            failures = busy = 0
                            if stop.is_set():
                                raise Stopped()
                    if total is None or part.stat().st_size >= total:
                        break
            except urllib.error.HTTPError as e:
                if e.code == 416 and have:  # nothing left to fetch
                    break
                _http_error(e, f.url)
                err = f"{e.code} {e.reason}"
                if e.code in (429, 503):
                    busy += 1
                    wait = retry_after(e)
                    if busy > THROTTLE_RETRIES:
                        raise FetchError(f"{f.name}: the server stayed busy ({err}); giving up (re-run to resume)")
            except NETWORK_ERRORS as e:
                err = str(getattr(e, "reason", e)) or type(e).__name__
            else:
                err = "connection closed early"
            if wait is None:
                failures += 1
                if failures >= RETRIES:
                    raise FetchError(f"{f.name}: {err}; giving up after {RETRIES} tries (re-run to resume)")
                wait = 5 * failures
            board.say(f"  {f.name}: {err}; retrying in {wait}s")
            if stop.wait(wait):
                break
    finally:
        if track is not None:
            board.finish(track)
    if stop.is_set():
        raise Stopped()


@dataclass
class Segment:
    pos: int
    end: int
    failures: int = 0
    busy: int = 0

    @property
    def left(self) -> int:
        return self.end - self.pos


class Segmented:
    """One file over several connections; see the module docstring."""

    def __init__(self, f: RemoteFile, part: Path, state: Path, ranges: list[list[int]], connections: int,
                 board: Board, stop: threading.Event):
        self.f, self.part, self.state, self.board, self.stop = f, part, state, board, stop
        self.size = f.size if f.size is not None else max((b for _, b in ranges), default=0)
        self.lock = threading.Lock()
        self.todo = [Segment(a, b) for a, b in ranges if b > a]
        self.active: list[Segment] = []
        while len(self.todo) < connections:  # one segment per connection to begin with
            big = max(self.todo, key=lambda s: s.left, default=None)
            if big is None or big.left < 2 * SPLIT_MIN:
                break
            self.todo.append(self._split(big))
        self.todo.sort(key=lambda s: s.pos)
        self.limit = max(1, connections)
        self.running = 0
        self.not_before = 0.0
        self.error: Exception | None = None
        have = self.size - sum(s.left for s in self.todo)
        self.track = board.start(f.name, self.size, have, min(self.limit, max(1, len(self.todo))))

    @staticmethod
    def _split(s: Segment) -> Segment:
        mid = s.pos + s.left // 2
        new = Segment(mid, s.end)
        s.end = mid
        return new

    def missing(self) -> list[list[int]]:
        with self.lock:
            return sorted([s.pos, s.end] for s in self.todo + self.active if s.left > 0)

    def save(self):
        _save_state(self.state, self.size, self.missing())

    def take(self) -> Segment | None:
        """The next segment for a connection: a waiting one, or the second half
        of the largest one in progress."""
        with self.lock:
            if self.running > self.limit:
                return None
            if self.todo:
                s = self.todo.pop(0)
            else:
                big = max(self.active, key=lambda s: s.left, default=None)
                if big is None or big.left < 2 * SPLIT_MIN:
                    return None
                s = self._split(big)
            self.active.append(s)
            return s

    def put_back(self, s: Segment):
        with self.lock:
            if s in self.active:
                self.active.remove(s)
            if s.left > 0:
                self.todo.insert(0, s)

    def done_with(self, s: Segment):
        with self.lock:
            if s in self.active:
                self.active.remove(s)

    def fail(self, e: Exception):
        with self.lock:
            if self.error is None:
                self.error = e

    def run(self):
        if not self.todo:
            self.board.finish(self.track)
            return
        self.part.touch()
        workers = [threading.Thread(target=self.worker, daemon=True, name=f"stew-{i}")
                   for i in range(min(self.limit, len(self.todo) * 64))]
        with self.lock:
            self.running = len(workers)
        for w in workers:
            w.start()
        saved = time.monotonic()
        try:
            while any(w.is_alive() for w in workers):
                for w in workers:
                    w.join(0.25)
                if time.monotonic() - saved >= SAVE_EVERY:
                    self.save()
                    saved = time.monotonic()
        except BaseException:
            self.stop.set()
            raise
        finally:
            self.save()
            self.board.finish(self.track)
        if self.error:
            raise self.error
        if self.stop.is_set():
            raise Stopped()
        if self.missing():
            raise FetchError(f"{self.f.name}: some parts are still missing (re-run to resume)")

    def worker(self):
        retired = False
        try:
            with self.part.open("r+b") as out:
                while not self.stop.is_set() and self.error is None:
                    s = self.take()
                    if s is None:
                        return
                    if not self.fetch_segment(s, out):
                        return
        except _Retire:
            retired = True
        except Exception as e:  # noqa: BLE001 - reported by run()
            self.fail(e if isinstance(e, FetchError) else FetchError(f"{self.f.name}: {e}"))
        finally:
            if not retired:
                with self.lock:
                    self.running -= 1

    def fetch_segment(self, s: Segment, out) -> bool:
        """Fetches `s` until it is done (True) or this connection should stop (False)."""
        name = self.f.name
        while s.left > 0 and not self.stop.is_set() and self.error is None:
            pause = self.not_before - time.monotonic()
            if pause > 0 and self.stop.wait(pause):
                break
            wait = None
            try:
                with request(self.f.url, headers={"Range": f"bytes={s.pos}-{s.end - 1}"}) as r:
                    if r.status != 206:
                        raise FetchError(f"{name}: the server stopped answering ranged requests; delete "
                                         f"{self.part} and {self.state} to download it from the start")
                    while not self.stop.is_set():
                        chunk = r.read(CHUNK)
                        if not chunk:
                            break
                        with self.lock:  # write under the lock: the segment may be split meanwhile
                            n = min(len(chunk), s.left)
                            if n > 0:
                                out.seek(s.pos)
                                out.write(chunk[:n])
                                out.flush()
                                s.pos += n
                            s.failures = s.busy = 0
                        self.board.update(self.track, n)
                        if s.left <= 0:
                            break
                if s.left <= 0 or self.stop.is_set():
                    break
                err = "connection closed early"
            except urllib.error.HTTPError as e:
                _http_error(e, self.f.url)
                if e.code == 416:
                    raise FetchError(f"{name}: the server says bytes {s.pos}-{s.end - 1} are past the end of "
                                     "the file; it may have changed")
                err = f"{e.code} {e.reason}"
                if e.code in (429, 503):
                    s.busy += 1
                    if s.busy > THROTTLE_RETRIES:
                        raise FetchError(f"{name}: the server stayed busy ({err}); giving up (re-run to resume)")
                    wait = retry_after(e)
                    with self.lock:
                        self.not_before = max(self.not_before, time.monotonic() + wait)
                        self.limit = max(1, self.limit - 1)
                        self.track.connections = self.limit
                        retire = self.running > self.limit
                        if retire:  # this connection closes; decided here so the last one never does
                            self.running -= 1
                    self.board.say(f"  {name}: server busy ({err}); waiting {wait}s, then using "
                                   f"{self.limit} connection(s)")
                    if retire:
                        self.put_back(s)
                        raise _Retire()
                    continue
            except NETWORK_ERRORS as e:
                err = str(getattr(e, "reason", e)) or type(e).__name__
            s.failures += 1
            if s.failures >= RETRIES:
                raise FetchError(f"{name}: {err}; giving up after {RETRIES} tries (re-run to resume)")
            wait = 5 * s.failures
            self.board.say(f"  {name}: {err}; retrying bytes {s.pos}- in {wait}s")
            if self.stop.wait(wait):
                break
        if s.left > 0:
            self.put_back(s)
            return False
        self.done_with(s)
        return True


# --- several files -----------------------------------------------------------------------

def fetch_many(files: list[RemoteFile], dest_dir: Path, jobs: int = JOBS, connections: int = CONNECTIONS,
               max_connections: int = MAX_CONNECTIONS) -> None:
    """Downloads `files` into `dest_dir`, up to `jobs` at a time, using at most
    `max_connections` connections altogether (but at least `connections` for a
    single file). Files that fail don't stop the others; the errors are raised
    at the end."""
    board = Board()
    n = len(files)
    cap = max(max_connections, connections) if jobs <= 1 or n == 1 else max_connections
    if jobs <= 1 or n == 1:
        for i, f in enumerate(files, 1):
            if n > 1:
                board.say(f"[{i}/{n}] {f.name}", sys.stdout)
            fetch(f, dest_dir, min(connections, cap), board)
        return
    stop = threading.Event()
    lock = threading.Lock()
    queue = list(enumerate(files, 1))
    errors: list[str] = []
    unfinished = [n]

    def runner():
        while not stop.is_set():
            with lock:
                if not queue:
                    return
                i, f = queue.pop(0)
                conns = max(1, min(connections, cap // min(jobs, unfinished[0])))
            board.say(f"[{i}/{n}] {f.name}", sys.stdout)
            try:
                fetch(f, dest_dir, conns, board, stop)
            except Stopped:
                return
            except FetchError as e:
                with lock:
                    errors.append(str(e))
                board.say(f"  {e}")
            except Exception as e:  # noqa: BLE001 - reported with the others
                with lock:
                    errors.append(f"{f.name}: {e}")
                board.say(f"  {f.name}: {e}")
            with lock:
                unfinished[0] -= 1

    threads = [threading.Thread(target=runner, daemon=True, name=f"stew-file-{i}") for i in range(min(jobs, n))]
    for t in threads:
        t.start()
    try:
        while any(t.is_alive() for t in threads):
            for t in threads:
                t.join(0.25)
    except BaseException:
        stop.set()
        deadline = time.monotonic() + 10  # let each file save where it got to
        for t in threads:
            t.join(max(0, deadline - time.monotonic()))
        raise
    if errors:
        more = f" (and {len(errors) - 1} more)" if len(errors) > 1 else ""
        raise FetchError(f"{len(errors)} of {n} file(s) failed: {errors[0]}{more}")


def check_space(dest_dir: Path, needed: int) -> None:
    probe = dest_dir
    while not probe.exists():
        probe = probe.parent
    free = shutil.disk_usage(probe).free
    if needed > free * 0.95:
        raise FetchError(f"not enough space: need {human_size(needed)}, {human_size(free)} free under {probe}")
