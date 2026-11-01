"""Downloads one file at a time, resuming partial downloads and checking size
and (where the source publishes one) the md5 checksum."""

from __future__ import annotations

import hashlib
import shutil
import sys
import time
import urllib.error
from pathlib import Path

from .sources import RemoteFile, SourceError, request
from .util import human_size

CHUNK = 1 << 20
RETRIES = 4


class FetchError(Exception):
    pass


class Progress:
    def __init__(self, name: str, total: int | None, start: int):
        self.name, self.total, self.done, self.start = name, total, start, start
        self.t0 = self.last = time.monotonic()
        self.tty = sys.stderr.isatty()

    def update(self, n: int, final: bool = False):
        self.done += n
        now = time.monotonic()
        if not final and now - self.last < (0.5 if self.tty else 30):
            return
        self.last = now
        rate = (self.done - self.start) / max(now - self.t0, 1e-6)
        pct = f"{100 * self.done / self.total:5.1f}% " if self.total else ""
        of = f" / {human_size(self.total)}" if self.total else ""
        msg = f"  {self.name}: {pct}{human_size(self.done)}{of}  {human_size(int(rate))}/s"
        if self.tty:
            sys.stderr.write("\r" + msg.ljust(78) + ("\n" if final else ""))
        else:
            sys.stderr.write(msg + "\n")
        sys.stderr.flush()


def digest(p: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(f: RemoteFile, dest_dir: Path) -> Path:
    """Downloads `f` into `dest_dir`. A finished file is left alone; a
    `<name>.part` left by an interrupted run is resumed."""
    dest = (dest_dir / f.name).resolve()
    if not dest.is_relative_to(dest_dir.resolve()):
        raise FetchError(f"refusing to write {f.name!r} outside {dest_dir}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    if dest.exists() and (f.size is None or dest.stat().st_size == f.size):
        print(f"  {f.name}: already here", file=sys.stderr)
        return dest
    for attempt in range(1, RETRIES + 1):
        have = part.stat().st_size if part.exists() else 0
        if f.size is not None and have == f.size:
            break
        headers = {"Range": f"bytes={have}-"} if have and f.url.startswith("http") else {}
        try:
            url = f.resolve() if f.resolve else f.url
            with request(url, headers=headers) as r:
                if have and r.status != 206:
                    have = 0  # the server ignored the range; start over
                total = f.size
                if total is None and r.headers.get("Content-Length", "").isdigit():
                    total = have + int(r.headers["Content-Length"])
                prog = Progress(f.name, total, have)
                with part.open("ab" if have else "wb") as out:
                    while chunk := r.read(CHUNK):
                        out.write(chunk)
                        prog.update(len(chunk))
                prog.update(0, final=True)
                if total is None or part.stat().st_size >= total:
                    break
        except urllib.error.HTTPError as e:
            if e.code == 416 and have:  # nothing left to fetch
                break
            if e.code in (401, 403, 404, 410):
                raise FetchError(f"{f.url} answered {e.code} {e.reason}")
            err = f"{e.code} {e.reason}"
        except (urllib.error.URLError, TimeoutError, OSError, SourceError) as e:
            err = str(getattr(e, "reason", e))
        else:
            err = "connection closed early"
        if attempt == RETRIES:
            raise FetchError(f"{f.name}: {err}; giving up after {RETRIES} tries (re-run to resume)")
        wait = 5 * attempt
        print(f"  {f.name}: {err}; retrying in {wait}s", file=sys.stderr)
        time.sleep(wait)
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


def check_space(dest_dir: Path, needed: int) -> None:
    probe = dest_dir
    while not probe.exists():
        probe = probe.parent
    free = shutil.disk_usage(probe).free
    if needed > free * 0.95:
        raise FetchError(f"not enough space: need {human_size(needed)}, {human_size(free)} free under {probe}")
