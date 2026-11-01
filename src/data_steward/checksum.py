"""File checksums: a manifest of every file in a dataset with its size, mtime
and hash, and a fingerprint of the whole dataset built from it.

Full hashes read every byte (BLAKE2b). Quick hashes read only the first and
last MiB plus the size, which is enough to spot copies and most changes in a
fraction of the time; fingerprints say which kind they are ("full:" or "quick:")
and only like are compared with like."""

from __future__ import annotations

import hashlib
import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

CHUNK = 4 << 20
EDGE = 1 << 20


def file_hash(path: Path, quick: bool) -> str:
    h = hashlib.blake2b(digest_size=20)
    size = path.stat().st_size
    with path.open("rb") as f:
        if quick:
            h.update(size.to_bytes(8, "little"))
            h.update(f.read(EDGE))
            if size > 2 * EDGE:
                f.seek(size - EDGE)
                h.update(f.read(EDGE))
        else:
            while chunk := f.read(CHUNK):
                h.update(chunk)
    return ("quick:" if quick else "full:") + h.hexdigest()


def list_files(root: Path) -> list[tuple[str, Path]]:
    """(relative path, path) for every regular file; a single-file dataset has
    relative path "" so copies under other names still match."""
    if root.is_file():
        return [("", root)]
    out = []
    for d, dirs, files in os.walk(root):
        dirs.sort()
        for f in sorted(files):
            p = Path(d) / f
            if not p.is_symlink() and not f.endswith(".part"):
                out.append((str(p.relative_to(root)), p))
    return out


def manifest(root: Path, quick: bool, jobs: int = 4, progress=None) -> list[tuple[str, int, float, str]]:
    """[(relpath, size, mtime, hash)] for every file under `root`."""
    files = list_files(root)

    def one(item):
        rel, p = item
        st = p.stat()
        h = file_hash(p, quick)
        if progress:
            progress(st.st_size)
        return rel, st.st_size, st.st_mtime, h

    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        return list(pool.map(one, files))


def fingerprint(entries) -> str:
    """One hash for the whole dataset, from its file names, sizes and hashes."""
    quick = any(e[3].startswith("quick:") for e in entries)
    h = hashlib.blake2b(digest_size=16)
    for rel, size, _mtime, fh in sorted(entries):
        h.update(f"{rel}\0{size}\0{fh}\n".encode())
    return ("quick:" if quick else "full:") + h.hexdigest()


@dataclass
class Diff:
    changed: list[str]
    missing: list[str]
    added: list[str]

    @property
    def ok(self) -> bool:
        return not (self.changed or self.missing or self.added)

    def summary(self) -> str:
        parts = []
        for label, items in (("changed", self.changed), ("missing", self.missing), ("new", self.added)):
            if items:
                shown = ", ".join(i or "(the file)" for i in items[:5]) + (" ..." if len(items) > 5 else "")
                parts.append(f"{len(items)} {label}: {shown}")
        return "; ".join(parts) or "unchanged"


def compare(old: dict[str, dict], new: list[tuple]) -> Diff:
    now = {e[0]: e for e in new}
    changed = sorted(r for r in now.keys() & old.keys() if (now[r][1], now[r][3]) != (old[r]["size"], old[r]["hash"]))
    return Diff(changed, sorted(old.keys() - now.keys()), sorted(now.keys() - old.keys()))
