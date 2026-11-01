from __future__ import annotations

import enum
import os
import pwd
import re
from pathlib import Path


def current_user() -> str:
    try:
        return pwd.getpwuid(os.getuid()).pw_name
    except KeyError:
        return os.environ.get("USER", "unknown")


def owner_of(p: Path) -> str | None:
    try:
        uid = p.stat().st_uid
    except OSError:
        return None
    try:
        return pwd.getpwuid(uid).pw_name
    except KeyError:
        return str(uid)


def disk_usage(p: Path) -> tuple[int, int]:
    """Total apparent size and number of regular files under a path."""
    if p.is_file():
        return p.stat().st_size, 1
    size = n = 0
    for root, _dirs, files in os.walk(p):
        for f in files:
            try:
                st = os.lstat(os.path.join(root, f))
            except OSError:
                continue
            if not os.path.islink(os.path.join(root, f)):
                size += st.st_size
                n += 1
    return size, n


class Presence(enum.Enum):
    READABLE = "readable"
    NO_ACCESS = "no access"
    MISSING = "missing"


def readable(p: Path) -> bool:
    if p.is_dir():
        return os.access(p, os.R_OK | os.X_OK)
    return os.access(p, os.R_OK)


def presence(p: Path) -> Presence:
    """Whether a path exists and the current user can read it. A path inside a
    folder we can't enter counts as no access, not missing."""
    try:
        p.stat()
    except PermissionError:
        return Presence.NO_ACCESS
    except OSError:
        return Presence.MISSING
    return Presence.READABLE if readable(p) else Presence.NO_ACCESS


def human_size(n: int | None) -> str:
    if n is None:
        return "-"
    units = ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]
    v, i = float(n), 0
    while v >= 1024 and i < len(units) - 1:
        v /= 1024
        i += 1
    return f"{n} B" if i == 0 else f"{v:.1f} {units[i]}"


def slug(s: str) -> str:
    """Lowercase, dash-separated identifier."""
    return re.sub(r"[^0-9a-z]+", "-", s.lower()).strip("-")


def levenshtein(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a):
        cur = [i + 1]
        for j, cb in enumerate(b):
            cur.append(min(prev[j] + (ca != cb), prev[j + 1] + 1, cur[j] + 1))
        prev = cur
    return prev[-1]


def link_path(p: Path) -> Path:
    """Absolute path with the parent resolved but the last component kept, so a
    symlink is named rather than followed."""
    p = Path(os.path.abspath(p))
    return p.parent.resolve(strict=True) / p.name
