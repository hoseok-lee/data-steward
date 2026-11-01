"""Where a pipeline output came from, read from the files the pipeline left
behind. 10x pipelines (cellranger, cellranger-arc, cellranger-atac,
spaceranger) write the command they were run with to `<id>/_cmdline`; the
FASTQ folders it names (directly, or through a --libraries CSV) and the
reference it used become relations in the catalogue."""

from __future__ import annotations

import csv
import os
import shlex
from dataclasses import dataclass, field
from pathlib import Path

INPUT_FLAGS = ("--fastqs",)
LIBRARY_FLAGS = ("--libraries", "--csv")
REFERENCE_FLAGS = ("--reference", "--transcriptome", "--genome")


@dataclass
class Provenance:
    command: str | None = None
    inputs: list[Path] = field(default_factory=list)
    references: list[Path] = field(default_factory=list)

    def empty(self) -> bool:
        return not (self.inputs or self.references)

    def as_details(self) -> dict:
        return {"command": self.command, "inputs": [str(p) for p in self.inputs],
                "references": [str(p) for p in self.references]}


def _flags(tokens: list[str]) -> list[tuple[str, str]]:
    out, i = [], 0
    while i < len(tokens):
        t = tokens[i]
        if t.startswith("--") and "=" in t:
            k, v = t.split("=", 1)
            out.append((k, v))
        elif t.startswith("--") and i + 1 < len(tokens) and not tokens[i + 1].startswith("--"):
            out.append((t, tokens[i + 1]))
            i += 1
        i += 1
    return out


def _libraries(path: Path) -> list[Path]:
    """FASTQ folders listed in a cellranger(-arc) libraries CSV."""
    try:
        with path.open(newline="") as f:
            rows = list(csv.reader(f))
    except OSError:
        return []
    # cellranger multi config: a [libraries] section; arc: a plain CSV with a fastqs column
    start = next((i + 1 for i, r in enumerate(rows) if r and r[0].strip().lower() == "[libraries]"), 0)
    rows = [r for r in rows[start:] if r and not r[0].startswith("[")]
    if not rows:
        return []
    header = [h.strip().lower() for h in rows[0]]
    if "fastqs" not in header:
        return []
    col = header.index("fastqs")
    return [Path(r[col].strip()) for r in rows[1:] if len(r) > col and r[col].strip()]


def read(dataset: Path) -> Provenance:
    """What a pipeline output folder says it was made from; empty if unknown."""
    p = Provenance()
    cmdline = dataset / "_cmdline"
    try:
        text = cmdline.read_text().strip()
    except OSError:
        return p
    try:
        tokens = shlex.split(text)
    except ValueError:
        return p
    if not tokens:
        return p
    p.command = " ".join(Path(t).name for t in tokens[:2])
    # The pipeline ran in the folder that holds <id>/, so relative paths start there.
    base = dataset.parent

    def where(v: str) -> Path:
        return Path(os.path.normpath(base / os.path.expanduser(v)))

    for k, v in _flags(tokens[1:]):
        if k in INPUT_FLAGS:
            p.inputs += [where(x) for x in v.split(",") if x]
        elif k in LIBRARY_FLAGS:
            lib = where(v)
            p.inputs += [x if x.is_absolute() else Path(os.path.normpath(lib.parent / x)) for x in _libraries(lib)]
        elif k in REFERENCE_FLAGS:
            p.references.append(where(v))
    p.inputs = list(dict.fromkeys(p.inputs))
    return p
