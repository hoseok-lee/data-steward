"""Recognises datasets on disk from their layout: instrument output bundles
(Xenium, Visium, CosMx, MERSCOPE, Cell Ranger), analysis objects (h5ad, h5mu,
zarr, rds, loom), raw reads and slide images."""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import util, vocab


@dataclass
class Candidate:
    path: Path
    name: str
    platform: str | None
    format: str
    modalities: set[str] = field(default_factory=set)
    organ: str | None = None
    species: str | None = None
    n_cells: int | None = None
    n_features: int | None = None
    details: dict = field(default_factory=dict)
    note: str | None = None
    kind: str = "dataset"          # or "reference": genomes, indices, annotations
    meta: dict = field(default_factory=dict)  # from a .stew-dataset file

    @classmethod
    def bare(cls, path: Path, platform: str | None, fmt: str, mods=()) -> "Candidate":
        return cls(path, default_name(path), platform, fmt, set(mods))

    @classmethod
    def new(cls, path: Path, platform: str | None, fmt: str, mods=()) -> "Candidate":
        c = cls.bare(path, platform, fmt, mods)
        c.infer_from(str(path))
        return c

    def infer_from(self, text: str) -> None:
        """Fill unset fields from words in a path or title, then add the
        modalities the platform implies."""
        toks = vocab.tokens(text)
        self.organ = self.organ or vocab.infer_organ(toks)
        self.species = self.species or vocab.infer_species(toks)
        self.platform = self.platform or vocab.infer_platform(toks)
        self.modalities |= vocab.infer_modalities(toks)
        self.add_platform_modalities()
        if "cite" in toks or "citeseq" in toks:
            self.modalities.add("rna")
        # Single-cell analysis objects are almost always expression data.
        if self.format in ("anndata", "anndata-zarr", "mudata", "seurat", "rds") and not (
            self.modalities & {"atac", "dna", "methylation"}
        ):
            self.modalities.add("rna")

    def read_contents(self) -> None:
        """Fill in what the files themselves say. Contents beat guesses from the
        path, so they overwrite organ, species and platform when they have them."""
        from . import contents
        c = contents.read(self.path, self.format)
        if c.error:
            self.note = c.error
        self.n_cells = c.n_cells or self.n_cells
        self.n_features = c.n_features or self.n_features
        self.species = c.species or self.species
        self.organ = c.organ or self.organ
        self.platform = c.platform or self.platform
        self.modalities |= c.modalities
        self.details.update(c.details)
        self.add_platform_modalities()

    def add_platform_modalities(self) -> None:
        if self.platform:
            self.modalities |= set(vocab.platform_modalities(self.platform))


def default_name(path: Path) -> str:
    """The folder name, or the file name without extensions."""
    base = path.name
    stem = base if path.is_dir() and not base.endswith(".zarr") else base.split(".")[0]
    return util.slug(stem) or "dataset"


def parent_name(path: Path) -> str | None:
    """`default_name` prefixed with the parent folder, for when the plain name is taken."""
    parent = util.slug(path.parent.name)
    return f"{parent}-{default_name(path)}" if parent else None


@dataclass
class ScanReport:
    found: list[Candidate] = field(default_factory=list)
    unreadable: list[Path] = field(default_factory=list)


def scan(root: Path, max_depth: int = 8) -> ScanReport:
    r = ScanReport()
    if root.is_file():
        c = detect_file(root)
        if c:
            r.found.append(c)
    else:
        _walk(root, 0, max_depth, r)
    return r


def _is_fastq(name: str) -> bool:
    return name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def _walk(d: Path, depth: int, max_depth: int, r: ScanReport) -> None:
    c = detect_dir(d)
    if c:
        r.found.append(c)
        return
    try:
        entries = sorted(os.scandir(d), key=lambda e: e.name)
    except OSError:
        r.unreadable.append(d)
        return
    fastq = bam = False
    for e in entries:
        name = e.name.lower()
        if name.startswith("."):
            continue
        path = Path(e.path)
        if e.is_dir(follow_symlinks=False):
            if depth < max_depth:
                _walk(path, depth + 1, max_depth, r)
        elif e.is_file(follow_symlinks=False):
            if _is_fastq(name):
                fastq = True
            elif name.endswith((".bam", ".cram")):
                bam = True
            elif not name.endswith(".part"):
                c = detect_file(path)
                if c:
                    r.found.append(c)
    # Reads are catalogued per folder rather than per file.
    if fastq or bam:
        fmt = "fastq+bam" if fastq and bam else "fastq" if fastq else "bam"
        r.found.append(Candidate.new(d, None, fmt))


MARKER = ".stew-dataset"

# Files that make a folder of annotations / genome support files a reference.
REFERENCE_EXTS = (".gtf", ".gtf.gz", ".gff", ".gff3", ".gff.gz", ".gff3.gz", ".bed", ".bed.gz", ".chrom.sizes",
                  ".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fai", ".2bit", ".tss", ".dict")
NOTES = (".org", ".md", ".txt", ".sh", ".readme")

# Plain tables, and the words that say a table holds measurements rather than notes.
TABLE_EXTS = (".csv", ".tsv", ".csv.gz", ".tsv.gz", ".txt.gz", ".tab", ".tab.gz")
MEASUREMENT_WORDS = ("count", "counts", "expr", "expression", "matrix", "umi", "tpm", "fpkm", "rpkm", "cpm",
                     "abundance", "quant", "dge")
SIDE_FILES = (".log", ".done", ".md5", ".xlsx", ".pdf", ".json", ".yaml", ".yml")


def _mtx_triplet(files: list[str]) -> bool:
    """matrix.mtx with barcodes and features/genes, as 10x writes them (gzipped
    or not, with or without a sample prefix like GSM123_)."""
    def has(stem):
        return any(f.endswith((f"{stem}.tsv", f"{stem}.tsv.gz", f"{stem}.txt", f"{stem}.txt.gz")) for f in files)
    return any(f.endswith(("matrix.mtx", "matrix.mtx.gz")) for f in files) and has("barcodes") and (
        has("features") or has("genes"))


def _count_tables(d: Path, files: list[str]) -> bool:
    """A folder of plain tables, at least one of which holds measurements
    (counts, expression, TPM, ...), and nothing else but notes and side files."""
    shown = [f for f in files if not f.startswith(".")]
    if not shown or any((d / f).is_dir() for f in shown):
        return False
    tables = [f for f in shown if f.endswith(TABLE_EXTS) or (f.endswith(".txt") and _measures(f))]
    rest = [f for f in shown if f not in tables and not f.endswith(NOTES + SIDE_FILES)]
    return bool(tables) and not rest and any(_measures(f) for f in tables)


def _measures(name: str) -> bool:
    words = set(re.split(r"[^a-z0-9]+", name.lower()))
    return bool(words & set(MEASUREMENT_WORDS)) or "matrix" in name.lower()


def _from_marker(d: Path) -> Candidate | None:
    """A `.stew-dataset` file says "this folder is one dataset". It may hold
    JSON with any of: name, organ, species, modalities, platform, format, kind,
    study, source, description, tags, steward."""
    marker = d / MARKER
    if not marker.is_file():
        return None
    try:
        text = marker.read_text().strip()
        meta = json.loads(text) if text else {}
    except (OSError, ValueError):
        meta = {}
    if not isinstance(meta, dict):
        meta = {}
    c = Candidate.bare(d, vocab.platform(meta["platform"]) if meta.get("platform") else None,
                       meta.get("format") or "directory")
    if meta.get("name"):
        c.name = util.slug(meta["name"]) or c.name
    c.organ = vocab.organ(meta["organ"]) if meta.get("organ") else None
    c.species = vocab.species(meta["species"]) if meta.get("species") else None
    mods = meta.get("modalities") or []
    for m in [mods] if isinstance(mods, str) else mods:
        try:
            c.modalities.add(vocab.modality(m))
        except vocab.VocabError:
            pass
    c.kind = "reference" if meta.get("kind") == "reference" else "dataset"
    c.meta = meta
    c.infer_from(str(d))
    return c


def _reference_genomes(d: Path) -> list[str]:
    try:
        return json.loads((d / "reference.json").read_text()).get("genomes") or []
    except (OSError, ValueError, AttributeError):
        return []


def detect_dir(d: Path) -> Candidate | None:
    """A folder that is a dataset as a whole, e.g. an instrument output bundle."""
    base = d.name.lower()

    c = _from_marker(d)
    if c:
        return c

    # 10x reference packages (cellranger/arc/atac/spaceranger mkref)
    if (d / "reference.json").is_file() and (d / "fasta").is_dir():
        c = Candidate.bare(d, None, "10x-reference")
        c.kind = "reference"
        genomes = _reference_genomes(d)
        c.details["genomes"] = genomes
        species = {vocab.infer_species(vocab.tokens(g)) for g in genomes}
        species.discard(None)
        if len(species) == 1:
            c.species = species.pop()
        elif not genomes:
            c.species = vocab.infer_species(vocab.tokens(d.name))
        return c

    xenium = d / "experiment.xenium"
    if xenium.is_file():
        c = Candidate.bare(d, "xenium", "xenium-bundle", ["rna", "spatial"])
        try:
            meta = json.loads(xenium.read_text())
        except (OSError, ValueError):
            meta = {}
        if meta.get("panel_organism"):
            c.species = vocab.species(meta["panel_organism"])
        if meta.get("panel_tissue_type"):
            c.organ = vocab.infer_organ(vocab.tokens(meta["panel_tissue_type"]))
        c.infer_from(str(d))
        return c

    if base.endswith(".zarr"):
        attrs = ""
        for f in (".zattrs", "zarr.json"):
            try:
                attrs += (d / f).read_text()
            except OSError:
                pass

        def has(p):
            return (d / p).exists()

        if "spatialdata" in attrs or (has("tables") and (has("shapes") or has("points"))):
            return Candidate.new(d, None, "spatialdata", ["spatial"])
        if has("obs") and has("var"):
            return Candidate.new(d, None, "anndata-zarr")
        return Candidate.new(d, None, "zarr")

    # 10x pipelines write their results to <sample>/outs
    for root in (d, d / "outs"):
        if not root.is_dir():
            continue

        def has(p, root=root):
            return (root / p).exists()

        if has("binned_outputs"):
            return Candidate.new(d, "visium-hd", "spaceranger")
        if any(has(f"spatial/{f}") for f in ("tissue_positions.csv", "tissue_positions_list.csv",
                                              "tissue_positions.parquet")):
            return Candidate.new(d, "visium", "spaceranger")
        matrix = has("filtered_feature_bc_matrix.h5") or has("filtered_feature_bc_matrix")
        peaks = has("filtered_peak_bc_matrix.h5") or has("filtered_peak_bc_matrix") or has("peaks.bed")
        if peaks and has("fragments.tsv.gz") and not matrix:
            return Candidate.new(d, "10x-atac", "cellranger-atac")
        if matrix and has("atac_fragments.tsv.gz"):
            return Candidate.new(d, "10x-multiome", "cellranger-arc")
        if matrix or has("per_sample_outs"):
            vdj = has("vdj_t") or has("vdj_b") or has("filtered_contig_annotations.csv")
            return Candidate.new(d, "10x-chromium", "cellranger", ["vdj"] if vdj else [])

    try:
        files = [f.lower() for f in os.listdir(d)]
    except OSError:
        return None
    if any(f.endswith(("exprmat_file.csv", "exprmat_file.csv.gz")) for f in files):
        return Candidate.new(d, "cosmx", "cosmx-flatfiles")
    if "cell_by_gene.csv" in files and ("detected_transcripts.csv" in files or "cell_metadata.csv" in files):
        return Candidate.new(d, "merscope", "merscope-output")
    if any(f in ("fragments.tsv.gz", "atac_fragments.tsv.gz") or f.endswith("_fragments.tsv.gz") for f in files):
        return Candidate.new(d, None, "atac-fragments", ["atac"])
    if _mtx_triplet(files):
        return Candidate.new(d, None, "mtx", ["rna"])
    if _count_tables(d, files):
        return Candidate.new(d, None, "count-tables")
    data = [f for f in files if not f.startswith(".") and not f.endswith(NOTES)]
    subdirs = any((d / f).is_dir() for f in data)
    if data and not subdirs and all(f.endswith(REFERENCE_EXTS) for f in data):
        c = Candidate.bare(d, None, "annotations")
        c.kind = "reference"
        species = vocab.all_species(vocab.tokens(" ".join([str(d)] + data)))
        c.species = species.pop() if len(species) == 1 else None
        c.details["files"] = sorted(data)[:50]
        return c
    return None


_FILE_FORMATS = (
    ((".h5ad",), "anndata", []),
    ((".h5mu",), "mudata", []),
    ((".loom",), "loom", ["rna"]),
    ((".h5seurat",), "seurat", []),
    ((".rds",), "rds", []),
    (("feature_bc_matrix.h5",), "10x-h5", ["rna"]),
    ((".svs", ".ndpi", ".mrxs", ".scn"), "slide-image", ["histology"]),
    ((".ome.tif", ".ome.tiff"), "ome-tiff", ["imaging"]),
)


def detect_file(path: Path) -> Candidate | None:
    """A single file that is a dataset on its own."""
    name = path.name.lower()
    for exts, fmt, mods in _FILE_FORMATS:
        if name.endswith(exts):
            return Candidate.new(path, None, fmt, mods)
    return None


def detect_any(path: Path) -> Candidate:
    """Best guess for an arbitrary path someone asked to register."""
    if path.is_dir():
        found = detect_dir(path)
        if found:
            return found
        try:
            names = [n.lower() for n in os.listdir(path)]
        except OSError:
            names = []
        fastq = any(_is_fastq(n) for n in names)
        bam = any(n.endswith((".bam", ".cram")) for n in names)
        if fastq or bam:
            return Candidate.new(path, None, "fastq+bam" if fastq and bam else "fastq" if fastq else "bam")
        return Candidate.new(path, None, "directory")
    c = detect_file(path)
    if c:
        return c
    suffix = path.name.lower().split(".", 1)[1] if "." in path.name else "file"
    return Candidate.new(path, None, suffix)
