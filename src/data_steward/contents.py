"""Reads what is inside a dataset: cell and feature counts, organism, tissue,
assay and modalities, from AnnData/MuData/10x/loom HDF5 files, zarr stores,
Xenium and Space Ranger outputs. Only headers, indices and small categorical
columns are read, never the expression matrix itself."""

from __future__ import annotations

import csv
import json
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from . import vocab

try:
    import h5py
except ImportError:  # stew still works; it just can't look inside HDF5 files
    h5py = None

# obs columns that name these things in CELLxGENE and common pipelines
ORGANISM_COLS = ("organism", "species")
TISSUE_COLS = ("tissue", "organ", "tissue_general", "tissue_type", "anatomical_site")
ASSAY_COLS = ("assay", "technology", "platform")
DISEASE_COLS = ("disease", "condition")
CELLTYPE_COLS = ("cell_type", "celltype", "cell_types", "annotation", "cell_type_annotation")

MOD_NAMES = {"rna": "rna", "gex": "rna", "gene expression": "rna", "atac": "atac", "peaks": "atac",
             "prot": "protein", "protein": "protein", "adt": "protein", "antibody capture": "protein",
             "cite": "protein", "spatial": "spatial", "vdj": "vdj"}

ASSAY_PLATFORMS = (("visium hd", "visium-hd"), ("visium", "visium"), ("xenium", "xenium"), ("cosmx", "cosmx"),
                   ("merfish", "merscope"), ("slide-seq", "slide-seq"), ("stereo-seq", "stereo-seq"),
                   ("multiome", "10x-multiome"), ("smart-seq", "smart-seq"), ("10x", "10x-chromium"))

SMALL = 20000  # don't read categorical columns with more categories than this


@dataclass
class Contents:
    n_cells: int | None = None
    n_features: int | None = None
    species: str | None = None
    organ: str | None = None
    platform: str | None = None
    modalities: set[str] = field(default_factory=set)
    details: dict = field(default_factory=dict)
    error: str | None = None

    def empty(self) -> bool:
        return not (self.n_cells or self.n_features or self.species or self.organ or self.modalities
                    or self.details)


def _s(x) -> str:
    return x.decode("utf-8", "replace") if isinstance(x, bytes) else str(x)


# --- HDF5 helpers ----------------------------------------------------------------------

def _column(group, name: str, limit: int = SMALL):
    """Values of an AnnData obs/var column, as category names (categoricals) or a
    sample of values. Handles the current and the pre-0.8 on-disk layouts."""
    if name not in group:
        cats = group.get("__categories")
        if cats is not None and name in cats:
            return [_s(v) for v in cats[name][:limit]]
        return None
    node = group[name]
    if isinstance(node, h5py.Group):
        if "categories" in node:
            if node["categories"].shape[0] > limit:
                return None
            return [_s(v) for v in node["categories"][:]]
        return None
    if node.ndim != 1:
        return None
    if "__categories" in group and name in group["__categories"]:
        return [_s(v) for v in group["__categories"][name][:limit]]
    return [_s(v) for v in node[:min(limit, node.shape[0])]]


def _index(group) -> tuple[int | None, list[str]]:
    """Length and first names of an obs/var dataframe's index."""
    idx = group.attrs.get("_index")
    idx = _s(idx) if idx is not None else ("_index" if "_index" in group else "index" if "index" in group else None)
    if idx is None or idx not in group:
        return None, []
    ds = group[idx]
    return ds.shape[0], [_s(v) for v in ds[:min(500, ds.shape[0])]]


def _columns(group) -> list[str]:
    order = group.attrs.get("column-order")
    if order is not None:
        return [_s(c) for c in order]
    return [k for k in group.keys() if not k.startswith("_")]


def species_from_genes(names: list[str]) -> str | None:
    """Human gene symbols are upper case (ACTB), mouse ones capitalised (Actb);
    Ensembl ids say so outright."""
    if not names:
        return None
    ens = Counter(re.match(r"^(ENS[A-Z]*?)G\d", n).group(1) for n in names if re.match(r"^ENS[A-Z]*?G\d", n))
    if ens:
        prefix = ens.most_common(1)[0][0]
        return {"ENS": "human", "ENSMUS": "mouse", "ENSRNO": "rat", "ENSDAR": "zebrafish", "ENSMMU": "macaque",
                "ENSSSC": "pig"}.get(prefix)
    words = [n for n in names if re.match(r"^[A-Za-z][A-Za-z0-9-]{1,}$", n) and any(c.isalpha() for c in n[1:])]
    if len(words) < 20:
        return None
    upper = sum(n == n.upper() for n in words) / len(words)
    title = sum(n[0].isupper() and n[1:] == n[1:].lower() for n in words) / len(words)
    if upper > 0.8:
        return "human"
    if title > 0.6:
        return "mouse"
    return None


def _organ_from(values: list[str]) -> tuple[str | None, list[str]]:
    """One organ if all tissue values name the same known organ ("lung",
    "lung parenchyma"), or there is a single value."""
    organs = {vocab.infer_organ(vocab.tokens(v)) for v in values if v}
    organs.discard(None)
    if len(organs) == 1:
        return organs.pop(), sorted(values)[:20]
    if not organs and len(values) == 1:
        return vocab.organ(values[0]), values
    return None, sorted(values)[:20]


def _platform_from_assay(assays: list[str]) -> str | None:
    found = set()
    for a in assays:
        a = a.lower()
        for needle, p in ASSAY_PLATFORMS:
            if needle in a:
                found.add(p)
                break
    return found.pop() if len(found) == 1 else None


def _read_anndata_group(g, c: Contents) -> None:
    """Fills `c` from an AnnData group (a whole .h5ad, or one modality of a .h5mu)."""
    if "obs" in g:
        n, _ = _index(g["obs"])
        c.n_cells = c.n_cells or n
    if "var" in g:
        n, names = _index(g["var"])
        c.n_features = c.n_features or n
        sym_col = next((k for k in ("feature_name", "gene_symbols", "gene_symbol", "gene_name", "symbol")
                        if k in g["var"]), None)
        symbols = (_column(g["var"], sym_col, 500) if sym_col else None) or names
        c.species = c.species or species_from_genes(symbols)
        ftypes = _column(g["var"], "feature_types") or _column(g["var"], "feature_type") or []
        for ft in set(ftypes):
            m = MOD_NAMES.get(ft.lower())
            if m:
                c.modalities.add(m)
    if "X" in g and (c.n_cells is None or c.n_features is None):
        shape = g["X"].attrs.get("shape") if isinstance(g["X"], h5py.Group) else g["X"].shape
        if shape is not None and len(shape) == 2:
            c.n_cells, c.n_features = c.n_cells or int(shape[0]), c.n_features or int(shape[1])
    if "obs" in g:
        obs = g["obs"]
        cols = _columns(obs)
        c.details["obs_columns"] = cols[:40]
        lower = {k.lower(): k for k in cols + list(obs.keys())}

        def values(names):
            for n in names:
                if n in lower:
                    v = _column(obs, lower[n])
                    if v:
                        return sorted(set(v))
            return None

        org = values(ORGANISM_COLS)
        if org:
            sp = {vocab.species(o) for o in org}
            c.details["organism"] = org[:10]
            if len(sp) == 1:
                c.species = sp.pop()
        tis = values(TISSUE_COLS)
        if tis:
            c.organ, c.details["tissues"] = _organ_from(tis)
        assay = values(ASSAY_COLS)
        if assay:
            c.details["assays"] = assay[:10]
            c.platform = c.platform or _platform_from_assay(assay)
        dis = values(DISEASE_COLS)
        if dis:
            c.details["diseases"] = dis[:20]
        ct = values(CELLTYPE_COLS)
        if ct:
            c.details["n_cell_types"] = len(ct)
            c.details["cell_types"] = ct[:30]
    if "obsm" in g:
        keys = list(g["obsm"].keys())
        c.details["obsm"] = keys[:20]
        if any(k.lower() in ("spatial", "x_spatial", "spatial_coords") for k in keys):
            c.modalities.add("spatial")
    if "layers" in g:
        c.details["layers"] = list(g["layers"].keys())[:20]
    if "uns" in g and "spatial" in g["uns"]:
        c.modalities.add("spatial")


def inspect_h5ad(path: Path) -> Contents:
    c = Contents()
    with h5py.File(path, "r") as f:
        _read_anndata_group(f, c)
    if not c.modalities & {"atac", "protein"}:
        c.modalities.add("rna")
    return c


def inspect_h5mu(path: Path) -> Contents:
    c = Contents()
    with h5py.File(path, "r") as f:
        if "obs" in f:
            c.n_cells = _index(f["obs"])[0]
        mods = {}
        for name in f.get("mod", {}):
            sub = Contents()
            _read_anndata_group(f["mod"][name], sub)
            mods[name] = {"n_cells": sub.n_cells, "n_features": sub.n_features}
            c.modalities.add(MOD_NAMES.get(name.lower(), name.lower()))
            c.modalities |= sub.modalities
            c.species = c.species or sub.species
            c.organ = c.organ or sub.organ
            for k, v in sub.details.items():
                c.details.setdefault(k, v)
            c.n_cells = c.n_cells or sub.n_cells
        c.details["modalities"] = mods
        c.n_features = sum(m["n_features"] or 0 for m in mods.values()) or None
    c.modalities &= {t.name for t in vocab.MODALITIES}
    return c


def inspect_10x_h5(path: Path) -> Contents:
    c = Contents()
    with h5py.File(path, "r") as f:
        if "matrix" in f:  # Cell Ranger 3+
            m = f["matrix"]
            c.n_features, c.n_cells = (int(x) for x in m["shape"][:])
            feats = m.get("features")
            if feats is not None:
                types = Counter(_s(t) for t in feats["feature_type"][:]) if "feature_type" in feats else {}
                c.details["feature_types"] = dict(types)
                for t in types:
                    mod = MOD_NAMES.get(t.lower())
                    if mod:
                        c.modalities.add(mod)
                genomes = sorted({_s(g) for g in feats["genome"][:]}) if "genome" in feats else []
                if genomes:
                    c.details["genome"] = genomes
                names = [_s(n) for n in feats["name"][:500]] if "name" in feats else []
                c.species = species_from_genomes(genomes) or species_from_genes(names)
        else:  # Cell Ranger 2: one group per genome
            genomes = list(f.keys())
            c.details["genome"] = genomes
            if genomes and "shape" in f[genomes[0]]:
                c.n_features, c.n_cells = (int(x) for x in f[genomes[0]]["shape"][:])
            c.species = species_from_genomes(genomes)
            c.modalities.add("rna")
    return c


def species_from_genomes(genomes: list[str]) -> str | None:
    sp = {vocab.infer_species(vocab.tokens(g)) for g in genomes}
    sp.discard(None)
    return sp.pop() if len(sp) == 1 else None


def inspect_loom(path: Path) -> Contents:
    c = Contents(modalities={"rna"})
    with h5py.File(path, "r") as f:
        if "matrix" in f:
            c.n_features, c.n_cells = (int(x) for x in f["matrix"].shape)
        genes = f.get("row_attrs/Gene")
        if genes is not None:
            c.species = species_from_genes([_s(g) for g in genes[:500]])
    return c


# --- zarr (JSON metadata only) ---------------------------------------------------------------

def _zarr_shape(node: Path) -> list[int] | None:
    for meta in (node / ".zarray", node / "zarr.json"):
        try:
            return json.loads(meta.read_text()).get("shape")
        except (OSError, ValueError):
            continue
    try:
        return json.loads((node / ".zattrs").read_text()).get("shape")
    except (OSError, ValueError):
        return None


def _zarr_anndata(root: Path, c: Contents) -> None:
    shape = _zarr_shape(root / "X")
    if shape and len(shape) == 2:
        c.n_cells, c.n_features = shape
    if (root / "obsm" / "spatial").exists():
        c.modalities.add("spatial")


def inspect_zarr(path: Path) -> Contents:
    c = Contents()
    if (path / "tables").is_dir():  # SpatialData
        c.modalities.add("spatial")
        tables = [t for t in sorted((path / "tables").iterdir()) if t.is_dir()]
        c.details["tables"] = [t.name for t in tables]
        for t in tables[:1]:
            _zarr_anndata(t, c)
        for elem in ("images", "labels", "points", "shapes"):
            if (path / elem).is_dir():
                c.details[elem] = sorted(p.name for p in (path / elem).iterdir() if p.is_dir())[:20]
    else:
        _zarr_anndata(path, c)
    if c.n_cells:
        c.modalities.add("rna")
    return c


# --- instrument outputs ----------------------------------------------------------------------

def _metrics_csv(path: Path) -> dict:
    try:
        with path.open(newline="") as f:
            rows = list(csv.reader(f))
    except OSError:
        return {}
    if len(rows) >= 2 and len(rows[0]) == len(rows[1]):
        return dict(zip(rows[0], rows[1]))
    return {}


def _number(v) -> int | None:
    try:
        return int(float(str(v).replace(",", "")))
    except ValueError:
        return None


def inspect_10x_dir(path: Path) -> Contents:
    outs = path / "outs" if (path / "outs").is_dir() else path
    c = Contents()
    for h5 in ("filtered_feature_bc_matrix.h5", "binned_outputs/square_008um/filtered_feature_bc_matrix.h5"):
        if (outs / h5).is_file() and h5py:
            try:
                c = inspect_10x_h5(outs / h5)
            except (OSError, KeyError, ValueError):
                pass
            break
    metrics = _metrics_csv(outs / "metrics_summary.csv")
    if metrics:
        c.details["metrics"] = dict(list(metrics.items())[:30])
        for k in ("Estimated Number of Cells", "Number of Spots Under Tissue", "Number of Cells"):
            if k in metrics:
                c.n_cells = _number(metrics[k]) or c.n_cells
                break
        ref = metrics.get("Reference Path") or metrics.get("Genome") or ""
        if ref:
            c.species = c.species or species_from_genomes([Path(ref).name])
    return c


def inspect_xenium(path: Path) -> Contents:
    c = Contents(modalities={"rna", "spatial"})
    try:
        meta = json.loads((path / "experiment.xenium").read_text())
    except (OSError, ValueError):
        meta = {}
    keep = ("run_name", "region_name", "panel_name", "panel_design_id", "panel_organism", "panel_tissue_type",
            "panel_num_targets_predesigned", "panel_num_targets_custom", "instrument_sn",
            "analysis_sw_version", "num_cells", "transcripts_per_cell")
    c.details["xenium"] = {k: meta[k] for k in keep if k in meta}
    c.n_cells = _number(meta.get("num_cells", "")) if meta.get("num_cells") is not None else None
    targets = sum(_number(meta.get(k, 0)) or 0 for k in ("panel_num_targets_predesigned", "panel_num_targets_custom"))
    c.n_features = targets or None
    if meta.get("panel_organism"):
        c.species = vocab.species(meta["panel_organism"])
    if meta.get("panel_tissue_type"):
        c.organ = vocab.infer_organ(vocab.tokens(meta["panel_tissue_type"]))
    if h5py and (path / "cell_feature_matrix.h5").is_file():
        try:
            m = inspect_10x_h5(path / "cell_feature_matrix.h5")
            c.n_cells = c.n_cells or m.n_cells
            c.n_features = c.n_features or m.n_features
            c.details["feature_types"] = m.details.get("feature_types")
        except (OSError, KeyError, ValueError):
            pass
    return c


def inspect_reads(path: Path) -> Contents:
    names = sorted(p.name for p in path.iterdir() if p.is_file())
    fq = [n for n in names if n.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq"))]
    c = Contents()
    if fq:
        reads = Counter(m.group(1) for n in fq if (m := re.search(r"_(R[12]|I[12])[_.]", n)))
        runs = {re.split(r"_(?:R[12]|I[12]|[12])[_.]", n)[0] for n in fq}
        c.details["fastq"] = {"files": len(fq), "reads": dict(reads), "samples_or_runs": len(runs)}
    bams = [n for n in names if n.endswith((".bam", ".cram"))]
    if bams:
        c.details["alignments"] = len(bams)
    return c


def inspect_10x_atac(path: Path) -> Contents:
    """cellranger-atac: cells per genome from summary.json. A barnyard
    reference (GRCh38-and-mm10) still gives one species when nearly all cells
    map to one genome."""
    outs = path / "outs" if (path / "outs").is_dir() else path
    c = Contents(modalities={"atac"})
    try:
        s = json.loads((outs / "summary.json").read_text())
    except (OSError, ValueError):
        s = {}
    per_genome = {k[len("annotated_cells_"):]: v for k, v in s.items()
                  if k.startswith("annotated_cells_") and isinstance(v, (int, float))}
    c.n_cells = s.get("cells_detected") or (sum(per_genome.values()) or None)
    if per_genome:
        c.details["cells_per_genome"] = per_genome
        top, n = max(per_genome.items(), key=lambda kv: kv[1])
        if n >= 0.9 * sum(per_genome.values()):
            c.species = species_from_genomes([top])
    for k in ("cellranger-atac_version", "cellranger-arc_version", "reference_assembly"):
        if k in s:
            c.details[k] = s[k]
    if h5py and (outs / "filtered_peak_bc_matrix.h5").is_file():
        try:
            with h5py.File(outs / "filtered_peak_bc_matrix.h5", "r") as f:
                c.n_features = int(f["matrix"]["shape"][0])
        except (OSError, KeyError, ValueError):
            pass
    return c


def inspect_fragments(path: Path) -> Contents:
    c = Contents(modalities={"atac"})
    for name in ("barcodes.tsv", "singlecell.csv"):
        p = path / name
        if p.is_file() and p.stat().st_size < 50 << 20:
            with p.open() as f:
                c.n_cells = sum(1 for line in f if line.strip()) - (name == "singlecell.csv")
            break
    return c


def inspect_directory(path: Path, limit: int = 20000) -> Contents:
    """A plain folder: what kinds of files it holds."""
    c = Contents()
    kinds: Counter = Counter()
    n = 0
    for _root, dirs, files in os.walk(path):
        dirs[:] = [x for x in dirs if not x.startswith(".")]
        for f in files:
            if f.startswith("."):
                continue
            low = f.lower()
            ext = next((e for e in (".fastq.gz", ".fq.gz", ".tar.gz", ".tsv.gz", ".csv.gz", ".bed.gz", ".gtf.gz")
                        if low.endswith(e)), Path(low).suffix or "(none)")
            kinds[ext] += 1
            n += 1
            if n >= limit:
                break
        if n >= limit:
            c.details["note"] = f"stopped counting at {limit} files"
            break
    if kinds:
        c.details["file_types"] = dict(kinds.most_common(15))
    return c


def inspect_mtx(path: Path) -> Contents:
    """Features x cells from the MatrixMarket header; species from gene names."""
    import gzip
    c = Contents(modalities={"rna"})
    names = sorted(p.name for p in path.iterdir())

    def opener(name):
        p = path / name
        return gzip.open(p, "rt") if name.endswith(".gz") else p.open()

    mtx = next((n for n in names if n.endswith(("matrix.mtx", "matrix.mtx.gz"))), None)
    if mtx:
        with opener(mtx) as f:
            for line in f:
                if not line.startswith("%"):
                    rows, cols, nnz = (int(x) for x in line.split()[:3])
                    c.n_features, c.n_cells = rows, cols
                    c.details["nonzero"] = nnz
                    break
    feat = next((n for n in names if re.search(r"(features|genes)\.(tsv|txt)(\.gz)?$", n)), None)
    if feat:
        with opener(feat) as f:
            rows = [line.rstrip("\n").split("\t") for _, line in zip(range(500), f)]
        c.species = (species_from_genes([r[0] for r in rows])  # Ensembl ids, when present
                     or species_from_genes([r[1] for r in rows if len(r) > 1]))
        types = {r[2] for r in rows if len(r) > 2}
        c.modalities |= {MOD_NAMES[t.lower()] for t in types if t.lower() in MOD_NAMES}
    return c


def inspect_tables(path: Path) -> Contents:
    """Size of each table (columns from the header, rows if small)."""
    import gzip
    c = Contents()
    tables = {}
    for p in sorted(path.iterdir()):
        if not p.is_file() or p.name.startswith("."):
            continue
        try:
            with (gzip.open(p, "rt") if p.name.endswith(".gz") else p.open()) as f:
                header = f.readline()
        except (OSError, UnicodeDecodeError, EOFError):
            continue
        sep = "\t" if header.count("\t") >= header.count(",") else ","
        tables[p.name] = {"columns": len(header.split(sep))}
    if tables:
        c.details["tables"] = dict(list(tables.items())[:20])
    return c


INSPECTORS = {
    "anndata": inspect_h5ad, "mudata": inspect_h5mu, "10x-h5": inspect_10x_h5, "loom": inspect_loom,
    "anndata-zarr": inspect_zarr, "spatialdata": inspect_zarr, "zarr": inspect_zarr,
    "cellranger": inspect_10x_dir, "cellranger-arc": inspect_10x_dir, "spaceranger": inspect_10x_dir,
    "xenium-bundle": inspect_xenium, "fastq": inspect_reads, "bam": inspect_reads, "fastq+bam": inspect_reads,
    "cellranger-atac": inspect_10x_atac, "atac-fragments": inspect_fragments, "directory": inspect_directory,
    "mtx": inspect_mtx, "count-tables": inspect_tables,
}
NEEDS_H5PY = {"anndata", "mudata", "10x-h5", "loom"}


def read(path: Path, fmt: str | None) -> Contents:
    """What can be learnt from inside `path`; never raises."""
    fn = INSPECTORS.get(fmt or "")
    if fn is None:
        return Contents()
    if fmt in NEEDS_H5PY and h5py is None:
        return Contents(error="h5py is not available, so HDF5 files can't be read")
    try:
        return fn(path)
    except (OSError, KeyError, ValueError, TypeError, IndexError) as e:
        return Contents(error=f"could not read {path.name}: {e}")
