#!/usr/bin/env python3
"""Builds a fake data tree that mimics real instrument outputs and analysis
objects, for testing data-steward.

Analysis objects (h5ad, h5mu, loom, 10x h5) are small but real HDF5 files with
CELLxGENE-style annotations, so `stew inspect` has something to read. Bulky
things (raw reads, images, transcripts) are sparse files: their sizes look
realistic but they take almost no disk space.

    python scripts/make-dummy-data.py /tmp/dummy-data
"""

import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np

RNG = np.random.default_rng(0)
STR = h5py.string_dtype()

HUMAN_GENES = ["ACTB", "GAPDH", "EPCAM", "KRT19", "KRT8", "PTPRC", "CD3E", "CD4", "CD8A", "MS4A1", "CD19", "LYZ",
               "CD68", "COL1A1", "COL1A2", "DCN", "ACTA2", "PECAM1", "VWF", "SFTPC", "SFTPB", "AGER", "ALB", "APOA1",
               "INS", "GCG", "PRSS1", "CPA1", "KRT5", "TP63", "MKI67", "TOP2A", "VIM", "FN1", "SPARC", "POSTN",
               "FAP", "PDGFRB", "RGS5", "MYH11", "NKG7", "GNLY", "FOXP3", "IL7R", "CCR7", "HLA-DRA", "CD74", "C1QA"]


def genes(n, species="human", ensembl=False):
    if ensembl:
        prefix = "ENSG" if species == "human" else "ENSMUSG"
        return [f"{prefix}{i:011d}" for i in range(n)]
    base = HUMAN_GENES + [f"GENE{i}" for i in range(n)]
    names = base[:n]
    return [g.capitalize() if species == "mouse" else g for g in names]


def sparse_file(p: Path, size: str = "", same_as: Path | None = None) -> None:
    """A file of `size` that takes almost no disk. It starts with a header
    derived from its path (or from `same_as`, to make an exact copy), so only
    deliberate copies have identical contents."""
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "wb") as f:
        if size:
            units = {"K": 1 << 10, "M": 1 << 20, "G": 1 << 30}
            n = int(float(size[:-1]) * units[size[-1]])
            f.write((hashlib.sha1(str(same_as or p).encode()).digest() * 52)[:min(n, 1024)])
            f.truncate(n)


def write_strings(g, name, values):
    g.create_dataset(name, data=np.array(values, dtype=object), dtype=STR)


def write_df(g, name, index, columns):
    """An AnnData dataframe; columns given as lists become categoricals."""
    df = g.create_group(name)
    df.attrs.update({"encoding-type": "dataframe", "encoding-version": "0.2.0", "_index": "_index"})
    df.attrs.create("column-order", np.array(list(columns), dtype=STR))
    write_strings(df, "_index", index)
    for col, values in columns.items():
        values = list(values)
        if isinstance(values[0], str):
            cats = sorted(set(values))
            c = df.create_group(col)
            c.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
            write_strings(c, "categories", cats)
            code = {v: i for i, v in enumerate(cats)}
            c.create_dataset("codes", data=np.array([code[v] for v in values],
                                                    dtype=np.int8 if len(cats) < 128 else np.int32))
        else:
            df.create_dataset(col, data=np.array(values))


def write_csr(g, name, n_rows, n_cols, density=0.05):
    nnz_per_row = max(1, int(n_cols * density))
    indices = np.concatenate([np.sort(RNG.choice(n_cols, nnz_per_row, replace=False)) for _ in range(n_rows)])
    x = g.create_group(name)
    x.attrs.update({"encoding-type": "csr_matrix", "encoding-version": "0.1.0", "shape": (n_rows, n_cols)})
    x.create_dataset("data", data=RNG.poisson(3, len(indices)).astype(np.float32) + 1)
    x.create_dataset("indices", data=indices.astype(np.int32))
    x.create_dataset("indptr", data=np.arange(0, (n_rows + 1) * nnz_per_row, nnz_per_row, dtype=np.int64))


def write_anndata(g, n_obs, var_names, obs_cols, obsm=("X_umap",), extra_var=None):
    g.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
    write_csr(g, "X", n_obs, len(var_names))
    write_df(g, "obs", [f"cell{i}" for i in range(n_obs)], obs_cols)
    write_df(g, "var", var_names, extra_var or {})
    om = g.create_group("obsm")
    for k in obsm:
        om.create_dataset(k, data=RNG.normal(size=(n_obs, 2)).astype(np.float32))
    g.create_group("layers")
    g.create_group("uns")


def pick(n, *values):
    return [values[i % len(values)] for i in range(n)]


def h5ad(path: Path, n_obs, species="human", tissue=("lung",), ensembl=False, assay="10x 3' v3",
         disease="normal", cell_types=("T cell", "B cell", "macrophage", "fibroblast"), spatial=False,
         n_vars=300, cxg=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    organism = {"human": "Homo sapiens", "mouse": "Mus musculus"}[species]
    names = genes(n_vars, species, ensembl)
    obs = {"cell_type": pick(n_obs, *cell_types), "donor_id": pick(n_obs, "D1", "D2", "D3")}
    if cxg:
        obs.update({"organism": [organism] * n_obs, "tissue": pick(n_obs, *tissue), "assay": [assay] * n_obs,
                    "disease": pick(n_obs, *disease) if isinstance(disease, tuple) else [disease] * n_obs})
    var_extra = {"feature_name": genes(n_vars, species)} if ensembl else None
    with h5py.File(path, "w") as f:
        write_anndata(f, n_obs, names, obs, ("X_umap", "spatial") if spatial else ("X_umap",), var_extra)


def tenx_h5(path: Path, n_cells, genome="GRCh38", species="human", extra=()):
    """A Cell Ranger 3+ filtered_feature_bc_matrix.h5."""
    path.parent.mkdir(parents=True, exist_ok=True)
    names = genes(200, species)
    types = ["Gene Expression"] * len(names)
    for ftype, n in extra:
        names += [f"{ftype.split()[0].upper()}_{i}" for i in range(n)]
        types += [ftype] * n
    with h5py.File(path, "w") as f:
        m = f.create_group("matrix")
        write_strings(m, "barcodes", [f"AAAC{i:08d}-1" for i in range(n_cells)])
        m.create_dataset("shape", data=np.array([len(names), n_cells], dtype=np.int32))
        nnz = 5
        m.create_dataset("data", data=RNG.poisson(2, n_cells * nnz).astype(np.int32) + 1)
        m.create_dataset("indices", data=np.tile(np.arange(nnz), n_cells).astype(np.int64))
        m.create_dataset("indptr", data=np.arange(0, (n_cells + 1) * nnz, nnz, dtype=np.int64))
        ft = m.create_group("features")
        write_strings(ft, "id", [f"F{i}" for i in range(len(names))])
        write_strings(ft, "name", names)
        write_strings(ft, "feature_type", types)
        write_strings(ft, "genome", [genome] * len(names))


def metrics(path: Path, rows: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows)
    path.write_text(",".join(keys) + "\n" + ",".join(f'"{rows[k]}"' for k in keys) + "\n")


def xenium(d: Path, organism, tissue, n_cells, panel, n_targets, name=None):
    d.mkdir(parents=True, exist_ok=True)
    species = "human" if organism == "Human" else "mouse"
    sparse_file(d / "transcripts.parquet", "4G")
    sparse_file(d / "cells.parquet", "300M")
    sparse_file(d / "morphology_focus/morphology_focus_0000.ome.tif", "6G")
    tenx_h5(d / "cell_feature_matrix.h5", min(n_cells, 3000), "Unknown", species,
            extra=(("Negative Control Probe", 20), ("Negative Control Codeword", 20)))
    (d / "experiment.xenium").write_text(json.dumps({
        "major_version": 5, "run_name": name or d.name, "region_name": "Region_1", "panel_name": panel,
        "panel_organism": organism, "panel_tissue_type": tissue, "panel_num_targets_predesigned": n_targets,
        "panel_num_targets_custom": 100, "instrument_sn": "XETG00248", "analysis_sw_version": "xenium-3.0.0.15",
        "num_cells": n_cells, "transcripts_per_cell": 212}, indent=1))


def visium(d: Path, species, n_spots, genome, hd=False):
    outs = d / "outs"
    if hd:
        tenx_h5(outs / "binned_outputs/square_008um/filtered_feature_bc_matrix.h5", 5000, genome, species)
        sparse_file(outs / "spatial/tissue_positions.parquet", "50M")
    else:
        tenx_h5(outs / "filtered_feature_bc_matrix.h5", n_spots, genome, species)
        sparse_file(outs / "spatial/tissue_positions.csv", "2M")
        sparse_file(outs / "spatial/tissue_hires_image.png", "10M")
    metrics(outs / "metrics_summary.csv", {"Sample ID": d.name, "Number of Spots Under Tissue": f"{n_spots:,}",
                                           "Mean Reads per Spot": "61,253", "Median Genes per Spot": "4,815",
                                           "Reference Path": f"/refs/refdata-gex-{genome}-2020-A"})


def ref10x(d: Path, genomes: list[str], same_as: Path | None = None):
    """A 10x reference package (mkref output)."""
    d.mkdir(parents=True, exist_ok=True)
    (d / "reference.json").write_text(json.dumps({"genomes": genomes, "mkref_version": "4.0.0",
                                                 "version": "2020-A"}, indent=1))
    for f, size in (("fasta/genome.fa", "300M"), ("fasta/genome.fa.fai", "1K"), ("genes/genes.gtf.gz", "40M"),
                    ("star/SA", "200M"), ("star/Genome", "300M")):
        sparse_file(d / f, size, same_as / f if same_as else None)
    sparse_file(d.parent / f"{d.name}.tar.gz", "80M")


def mirror(m: Path):
    """A typical lab data folder, in miniature: same
    folder names, same kinds of files, two samples instead of eighteen."""
    shutil.rmtree(m, ignore_errors=True)
    refs = m / "10x_refs"
    ref10x(refs / "cellranger/refdata-gex-GRCh38-2020-A", ["GRCh38"])
    ref10x(refs / "cellranger/refdata-gex-mm10-2020-A", ["mm10"])
    ref10x(refs / "spaceranger/refdata-gex-GRCh38-2020-A", ["GRCh38"],
           same_as=refs / "cellranger/refdata-gex-GRCh38-2020-A")  # the same reference twice
    arc = refs / "cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
    ref10x(arc, ["GRCh38"])
    atac_ref = refs / "cellranger-atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0"
    ref10x(atac_ref, ["GRCh38", "mm10"])
    (refs / "untar.sh").write_text("for f in `find . -name '*.tar.gz'`; do tar -xzf $f -C $(dirname $f); done\n")

    for f in ("hg38.knownGene.gtf.gz", "hg38.ensGene.bed.gz", "homer_hg38.tss"):
        sparse_file(m / "genes/hg38" / f, "20M")
    (m / "genes/hg38/info.org").write_text("* Location\nhttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/\n")
    for f in ("hg19.chrom.sizes", "hg38.chrom.sizes"):
        sparse_file(m / "genome_sizes" / f, "10K")
    for f in ("ENCFF356LFX_encode4_GRCh38_exclusion.bed.gz", "hg19-blacklist.v2.bed.gz", "mm10-blacklist.v2.bed.gz"):
        sparse_file(m / "exclusion_genomic_regions" / f, "10K")
    sparse_file(m / "logos/Logo Assets/Web Finals/logo.png", "200K")
    sparse_file(m / "logos/Logo Assets-20220915T141649Z-001.zip", "20M")

    for cond, n in (("untreated", 3887), ("sustained", 3120)):
        d = m / "faryabi_t-all_scATAC-seq/data" / cond
        d.mkdir(parents=True, exist_ok=True)
        (d / "barcodes.tsv").write_text("".join(f"AAAC{i:08d}-1\n" for i in range(n)))
        sparse_file(d / "fragments.tsv.gz", "1G")

    lp = m / "example_breast_multiome"
    atac_run = lp / "231005_A00000_0001_AHXXXXXXXX_ATAC"
    gex_run = lp / "231014_A00000_0002_BHXXXXXXXX_GEX"
    samples = (("S1_P1", 8), ("S2_P5", 9))
    for sample, s in samples:
        for lane in ("L001", "L002"):
            for read in ("I1", "R1", "R2", "R3"):
                name = f"{sample}_S{s}_{lane}_{read}_001.fastq.gz"
                sparse_file(atac_run / name, "2G")
                sparse_file(lp / "copied_atac" / sample / name, "2G", same_as=atac_run / name)  # a full copy
            for read in ("I1", "I2", "R1", "R2"):
                sparse_file(gex_run / f"{sample}_S{s}_{lane}_{read}_001.fastq.gz", "1.5G")
        lib = lp / "libraries" / f"{sample}.csv"
        lib.parent.mkdir(parents=True, exist_ok=True)
        lib.write_text(f"fastqs,sample,library_type\n{atac_run},{sample},Chromatin Accessibility\n"
                       f"{gex_run},{sample},Gene Expression\n")
        run = lp / "counts" / sample
        tenx_h5(run / "outs/filtered_feature_bc_matrix.h5", 900 + 100 * s, "GRCh38", "human", extra=(("Peaks", 300),))
        sparse_file(run / "outs/atac_fragments.tsv.gz", "3G")
        sparse_file(run / "outs/atac_possorted_bam.bam", "8G")
        (run / "_cmdline").write_text(f"cellranger-arc count --id={sample} --reference={arc} "
                                      f"--libraries=../libraries/{sample}.csv --localcores=42 --localmem=300\n")
        run = lp / "counts_atac" / sample
        outs = run / "outs"
        outs.mkdir(parents=True, exist_ok=True)
        sparse_file(outs / "fragments.tsv.gz", "3G")
        sparse_file(outs / "possorted_bam.bam", "8G")
        (outs / "peaks.bed").write_text("chr1\t100\t600\n")
        (outs / "summary.json").write_text(json.dumps({"cells_detected": 870 + 10 * s, "annotated_cells_GRCh38": 860,
                                                       "annotated_cells_mm10": 5 + s,
                                                       "cellranger-atac_version": "cellranger-atac-2.1.0"}))
        with h5py.File(outs / "filtered_peak_bc_matrix.h5", "w") as f:
            f.create_group("matrix").create_dataset("shape", data=np.array([4200, 870 + 10 * s], dtype=np.int32))
        (run / "_cmdline").write_text(f"cellranger-atac count --id={sample} --reference={atac_ref} "
                                      f"--chemistry=ARC-v1 --fastqs=../copied_atac/{sample} --localcores=42\n")
    (lp / "counts/count.sh").write_text("for libraryFile in ../libraries/*.csv; do cellranger-arc count ...; done\n")
    (lp / "info.org").write_text("* General information\nexample scMultiome project\n")
    sparse_file(lp / "SampleSummary.xlsx", "10K")

    sc = m / "scDEED_scRNA-seq"
    for name, files in (("Alveolar", ("ClassfiedCells.RData", "data_processing_Alveolar.Rmd", "loadings_n6_m0.7.Rds",
                                      "nneighbors6_mindist0.7_ClassifiedCells.Rds")),
                        ("CART", ("CART_STEAM.ipynb", "perplexity750_ClassifiedCells.Rds", "tSNEloadings_perplexity170.Rds")),
                        ("Hydra", ("Hydra_Seurat_Whole_Transcriptome.rds",)),
                        ("Samusik", ("ClassfiedCells.RData", "data_processing_samusik01.Rmd", "loadings_nneighbors160.Rds"))):
        for f in files:
            sparse_file(sc / name / f, "300M")

    xenium(m / "Xenium_V1_humanLung_Cancer_FFPE", "Human", "Multi", 162254, "hMulti_v1", 377,
           name="Human Multi-Tissue Cancer (FFPE)")


def main(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(root)
    try:
        build(root)
    finally:
        os.chdir(cwd)


def build(root: Path):

    # --- spatial transcriptomics
    xenium(Path("spatial/xenium/human_lung_cancer_ffpe_rep1"), "Human", "Lung", 162254, "hLung_v1", 289)
    xenium(Path("spatial/xenium/human_lung_cancer_ffpe_rep2"), "Human", "Lung", 158730, "hLung_v1", 289)
    xenium(Path("spatial/xenium/human_breast_idc"), "Human", "Breast", 167780, "hBreast_v1", 280)
    xenium(Path("spatial/xenium/mouse_brain_coronal"), "Mouse", "Brain", 36602, "mBrain_v1", 248)
    visium(Path("spatial/visium/mouse_brain_sagittal"), "mouse", 2696, "mm10")
    visium(Path("spatial/visium/human_lymph_node"), "human", 4035, "GRCh38")
    visium(Path("spatial/visium_hd/human_crc_p1"), "human", 0, "GRCh38", hd=True)
    for s in ("pdac_slide1", "pdac_slide2"):
        for f, size in (("Pancreas_exprMat_file.csv", "2G"), ("Pancreas_metadata_file.csv", "300M"),
                        ("Pancreas_fov_positions_file.csv", "10K")):
            sparse_file(Path(f"spatial/cosmx/{s}/{f}"), size)
    for f, size in (("cell_by_gene.csv", "400M"), ("cell_metadata.csv", "100M"), ("detected_transcripts.csv", "8G")):
        sparse_file(Path(f"spatial/merscope/mouse_liver_s1/{f}"), size)
    z = Path("spatial/spatialdata/human_lung_xenium_rep1.zarr")
    shutil.rmtree(z, ignore_errors=True)
    for sub in ("images/morphology_focus", "points/transcripts", "shapes/cell_boundaries", "tables/table/X"):
        (z / sub).mkdir(parents=True, exist_ok=True)
    (z / ".zattrs").write_text(json.dumps({"spatialdata_attrs": {"version": "0.2"}}))
    (z / "tables/table/X/.zarray").write_text(json.dumps({"shape": [162254, 389], "chunks": [1000, 389]}))

    # --- single cell
    h5ad(Path("single-cell/hlca_lung_atlas/hlca_core.h5ad"), 3000, tissue=("lung", "lung parenchyma"),
         assay="10x 3' v3", cell_types=("alveolar type 2 cell", "alveolar macrophage", "T cell", "fibroblast",
                                        "ciliated cell", "club cell"))
    h5ad(Path("single-cell/hlca_lung_atlas/hlca_core_epithelial.h5ad"), 900, tissue=("lung",),
         cell_types=("alveolar type 2 cell", "ciliated cell", "club cell"))
    Path("single-cell/hlca_lung_atlas_copy").mkdir(parents=True, exist_ok=True)
    shutil.copy2("single-cell/hlca_lung_atlas/hlca_core.h5ad", "single-cell/hlca_lung_atlas_copy/hlca_core.h5ad")
    h5ad(Path("single-cell/tabula_sapiens/ts_liver.h5ad"), 1200, tissue=("liver",), ensembl=True,
         assay="10x 3' v3", cell_types=("hepatocyte", "Kupffer cell", "endothelial cell"))
    h5ad(Path("single-cell/tabula_sapiens/ts_pancreas.h5ad"), 800, tissue=("exocrine pancreas", "pancreas"),
         ensembl=True, cell_types=("acinar cell", "pancreatic ductal cell", "type B pancreatic cell"))
    tenx_h5(Path("single-cell/pbmc_10k_v3/outs/filtered_feature_bc_matrix.h5"), 11769, "GRCh38", "human")
    sparse_file(Path("single-cell/pbmc_10k_v3/outs/web_summary.html"), "4M")
    metrics(Path("single-cell/pbmc_10k_v3/outs/metrics_summary.csv"),
            {"Estimated Number of Cells": "11,769", "Mean Reads per Cell": "54,502", "Median Genes per Cell": "3,051"})
    tenx_h5(Path("single-cell/mouse_kidney_multiome/outs/filtered_feature_bc_matrix.h5"), 4000, "mm10", "mouse",
            extra=(("Peaks", 300),))
    sparse_file(Path("single-cell/mouse_kidney_multiome/outs/atac_fragments.tsv.gz"), "3G")
    p = Path("single-cell/pdac_citeseq/pdac_cite.h5mu")
    p.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(p, "w") as f:
        write_df(f, "obs", [f"cell{i}" for i in range(1500)], {"sample": pick(1500, "P1", "P2")})
        mod = f.create_group("mod")
        write_anndata(mod.create_group("rna"), 1500, genes(300), {
            "organism": ["Homo sapiens"] * 1500, "tissue": ["pancreas"] * 1500,
            "cell_type": pick(1500, "ductal cell", "CAF", "T cell")})
        write_anndata(mod.create_group("prot"), 1500, [f"CD{i}" for i in range(1, 41)], {"sample": pick(1500, "a")})
    sparse_file(Path("single-cell/glioblastoma_snrna/seurat_obj.rds"), "3G")
    p = Path("single-cell/mouse_retina/retina_drop_seq.loom")
    p.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(p, "w") as f:
        f.create_dataset("matrix", data=RNG.poisson(1, (300, 1000)).astype(np.float32))
        ra = f.create_group("row_attrs")
        write_strings(ra, "Gene", genes(300, "mouse"))
        ca = f.create_group("col_attrs")
        write_strings(ca, "CellID", [f"c{i}" for i in range(1000)])

    # --- raw reads and images
    for r in ("R1", "R2", "I1"):
        sparse_file(Path(f"raw/GSE274673_pancreas_scrnaseq/S1_L001_{r}_001.fastq.gz"), "5G")
    sparse_file(Path("raw/breast_wgs/tumour.bam"), "60G")
    sparse_file(Path("raw/breast_wgs/normal.bam"), "40G")
    sparse_file(Path("histology/breast_he/slide_001.svs"), "2G")
    sparse_file(Path("histology/breast_he/slide_002.svs"), "2G")

    mirror(root / "mirror")

    # A folder only its owner can read, to test access checks in `get`.
    p = Path("restricted/restricted_lung")
    if p.exists():
        os.chmod(p, 0o700)
    h5ad(p / "patient_lung_raw.h5ad", 500, tissue=("lung",), cxg=False)
    os.chmod(p, 0o700)

    files = [q for q in root.rglob("*") if q.is_file() and not q.is_relative_to(root / ".data-steward")
             and not q.is_relative_to(root / "downloads")]
    apparent = sum(q.stat().st_size for q in files)
    used = sum(q.stat().st_blocks * 512 for q in files)
    print(f"Dummy data written to {root}: {len(files)} files, {apparent / 2**30:.0f} GiB apparent, "
          f"{used / 2**20:.0f} MiB on disk")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]).resolve())
