"""A typical lab data folder, in miniature (scripts/make-dummy-data.py mirror()): references, cellranger-atac outputs, fragment folders, project folders,
provenance links, markers, and moving datasets with `stew mv`."""

import importlib.util
import json
import os
from pathlib import Path

import pytest

from data_steward import cli, detect, provenance, vocab
from data_steward.db import Catalog, Filter

pytest.importorskip("h5py")

SCRIPT = Path(__file__).parent.parent / "scripts" / "make-dummy-data.py"


@pytest.fixture(scope="module")
def mirror(tmp_path_factory):
    spec = importlib.util.spec_from_file_location("make_dummy_data", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    root = (tmp_path_factory.mktemp("lab") / "data").resolve()
    mod.mirror(root)
    catalog = root / ".data-steward/catalog.db"  # the data folder itself, so nothing is moved
    Catalog.create(catalog)
    # recommended order: project folders as wholes first, then everything else
    assert stew(catalog, "scan", str(root / "scDEED_scRNA-seq"), "--each", "--yes") == 0
    assert stew(catalog, "scan", str(root), "--yes") == 0
    return root, catalog


def stew(catalog, *args):
    return cli.main(["--catalog", str(catalog), *args])


def test_joined_words_are_split():
    toks = vocab.tokens("Xenium_V1_humanLung_Cancer_FFPE")
    assert vocab.infer_organ(toks) == "lung" and vocab.infer_species(toks) == "human"
    assert "sc" in vocab.tokens("scATAC") and "atac" in vocab.tokens("scATAC")


def test_everything_is_recognised(mirror):
    root, catalog = mirror
    ds = {Path(d["path"]).relative_to(root).as_posix(): d for d in Catalog.open(catalog).list(Filter())}
    kinds = {p: (d["kind"], d["format"]) for p, d in ds.items()}
    refs = {p for p, (k, _) in kinds.items() if k == "reference"}
    assert refs == {
        "10x_refs/cellranger/refdata-gex-GRCh38-2020-A", "10x_refs/cellranger/refdata-gex-mm10-2020-A",
        "10x_refs/spaceranger/refdata-gex-GRCh38-2020-A",
        "10x_refs/cellranger-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0",
        "10x_refs/cellranger-atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0",
        "genes/hg38", "genome_sizes", "exclusion_genomic_regions"}
    assert ds["10x_refs/cellranger/refdata-gex-mm10-2020-A"]["species"] == "mouse"
    assert ds["10x_refs/cellranger-atac/refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0"]["species"] is None
    assert ds["genes/hg38"]["species"] == "human" and ds["exclusion_genomic_regions"]["species"] is None

    lp = "example_breast_multiome"
    arc = ds[f"{lp}/counts/S1_P1"]
    assert (arc["format"], arc["organ"], arc["species"], arc["n_cells"]) == ("cellranger-arc", "breast", "human", 1700)
    atac = ds[f"{lp}/counts_atac/S1_P1"]
    assert (atac["format"], atac["platform"], atac["species"], atac["n_cells"]) == \
        ("cellranger-atac", "10x-atac", "human", 950)
    assert atac["modalities"] == ["atac"]
    assert not any(p.startswith(f"{lp}/counts_atac/S1_P1/") for p in ds)  # no stray BAM entries

    frag = ds["faryabi_t-all_scATAC-seq/data/untreated"]
    assert (frag["format"], frag["modalities"], frag["n_cells"]) == ("atac-fragments", ["atac"], 3887)

    xen = ds["Xenium_V1_humanLung_Cancer_FFPE"]
    assert (xen["organ"], xen["species"]) == ("lung", "human")

    for p in ("Alveolar", "CART", "Hydra", "Samusik"):
        d = ds[f"scDEED_scRNA-seq/{p}"]
        assert d["format"] == "directory" and ".rds" in d["details"]["file_types"]
    assert not any(p.endswith((".Rds", ".rds")) for p in ds)  # inside the project folders, so skipped
    assert not any(p.startswith("logos") for p in ds)


def test_pipeline_outputs_are_linked_to_their_inputs(mirror):
    root, catalog = mirror
    cat = Catalog.open(catalog)
    lp = root / "example_breast_multiome"
    arc = cat.by_path(lp / "counts/S1_P1")
    up, _ = cat.relations(arc)
    got = {(k, Path(x["path"]).name) for k, x in up}
    assert got == {("derived-from", "231005_A00000_0001_AHXXXXXXXX_ATAC"),
                   ("derived-from", "231014_A00000_0002_BHXXXXXXXX_GEX"),
                   ("uses", "refdata-cellranger-arc-GRCh38-2020-A-2.0.0")}
    atac = cat.by_path(lp / "counts_atac/S2_P5")
    up, _ = cat.relations(atac)
    assert {(k, Path(x["path"]).name) for k, x in up} == {
        ("derived-from", "S2_P5"), ("uses", "refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0")}
    assert atac["details"]["pipeline"]["command"] == "cellranger-atac count"
    # the FASTQ folders take the species of what was made from them
    assert cat.by_path(lp / "copied_atac/S2_P5")["species"] == "human"
    assert cat.by_path(lp / "231014_A00000_0002_BHXXXXXXXX_GEX")["species"] == "human"


def test_provenance_reader(tmp_path):
    run = tmp_path / "counts/S1"
    run.mkdir(parents=True)
    lib = tmp_path / "libraries/S1.csv"
    lib.parent.mkdir()
    lib.write_text("fastqs,sample,library_type\n/abs/gex,S1,Gene Expression\nrel/atac,S1,Chromatin Accessibility\n")
    (run / "_cmdline").write_text("cellranger-arc count --id=S1 --reference /refs/arc --libraries=../libraries/S1.csv\n")
    p = provenance.read(run)
    assert p.inputs == [Path("/abs/gex"), tmp_path / "libraries/rel/atac"]
    assert p.references == [Path("/refs/arc")]
    assert provenance.read(tmp_path).empty()


def test_references_hidden_from_ls_by_default(mirror, capsys):
    _, catalog = mirror
    capsys.readouterr()
    stew(catalog, "ls", "--species", "human")
    out = capsys.readouterr()
    assert "refdata" not in out.out and "reference(s) not shown" in out.err
    stew(catalog, "ls", "--kind", "reference", "--paths")
    assert "refdata-gex-GRCh38-2020-A" in capsys.readouterr().out


def test_duplicates_found(mirror, capsys):
    root, catalog = mirror
    assert stew(catalog, "checksum", "--all", "--quick") == 0
    capsys.readouterr()
    stew(catalog, "duplicates")
    out = capsys.readouterr().out
    assert "refdata-gex-grch38-2020-a" in out.lower()  # cellranger and spaceranger copies
    assert "s1-p1:S1_P1_S8_" in out  # copied_atac duplicates the run folder


def test_marker_file_makes_a_folder_one_dataset(tmp_path):
    d = tmp_path / "project"
    (d / "sub").mkdir(parents=True)
    (d / "sub/a.h5ad").write_bytes(b"x")
    (d / detect.MARKER).write_text(json.dumps({"name": "My Project", "organ": "lungs", "species": "Homo sapiens",
                                               "modalities": ["scrna"], "study": "GSE1", "tags": ["pilot"]}))
    found = detect.scan(tmp_path).found
    assert len(found) == 1
    c = found[0]
    assert (c.name, c.organ, c.species, c.modalities) == ("my-project", "lung", "human", {"rna"})
    catalog = tmp_path / "c/.data-steward/catalog.db"
    Catalog.create(catalog)
    assert stew(catalog, "scan", str(tmp_path), "--yes") == 0
    got = Catalog.open(catalog).find("my-project")
    assert (got["study"], got["tags"]) == ("GSE1", ["pilot"])


def test_mv_moves_folder_updates_catalogue_and_links(mirror, tmp_path, capsys):
    root, catalog = mirror
    cat = Catalog.open(catalog)
    lp = root / "example_breast_multiome"
    arc = cat.by_path(lp / "counts/S1_P1")
    work = tmp_path / "work"
    work.mkdir()
    assert stew(catalog, "get", str(arc["id"]), "--to", str(work)) == 0
    n_under = len(cat.under(str(lp)))
    target = root / "projects/example_breast_multiome"

    assert stew(catalog, "mv", str(lp), str(target), "--parents", "--dry-run") == 0
    assert lp.is_dir() and not target.exists()

    assert stew(catalog, "mv", str(lp), str(target), "--parents", "--leave-link", "--yes") == 0
    assert (target / "counts/S1_P1/outs").is_dir()
    assert lp.is_symlink() and os.readlink(lp) == str(target)
    moved = cat.by_id(arc["id"])
    assert moved["path"] == str(target / "counts/S1_P1")
    assert len(cat.under(str(target))) == n_under and cat.under(str(lp)) == []
    assert os.readlink(work / "S1_P1") == moved["path"]  # the handed-out link follows
    up, _ = cat.relations(moved)
    assert len(up) == 3  # relations are by id, so they survive

    assert stew(catalog, "mv", str(target), str(root / "other" / "x")) == 1
    assert "--parents" in capsys.readouterr().err

    # a trailing slash means "into this folder", even when it doesn't exist yet
    xen = root / "Xenium_V1_humanLung_Cancer_FFPE"
    assert stew(catalog, "mv", str(xen), str(root / "public/10x") + "/", "--parents", "--yes") == 0
    assert (root / "public/10x/Xenium_V1_humanLung_Cancer_FFPE/experiment.xenium").is_file()


def test_mv_refuses_without_write_permission(tmp_path, capsys):
    root = tmp_path / "d"
    locked = root / "locked"
    locked.mkdir(parents=True)
    (locked / "a.h5ad").write_bytes(b"x")
    catalog = tmp_path / "c/.data-steward/catalog.db"
    Catalog.create(catalog)
    assert stew(catalog, "add", str(locked), "--whole", "--yes") == 0
    os.chmod(locked, 0o555)
    try:
        if os.access(locked, os.W_OK):
            pytest.skip("running as root")
        assert stew(catalog, "mv", str(locked), str(root / "elsewhere") + "/", "--parents", "--yes") == 1
        assert "no permission to move it" in capsys.readouterr().err
        assert locked.is_dir()
    finally:
        os.chmod(locked, 0o755)


def test_add_moves_into_the_data_folder_and_leaves_a_link(tmp_path, capsys):
    home = tmp_path / "data"
    catalog = home / ".data-steward/catalog.db"
    Catalog.create(catalog)
    mine = tmp_path / "me/project/GSE123"
    mine.mkdir(parents=True)
    (mine / "counts.h5ad").write_bytes(b"x")
    assert stew(catalog, "add", str(mine), "--whole", "--source", "GSE123", "--dry-run") == 0
    assert "public/geo/GSE123" in capsys.readouterr().out and mine.is_dir() and not mine.is_symlink()
    assert stew(catalog, "add", str(mine), "--whole", "--source", "GSE123", "--yes") == 0
    dest = home / "public/geo/GSE123"
    assert (dest / "counts.h5ad").is_file()
    assert mine.is_symlink() and os.readlink(mine) == str(dest)
    assert (mine / "counts.h5ad").read_bytes() == b"x"  # the old path still works
    assert Catalog.open(catalog).find("gse123")["path"] == str(dest)
    import stat
    mode = (home / "public/geo").stat().st_mode
    assert mode & stat.S_IWGRP  # shared
    probe = tmp_path / "probe"
    probe.mkdir()
    try:
        os.chmod(probe, 0o3775)
        special = True
    except OSError:
        special = False  # e.g. the Nix build sandbox
    if special:
        assert mode & stat.S_ISVTX and mode & stat.S_ISGID  # sticky, inherits the group

    # no source: projects/<name>; --to picks the folder; --in-place stays
    other = tmp_path / "me/pilot.h5ad"
    other.write_bytes(b"y")
    assert stew(catalog, "add", str(other), "--yes") == 0
    assert (home / "projects/pilot.h5ad").is_file() and other.is_symlink()
    third = tmp_path / "me/third"
    third.mkdir()
    (third / "a.h5ad").write_bytes(b"z")
    assert stew(catalog, "add", str(third), "--to", "projects/study_x", "--yes") == 0
    assert (home / "projects/study_x/a.h5ad").is_file()  # the dataset found inside the folder
    fourth = tmp_path / "me/fourth.h5ad"
    fourth.write_bytes(b"w")
    assert stew(catalog, "add", str(fourth), "--in-place", "--yes") == 0
    assert not fourth.is_symlink()

    # data already in the data folder is not moved
    inside = home / "projects/already"
    inside.mkdir(parents=True)
    (inside / "b.h5ad").write_bytes(b"v")
    assert stew(catalog, "add", str(inside), "--yes") == 0
    assert not inside.is_symlink()


def test_scan_moves_each_dataset_and_keeps_the_rest(tmp_path):
    home = tmp_path / "data"
    catalog = home / ".data-steward/catalog.db"
    Catalog.create(catalog)
    proj = tmp_path / "me/proj/data"
    (proj / "atlas").mkdir(parents=True)
    (proj / "atlas/lung.h5ad").write_bytes(b"x")
    (proj / "download.sh").write_text("curl ...\n")
    assert stew(catalog, "scan", str(proj), "--yes") == 0
    assert (home / "projects/lung.h5ad").is_file()
    assert (proj / "atlas/lung.h5ad").is_symlink()
    assert (proj / "download.sh").is_file() and not (proj / "download.sh").is_symlink()


def test_add_refuses_when_it_cannot_move(tmp_path, capsys):
    home = tmp_path / "data"
    catalog = home / ".data-steward/catalog.db"
    Catalog.create(catalog)
    (home / "projects").mkdir()
    os.chmod(home / "projects", 0o555)
    f = tmp_path / "me.h5ad"
    f.write_bytes(b"x")
    try:
        if os.access(home / "projects", os.W_OK):
            pytest.skip("running as root")
        assert stew(catalog, "add", str(f), "--yes") == 1
        err = capsys.readouterr().err
        assert "permission" in err and "--in-place" in err
        assert f.is_file() and not f.is_symlink()
        assert Catalog.open(catalog).list(Filter()) == []
    finally:
        os.chmod(home / "projects", 0o755)


def test_mtx_triplets_and_count_tables_are_datasets(tmp_path):
    import gzip
    geo = tmp_path / "GSE1"
    s1 = geo / "sample1"
    s1.mkdir(parents=True)
    with gzip.open(s1 / "GSM1_matrix.mtx.gz", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n%\n300 1234 5000\n1 1 3\n")
    genes = "".join(f"ENSG{i:011d}\tG{i}\tGene Expression\n" for i in range(300))
    with gzip.open(s1 / "GSM1_features.tsv.gz", "wt") as f:
        f.write(genes)
    with gzip.open(s1 / "GSM1_barcodes.tsv.gz", "wt") as f:
        f.write("AAAC-1\n")
    tables = geo / "bulk"
    tables.mkdir()
    (tables / "GSE1_raw_counts.csv").write_text("gene,s1,s2\nACTB,1,2\n")
    (tables / "GSE1_sample_info.csv").write_text("sample,condition\ns1,ctrl\n")
    (tables / "README.txt").write_text("notes")
    notes = geo / "notes"
    notes.mkdir()
    (notes / "samples.csv").write_text("a,b\n")  # tables, but nothing says they hold measurements
    found = {c.path.name: c for c in detect.scan(geo).found}
    assert set(found) == {"sample1", "bulk"}
    assert found["sample1"].format == "mtx" and found["bulk"].format == "count-tables"
    from data_steward import contents
    m = contents.read(s1, "mtx")
    assert (m.n_features, m.n_cells, m.species) == (300, 1234, "human")
    assert contents.read(tables, "count-tables").details["tables"]["GSE1_raw_counts.csv"]["columns"] == 3


def test_add_decides_from_what_it_is_pointed_at(tmp_path, capsys):
    home = tmp_path / "data"
    catalog = home / ".data-steward/catalog.db"
    Catalog.create(catalog)
    # a folder of datasets: each registered separately, other files stay
    proj = tmp_path / "me/proj"
    for n in ("a", "b"):
        (proj / n).mkdir(parents=True)
        (proj / n / f"{n}.h5ad").write_bytes(b"x")
    (proj / "run.sh").write_text("echo\n")
    assert stew(catalog, "add", str(proj), "--name", "x", "--yes") == 1  # --name needs one dataset
    assert stew(catalog, "add", str(proj), "--yes") == 0
    assert (home / "projects/a.h5ad").is_file() and (home / "projects/b.h5ad").is_file()
    assert (proj / "a/a.h5ad").is_symlink() and not (proj / "run.sh").is_symlink()
    # --whole keeps a folder together
    whole = tmp_path / "me/analysis"
    for n in ("c", "d"):
        (whole / f"{n}.h5ad").parent.mkdir(parents=True, exist_ok=True)
        (whole / f"{n}.h5ad").write_bytes(b"y")
    assert stew(catalog, "add", str(whole), "--whole", "--yes") == 0
    assert (home / "projects/analysis/c.h5ad").is_file() and whole.is_symlink()
    # a folder with nothing stew recognises: the folder as one dataset
    misc = tmp_path / "me/misc"
    misc.mkdir()
    (misc / "notes.pdf").write_bytes(b"z")
    assert stew(catalog, "add", str(misc), "--yes") == 0
    assert (home / "projects/misc/notes.pdf").is_file()
    # a recognised dataset folder stays one dataset
    xen = tmp_path / "me/run1"
    xen.mkdir()
    (xen / "experiment.xenium").write_text("{}")
    assert stew(catalog, "add", str(xen), "--yes") == 0
    assert Catalog.open(catalog).find("run1")["format"] == "xenium-bundle"
    # `scan` is the same command
    # one dataset inside: that dataset is registered, the notes stay; `scan` is the same command
    more = tmp_path / "me/more"
    more.mkdir()
    (more / "e.h5ad").write_bytes(b"w")
    (more / "README.md").write_text("about e\n")
    assert stew(catalog, "scan", str(more), "--yes") == 0
    assert (home / "projects/e.h5ad").is_file() and not (more / "README.md").is_symlink()
