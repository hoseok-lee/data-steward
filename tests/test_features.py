"""Tests for reading file contents, relations, checksums, saved queries,
search, the web page, the extra download sources and background downloads."""

import csv
import importlib.util
import io
import json
import shutil
import threading
import urllib.request
from pathlib import Path

import pytest

from data_steward import checksum, cli, contents, sources, web
from data_steward.db import Catalog

from test_stew import Server, touch  # noqa: F401

h5py = pytest.importorskip("h5py")

SCRIPT = Path(__file__).parent.parent / "scripts" / "make-dummy-data.py"


@pytest.fixture(scope="module")
def dummy():
    spec = importlib.util.spec_from_file_location("make_dummy_data", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def tree(tmp_path_factory, dummy):
    """The full dummy tree, catalogued once and shared by the tests below."""
    root = (tmp_path_factory.mktemp("lab") / "data").resolve()
    dummy.main(root)
    catalog = root / ".data-steward/catalog.db"
    Catalog.create(catalog)
    assert cli.main(["--catalog", str(catalog), "scan", str(root), "--yes"]) == 0
    return root, catalog


def stew(catalog, *args):
    return cli.main(["--catalog", str(catalog), *args])


# --- reading contents ---------------------------------------------------------------------------------

def test_reads_cellxgene_style_h5ad(tmp_path, dummy):
    p = tmp_path / "a.h5ad"
    dummy.h5ad(p, 300, tissue=("lung", "lung parenchyma"), cell_types=("T cell", "fibroblast"))
    c = contents.read(p, "anndata")
    assert (c.n_cells, c.n_features, c.species, c.organ, c.platform) == (300, 300, "human", "lung", "10x-chromium")
    assert c.modalities == {"rna"}
    assert c.details["cell_types"] == ["T cell", "fibroblast"] and c.details["diseases"] == ["normal"]


def test_species_from_gene_names_alone(tmp_path, dummy):
    p = tmp_path / "m.h5ad"
    dummy.h5ad(p, 50, species="mouse", cxg=False)
    assert contents.read(p, "anndata").species == "mouse"
    p = tmp_path / "e.h5ad"
    dummy.h5ad(p, 50, ensembl=True, cxg=False)
    assert contents.read(p, "anndata").species == "human"
    assert contents.species_from_genes(["ENSMUSG00000000001"] * 30) == "mouse"


def test_reads_other_formats(tree):
    root, _ = tree
    c = contents.read(root / "single-cell/pdac_citeseq/pdac_cite.h5mu", "mudata")
    assert c.n_cells == 1500 and c.modalities == {"rna", "protein"} and c.organ == "pancreas"
    c = contents.read(root / "single-cell/mouse_kidney_multiome", "cellranger-arc")
    assert (c.species, c.n_cells) == ("mouse", 4000) and {"rna", "atac"} <= c.modalities
    c = contents.read(root / "spatial/xenium/human_breast_idc", "xenium-bundle")
    assert (c.n_cells, c.n_features, c.species, c.organ) == (167780, 380, "human", "breast")
    c = contents.read(root / "spatial/visium/mouse_brain_sagittal", "spaceranger")
    assert (c.n_cells, c.species) == (2696, "mouse") and c.details["metrics"]["Mean Reads per Spot"] == "61,253"
    c = contents.read(root / "single-cell/mouse_retina/retina_drop_seq.loom", "loom")
    assert (c.n_cells, c.n_features, c.species) == (1000, 300, "mouse")
    c = contents.read(root / "spatial/spatialdata/human_lung_xenium_rep1.zarr", "spatialdata")
    assert (c.n_cells, c.n_features) == (162254, 389) and "spatial" in c.modalities
    c = contents.read(root / "raw/GSE274673_pancreas_scrnaseq", "fastq")
    assert c.details["fastq"]["reads"] == {"R1": 1, "R2": 1, "I1": 1}
    bad = root.parent / "broken.h5ad"
    bad.write_bytes(b"not hdf5")
    assert contents.read(bad, "anndata").error


def test_scan_records_contents(tree):
    _, catalog = tree
    cat = Catalog.open(catalog)
    d = cat.find("hlca-core")
    assert (d["n_cells"], d["species"], d["organ"]) == (3000, "human", "lung")
    assert d["details"]["assays"] == ["10x 3' v3"]
    # a path guess ("pbmc" -> blood) is kept, the genome fills in the species
    pbmc = cat.find("pbmc-10k-v3")
    assert (pbmc["species"], pbmc["organ"], pbmc["n_cells"]) == ("human", "blood", 11769)


def test_inspect_fills_only_missing_fields_unless_overwrite(tree, capsys):
    _, catalog = tree
    cat = Catalog.open(catalog)
    stew(catalog, "set", "ts-liver", "--organ", "kidney", "--species", "")
    assert stew(catalog, "inspect", "ts-liver") == 0
    d = cat.find("ts-liver")
    assert (d["organ"], d["species"]) == ("kidney", "human")
    assert stew(catalog, "inspect", "ts-liver", "--overwrite") == 0
    assert cat.find("ts-liver")["organ"] == "liver"
    capsys.readouterr()
    assert stew(catalog, "inspect", "--all") == 0
    assert "nothing to read" in capsys.readouterr().out  # e.g. the .svs slides


# --- relations ---------------------------------------------------------------------------------------

def test_relations(tree, capsys):
    _, catalog = tree
    assert stew(catalog, "relate", "hlca-core-epithelial", "hlca-core") == 0
    assert stew(catalog, "relate", "hlca-core", "hlca-core-epithelial") == 1  # would be a loop
    assert "loop" in capsys.readouterr().err
    capsys.readouterr()
    stew(catalog, "info", "hlca-core")
    assert "used by     #" in capsys.readouterr().out
    assert stew(catalog, "unrelate", "hlca-core-epithelial", "hlca-core") == 0
    assert stew(catalog, "unrelate", "hlca-core-epithelial", "hlca-core") == 1


# --- checksums ---------------------------------------------------------------------------------------

def test_checksums_duplicates_and_verify(tree, capsys):
    root, catalog = tree
    small = ["hlca-core", "hlca-core-2", "ts-liver", "pdac-cite"]
    names = {d["name"] for d in Catalog.open(catalog).list()}
    small = [n for n in small if n in names] + [n for n in names if n.startswith("hlca-lung-atlas-copy")]
    assert stew(catalog, "checksum", *small) == 0
    capsys.readouterr()
    stew(catalog, "duplicates")
    out = capsys.readouterr().out
    assert "Whole datasets stored more than once" in out and "hlca-core" in out
    assert stew(catalog, "verify", "hlca-core") == 0
    copy = root / "single-cell/hlca_lung_atlas_copy/hlca_core.h5ad"
    with open(copy, "r+b") as f:
        f.seek(100)
        f.write(b"changed")
    copy_name = next(n for n in small if n != "hlca-core" and Catalog.open(catalog).find(n)["path"] == str(copy))
    capsys.readouterr()
    assert stew(catalog, "verify", copy_name) == 1
    assert "1 changed" in capsys.readouterr().out


def test_quick_and_full_hashes_differ_in_kind(tmp_path):
    p = touch(tmp_path / "x.bin", b"a" * (3 << 20))
    q, f = checksum.file_hash(p, True), checksum.file_hash(p, False)
    assert q.startswith("quick:") and f.startswith("full:")
    m = checksum.manifest(tmp_path, quick=False)
    assert m[0][0] == "x.bin" and checksum.fingerprint(m).startswith("full:")
    assert checksum.manifest(p, quick=True)[0][0] == ""  # single files have no relative path


# --- search, saved queries, CSV -----------------------------------------------------------------------

def test_search_saved_queries_and_csv(tree, capsys):
    _, catalog = tree
    capsys.readouterr()
    assert stew(catalog, "search", "alveolar", "--paths") == 0
    found = capsys.readouterr().out.split()
    assert any("hlca_core.h5ad" in p for p in found)  # matched inside the cell types read from the file
    assert stew(catalog, "ls", "-m", "spatial", "--species", "human", "--save", "hs") == 0
    first = capsys.readouterr().out
    assert stew(catalog, "ls", "-q", "hs") == 0
    assert capsys.readouterr().out == first
    assert stew(catalog, "ls", "-q", "hs", "--organ", "lung", "--csv") == 0
    rows = list(csv.DictReader(io.StringIO(capsys.readouterr().out)))
    assert rows and {r["organ"] for r in rows} == {"lung"} and all("spatial" in r["modalities"] for r in rows)
    stew(catalog, "saved")
    assert "--species human --modality spatial" in capsys.readouterr().out
    assert stew(catalog, "saved", "--rm", "hs") == 0
    assert stew(catalog, "ls", "-q", "hs") == 1


# --- web page ----------------------------------------------------------------------------------------

def test_web_page_is_read_only_json(tree):
    _, catalog = tree
    web.Handler.store = web.Store(catalog)
    httpd = web.ThreadingHTTPServer(("127.0.0.1", 0), web.Handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"
    try:
        page = urllib.request.urlopen(base + "/").read().decode()
        assert "Dataset catalogue" in page and "<script>" in page
        data = json.load(urllib.request.urlopen(base + "/api/datasets"))
        assert len(data["datasets"]) > 20 and "summary" in data
        one = json.load(urllib.request.urlopen(f"{base}/api/datasets/{data['datasets'][0]['id']}"))
        assert "relations" in one and "exists" in one
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(base + "/api/datasets/99999")
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(urllib.request.Request(base + "/api/datasets", data=b"x", method="POST"))
        with pytest.raises(Exception):
            web.Handler.store.cat.conn.execute("DELETE FROM datasets")
    finally:
        httpd.shutdown()


# --- sources ------------------------------------------------------------------------------------------

@pytest.mark.parametrize("text,key", [
    ("SRP000001", "sra:SRP000001"), ("prjna1147630", "sra:PRJNA1147630"), ("ERR123456", "sra:ERR123456"),
    ("https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1147630", "sra:PRJNA1147630"),
    ("https://www.ebi.ac.uk/ena/browser/view/PRJEB1234", "sra:PRJEB1234"),
    ("E-MTAB-5061", "biostudies:E-MTAB-5061"),
    ("https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5061", "biostudies:E-MTAB-5061"),
    ("figshare:5968960", "figshare:5968960"), ("10.6084/m9.figshare.5968960.v2", "figshare:5968960"),
    ("https://figshare.com/articles/dataset/Single-cell_RNA-seq/5968960/2", "figshare:5968960"),
    ("https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293",
     "cellxgene:6f6d381a-7701-4781-935c-db10d30de293"),
    ("hca:74b6d569-3b11-42ef-b6b1-a0454522b4a0", "hca:74b6d569-3b11-42ef-b6b1-a0454522b4a0"),
    ("https://data.humancellatlas.org/hca-bio-networks/lung/projects/74b6d569-3b11-42ef-b6b1-a0454522b4a0",
     "hca:74b6d569-3b11-42ef-b6b1-a0454522b4a0"),
    ("GSE274673:raw", "geo:GSE274673:raw"),
])
def test_parse_more_sources(text, key):
    assert sources.parse(text).key == key


def test_title_naming_two_species_leaves_species_unset(lab, monkeypatch, capsys):
    monkeypatch.setitem(sources.RESOLVERS, "geo", lambda src: sources.Remote(
        title="A single-cell map of the human and mouse pancreas", files=[sources.RemoteFile("a.tar", "https://x", 1)]))
    assert stew(lab, "download", "GSE84133", "--dry-run") == 0
    out = capsys.readouterr().out
    assert "organ pancreas, species -" in out


def test_raw_flag():
    assert sources.parse("GSE1", raw=True).key == "geo:GSE1:raw"
    assert sources.parse("GSE1", raw=True).folder == "GSE1-raw"
    assert sources.parse("SRP1", raw=True).key == "sra:SRP1"  # SRA is always reads
    with pytest.raises(sources.SourceError):
        sources.parse("zenodo:1", raw=True)


def test_resolvers_parse_api_answers(monkeypatch):
    answers = {
        "ena": [{"run_accession": "SRR1", "fastq_ftp": "ftp.sra.ebi.ac.uk/a/SRR1_1.fastq.gz;ftp.sra.ebi.ac.uk/a/SRR1_2.fastq.gz",
                 "fastq_bytes": "10;20", "fastq_md5": "aa;bb", "scientific_name": "Homo sapiens",
                 "library_strategy": "RNA-Seq", "study_title": "T", "study_accession": "PRJNA1"}],
        "cellxgene": {"name": "Lung atlas", "doi": "10.1/x", "datasets": [
            {"dataset_id": "abcdef12-0000", "title": "Core", "organism": [{"label": "Homo sapiens"}],
             "tissue": [{"label": "lung"}], "assay": [{"label": "10x 3' v3"}], "spatial": None,
             "assets": [{"filetype": "H5AD", "url": "https://x/a.h5ad", "filesize": 5}]}]},
        "figshare": {"title": "F", "doi": "10.6084/m9.figshare.1", "files": [
            {"name": "a.csv", "size": 3, "download_url": "https://x/1", "computed_md5": "cc"},
            {"name": "link", "is_link_only": True, "download_url": "https://x/2"}]},
    }

    def fake_json(url):
        for k, v in answers.items():
            if k in url or (k == "ena" and "ebi.ac.uk/ena" in url):
                return v
        raise AssertionError(url)

    monkeypatch.setattr(sources, "get_json", fake_json)
    r = sources.resolve(sources.parse("SRP1"))
    assert [(f.name, f.size, f.md5) for f in r.files] == [("SRR1_1.fastq.gz", 10, "aa"), ("SRR1_2.fastq.gz", 20, "bb")]
    assert r.files[0].url.startswith("https://ftp.sra.ebi.ac.uk/") and r.modalities == {"rna"}
    r = sources.resolve(sources.parse("cellxgene:6f6d381a-7701-4781-935c-db10d30de293"))
    assert (r.species, r.organ, r.platform, r.files[0].name) == ("human", "lung", "10x-chromium", "core-abcdef12.h5ad")
    r = sources.resolve(sources.parse("figshare:1"))
    assert [f.name for f in r.files] == ["a.csv"]


def test_geo_raw_reads_go_through_sra(monkeypatch):
    monkeypatch.setattr(sources, "_geo_fields", lambda acc: {
        "Series_title": ["X"], "Series_sample_taxid": ["9606"], "Series_type": ["Expression profiling by high throughput sequencing"],
        "Series_relation": ["BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA9",
                            "SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRP9"]})
    monkeypatch.setattr(sources, "resolve_sra", lambda src: sources.Remote(
        files=[sources.RemoteFile("SRR9.fastq.gz", "https://x", 1)], modalities={"rna"}))
    r = sources.resolve(sources.parse("GSE9", raw=True))
    assert [f.name for f in r.files] == ["SRR9.fastq.gz"]
    assert set(r.parents) == {"sra:PRJNA9", "sra:SRP9"}


# --- downloads -----------------------------------------------------------------------------------------

@pytest.fixture
def lab(tmp_path):
    catalog = tmp_path / "data/.data-steward/catalog.db"
    Catalog.create(catalog)
    return catalog


def test_raw_download_blocked_when_its_sra_study_is_here(lab, monkeypatch, capsys):
    Catalog.open(lab).insert({"name": "reads", "path": "/r", "source": "SRP9", "source_key": "sra:SRP9",
                              "modalities": []})
    monkeypatch.setitem(sources.RESOLVERS, "geo", lambda src: sources.Remote(
        parents=["sra:SRP9"], files=[sources.RemoteFile("a.fastq.gz", "https://x", 1)]))
    assert stew(lab, "download", "GSE9", "--raw", "--yes") == 1
    assert "same as) SRP9" in capsys.readouterr().err


def test_download_registers_parts_of_a_record(lab, tmp_path, dummy, monkeypatch, capsys):
    served = tmp_path / "served"
    dummy.h5ad(served / "lung.h5ad", 100, tissue=("lung",))
    dummy.h5ad(served / "liver.h5ad", 120, tissue=("liver",))
    server = Server(served)
    monkeypatch.setitem(sources.RESOLVERS, "zenodo", lambda src: sources.Remote(
        title="Two organs", files=[sources.RemoteFile(n, f"{server.url}/{n}") for n in ("lung.h5ad", "liver.h5ad")]))
    try:
        assert stew(lab, "download", "zenodo:77", "--yes") == 0
    finally:
        server.httpd.shutdown()
    cat = Catalog.open(lab)
    whole = cat.find("zenodo-77")
    assert whole["source"] == "zenodo:77" and whole["description"] == "Two organs"
    _, parts = cat.relations(whole)
    assert sorted((p["name"], p["organ"], p["n_cells"]) for _, p in parts) == [("liver", "liver", 120),
                                                                               ("lung", "lung", 100)]
    assert all(k == "part-of" for k, _ in parts)


def test_background_download_hands_over_to_tmux(lab, tmp_path, monkeypatch, capsys):
    served = touch(tmp_path / "served/x.h5ad", b"x" * 10).parent
    server = Server(served)
    calls = []
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/tmux")
    monkeypatch.setattr(cli.subprocess, "run", lambda cmd, **kw: calls.append((cmd, kw)))
    try:
        assert stew(lab, "download", f"{server.url}/x.h5ad", "--yes", "--background", "--organ", "lung") == 0
    finally:
        server.httpd.shutdown()
    (cmd, kw), = calls
    assert cmd[:4] == ["tmux", "new-session", "-d", "-s"] and cmd[4].startswith("stew-dl-")
    assert "download" in cmd[5] and "--background" not in cmd[5] and "--yes" in cmd[5]
    assert kw["env"]["STEW_BACKGROUND_LOG"].endswith(".log")
    assert cmd[5].startswith("STEW_BACKGROUND_LOG=")
    assert Catalog.open(lab).downloads() == []  # nothing recorded until the tmux job starts


def test_finished_background_downloads_are_announced_once(lab, capsys):
    cat = Catalog.open(lab)
    cat.conn.execute("INSERT INTO downloads (source_key, source, user, host, pid, dest, status, started_at, background, "
                     "log) VALUES ('zenodo:5', 'zenodo:5', ?, 'h', 1, '/d', 'done', 'now', 1, '/l')", (cat.user,))
    stew(lab, "ls")
    assert "your background download of zenodo:5 finished" in capsys.readouterr().err
    stew(lab, "ls")
    assert "background download" not in capsys.readouterr().err


def test_sra_falls_back_to_ncbi_when_ena_has_nothing(monkeypatch):
    monkeypatch.setattr(sources, "get_json", lambda url: [])
    csv_text = ("Run,LibraryStrategy,SRAStudy,ScientificName\n"
                "SRR1,ATAC-seq,SRP9,Homo sapiens\nSRR2,ATAC-seq,SRP9,Homo sapiens\n")
    monkeypatch.setattr(sources, "get_text", lambda url: csv_text if "runinfo" in url else "")
    monkeypatch.setattr(sources, "head", lambda url: (1000, None))
    r = sources.resolve(sources.parse("SRP9"))
    assert [f.name for f in r.files] == ["SRR1.sra", "SRR2.sra"]
    assert r.files[0].url == "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR1/SRR1"
    assert (r.species, r.modalities, r.study) == ("human", {"atac"}, "SRP9") and "fasterq-dump" in r.note
    monkeypatch.setattr(sources, "get_text", lambda url: "Run,SRAStudy\n")
    with pytest.raises(sources.SourceError, match="neither ENA nor NCBI"):
        sources.resolve(sources.parse("PRJNA1"))


def test_links_prune(lab, tmp_path, capsys):
    f = touch(tmp_path / "data/x.h5ad")
    assert stew(lab, "add", str(f), "--yes") == 0
    work = tmp_path / "work"
    work.mkdir()
    assert stew(lab, "get", "x", "--to", str(work)) == 0
    assert stew(lab, "get", "x", "--to", str(work), "--as", "y") == 0
    (work / "x.h5ad").unlink()
    capsys.readouterr()
    assert stew(lab, "links", "--prune") == 0
    assert "1 link record(s) pruned" in capsys.readouterr().out
    assert [lk["link_path"] for lk in Catalog.open(lab).links()] == [str(work / "y")]


def test_relate_suggest(lab, tmp_path, capsys):
    cat = Catalog.open(lab)
    cat.insert({"name": "gse1", "path": "/d/gse1", "source": "GSE1", "source_key": "geo:GSE1", "modalities": []})
    cat.insert({"name": "vdx", "path": "/d/vdx", "description": "Bioconductor package of GSE1", "modalities": []})
    cat.insert({"name": "s1", "path": "/d/s1", "source": "GSM9", "source_key": "geo:GSM9", "study": "GSE1",
                "modalities": []})
    capsys.readouterr()
    assert stew(lab, "relate", "--suggest") == 0
    out = capsys.readouterr().out
    assert "vdx" in out and "derived-from" in out and "part-of" in out
    assert stew(lab, "relate", "--suggest", "--apply") == 0
    _, down = Catalog.open(lab).relations(Catalog.open(lab).find("gse1"))
    assert sorted((k, d["name"]) for k, d in down) == [("derived-from", "vdx"), ("part-of", "s1")]
    capsys.readouterr()
    stew(lab, "relate", "--suggest")
    assert "No new relations" in capsys.readouterr().out


def test_download_missing_tops_up_a_catalogued_record(lab, tmp_path, monkeypatch, capsys):
    served = tmp_path / "served"
    for n, body in (("a.zip", b"a" * 10), ("b.zip", b"b" * 20), ("c.zip", b"c" * 30)):
        touch(served / n, body)
    server = Server(served)
    monkeypatch.setitem(sources.RESOLVERS, "zenodo", lambda src: sources.Remote(
        title="Tiles", files=[sources.RemoteFile(n, f"{server.url}/{n}", s)
                              for n, s in (("a.zip", 10), ("b.zip", 20), ("c.zip", 30))]))
    try:
        assert stew(lab, "download", "zenodo:88", "--files", "a.zip", "--files", "b.zip", "--yes") == 0
        cat = Catalog.open(lab)
        d = cat.find("zenodo-88")
        assert d["n_files"] == 2
        capsys.readouterr()
        assert stew(lab, "download", "zenodo:88", "--yes") == 1          # still refused without --missing
        assert "already downloaded" in capsys.readouterr().err
        assert stew(lab, "download", "zenodo:88", "--missing", "--name", "x", "--yes") == 1
        assert "leave out --name" in capsys.readouterr().err
        assert stew(lab, "download", "zenodo:88", "--missing", "--dry-run") == 0
        out = capsys.readouterr().out
        assert "c.zip" in out and "a.zip" not in out.split("missing:")[1]
        assert not (Path(d["path"]) / "c.zip").exists()
        assert stew(lab, "download", "zenodo:88", "--missing", "--yes") == 0
        assert "Added 1 file(s)" in capsys.readouterr().out
        assert stew(lab, "download", "zenodo:88", "--missing", "--yes") == 0
        assert "Nothing to add" in capsys.readouterr().out
    finally:
        server.httpd.shutdown()
    cat = Catalog.open(lab)
    after = cat.find("zenodo-88")
    assert after["id"] == d["id"] and after["path"] == d["path"] and after["description"] == "Tiles"
    assert (after["n_files"], after["size_bytes"]) == (3, 60)
    assert (Path(d["path"]) / "c.zip").read_bytes() == b"c" * 30
    last = cat.downloads()[0]
    assert last["status"] == "done" and last["dataset_id"] == d["id"]
    assert last["total_bytes"] == 60  # the folder's files, so progress isn't 100% from the start


def test_download_missing_refuses_uncatalogued_or_changed_files(lab, tmp_path, monkeypatch, capsys):
    monkeypatch.setitem(sources.RESOLVERS, "zenodo", lambda src: sources.Remote(
        files=[sources.RemoteFile("a.zip", "https://x/a.zip", 99)]))
    assert stew(lab, "download", "zenodo:89", "--missing", "--yes") == 1
    assert "isn't catalogued yet" in capsys.readouterr().err
    folder = touch(tmp_path / "data/public/zenodo/89/a.zip", b"short").parent
    Catalog.open(lab).insert({"name": "z89", "path": str(folder), "source": "zenodo:89", "source_key": "zenodo:89",
                              "modalities": []})
    assert stew(lab, "download", "zenodo:89", "--missing", "--yes") == 1
    assert "differ in size" in capsys.readouterr().err
