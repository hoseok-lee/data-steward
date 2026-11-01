import functools
import hashlib
import http.server
import json
import sqlite3
import threading
from pathlib import Path

import pytest

import data_steward
from data_steward import cli, detect, sources, util, vocab
from data_steward.db import SCHEMA_V1, Catalog, Changes, Filter


def touch(p: Path, data: bytes = b"x") -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data)
    return p


def test_version_matches_pyproject():
    toml = (Path(__file__).parent.parent / "pyproject.toml").read_text()
    assert f'version = "{data_steward.__version__}"' in toml


# --- vocabulary and detection -------------------------------------------------------------

def test_vocab_normalises_aliases():
    assert vocab.organ("Lungs") == "lung"
    assert vocab.organ("PDAC") == "pancreas"
    assert vocab.organ("lymph node") == "lymph-node"
    assert vocab.organ("Gallbladder") == "gallbladder"
    assert vocab.species("Homo sapiens") == "human"
    assert vocab.modality("scRNA") == "rna"
    with pytest.raises(vocab.VocabError):
        vocab.modality("transcripts-ish")
    assert vocab.platform("MERFISH") == "merscope"


def test_vocab_infers_from_paths():
    toks = vocab.tokens("/data/spatial/Visium_HD/human_Lung5_Rep1/outs")
    assert vocab.infer_organ(toks) == "lung"
    assert vocab.infer_species(toks) == "human"
    assert vocab.infer_platform(toks) == "visium-hd"
    assert vocab.infer_modalities(toks) == {"spatial"}
    assert vocab.infer_modalities(vocab.tokens("cosmx/pdac_slide1")) == set()


def test_util():
    assert util.slug("Lung5_Rep1 (FFPE)") == "lung5-rep1-ffpe"
    assert util.human_size(512) == "512 B"
    assert util.human_size(3 * 1024**3) == "3.0 GiB"
    assert util.levenshtein("hcla-core", "hlca-core") == 2


def test_scan_finds_bundles_and_files(tmp_path):
    r = tmp_path
    (r / "xen/run1").mkdir(parents=True)
    (r / "xen/run1/experiment.xenium").write_text('{"panel_organism":"Human","panel_tissue_type":"Lung"}')
    touch(r / "mouse_brain/outs/spatial/tissue_positions.csv")
    touch(r / "mouse_brain/outs/filtered_feature_bc_matrix.h5")
    touch(r / "atlas/liver_atlas.h5ad")
    touch(r / "reads/S1_R1_001.fastq.gz")
    touch(r / "reads/S1_R2_001.fastq.gz")
    touch(r / "pdac_slide1/Pancreas_exprMat_file.csv")
    touch(r / "half/thing.h5ad.part")

    found = {c.name: c for c in detect.scan(r).found}
    assert sorted(found) == ["liver-atlas", "mouse-brain", "pdac-slide1", "reads", "run1"]
    x = found["run1"]
    assert (x.platform, x.organ, x.species) == ("xenium", "lung", "human")
    v = found["mouse-brain"]
    assert (v.platform, v.organ, v.species) == ("visium", "brain", "mouse")
    assert {"rna", "spatial"} <= v.modalities
    assert found["liver-atlas"].format == "anndata" and "rna" in found["liver-atlas"].modalities
    assert found["reads"].format == "fastq"
    assert (found["pdac-slide1"].platform, found["pdac-slide1"].organ) == ("cosmx", "pancreas")


# --- catalogue ------------------------------------------------------------------------------------

def ds(name, organ, mods, **kw):
    return {"name": name, "path": f"/data/{name}", "organ": organ, "modalities": mods, **kw}


def test_catalog_filters_and_links(tmp_path):
    cat = Catalog.create(tmp_path / "c/catalog.db")
    cat.insert(ds("a", "lung", ["rna", "spatial"]))
    cat.insert(ds("b", "lung", ["rna"], source="GSE1", source_key="geo:GSE1"))
    cat.insert(ds("c", "liver", ["rna", "spatial"]))

    def names(**kw):
        return [d["name"] for d in cat.list(Filter(**kw))]

    assert names(organs=["lung"]) == ["a", "b"]
    assert names(organs=["lung"], modalities=["rna", "spatial"]) == ["a"]
    assert len(names(under="/data")) == 3
    assert names(under="/dat") == []
    assert names(source_key="geo:GSE1") == ["b"]

    a = cat.find("a")
    assert a["modalities"] == ["rna", "spatial"]
    with pytest.raises(Exception, match="Did you mean"):
        cat.find("aa")
    assert cat.unique_name("a") == "a-2"

    cat.add_link(a, Path("/home/x/a"))
    link = cat.active_link(Path("/home/x/a"))
    cat.release_link(link)
    assert cat.links() == [] and len(cat.links(released=True)) == 1

    cat.update(a, Changes(add_tags=["ffpe"], fields=[("source", "zenodo:9"), ("source_key", "zenodo:9")]))
    assert cat.find("a")["source"] == "zenodo:9"
    cat.delete(cat.find("b"))
    assert cat.summary()["datasets"] == 2


def test_upgrades_a_version_1_catalogue(tmp_path):
    """Catalogues made by the Rust version (schema 1) are upgraded in place."""
    path = tmp_path / "old.db"
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA_V1)
    conn.execute("INSERT INTO meta VALUES ('schema_version', '1')")
    conn.execute("INSERT INTO datasets (name, path, added_by, added_at, updated_at) "
                 "VALUES ('old', '/x', 'me', 'then', 'then')")
    conn.commit()
    conn.close()
    cat = Catalog.open(path)
    assert cat.find("old")["source"] is None
    assert cat.conn.execute("SELECT value FROM meta").fetchone()[0] == "4"
    assert cat.downloads() == []


# --- sources ------------------------------------------------------------------------------------

@pytest.mark.parametrize("text,key", [
    ("GSE274673", "geo:GSE274673"),
    ("gsm8454446", "geo:GSM8454446"),
    ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274673", "geo:GSE274673"),
    ("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274673/suppl/", "geo:GSE274673"),
    ("zenodo:5979761", "zenodo:5979761"),
    ("10.5281/zenodo.5979761", "zenodo:5979761"),
    ("https://doi.org/10.5281/zenodo.5979761", "zenodo:5979761"),
    ("https://zenodo.org/records/5979761", "zenodo:5979761"),
    ("HTTPS://Example.org:443/a/b.h5ad/#frag", "url:https://example.org/a/b.h5ad"),
])
def test_parse_sources(text, key):
    assert sources.parse(text).key == key


def test_file_urls_know_their_collection():
    s = sources.parse("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274673/suppl/GSE274673_RAW.tar")
    assert s.kind == "url" and s.related == ["geo:GSE274673"]
    assert s.key == "url:https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274673/suppl/GSE274673_RAW.tar"
    s = sources.parse("https://zenodo.org/records/123/files/a.zip")
    assert s.related == ["zenodo:123"]
    with pytest.raises(sources.SourceError):
        sources.parse("my dataset")


# --- download --------------------------------------------------------------------------------------

class Server:
    """A local HTTP server that counts requests."""

    def __init__(self, root: Path):
        self.requests = []
        outer = self

        class Handler(http.server.SimpleHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def send_head(self):
                outer.requests.append(self.path)
                return super().send_head()

        self.httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Handler, directory=root))
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}"


@pytest.fixture
def env(tmp_path):
    served = tmp_path / "served"
    touch(served / "lung_atlas.h5ad", b"a" * 1000)
    server = Server(served)
    catalog = tmp_path / "data/.data-steward/catalog.db"
    Catalog.create(catalog)
    yield server, catalog, tmp_path
    server.httpd.shutdown()


def stew(catalog, *args):
    return cli.main(["--catalog", str(catalog), *args])


def test_download_url_then_refuses_before_any_request(env, capsys):
    server, catalog, root = env
    url = f"{server.url}/lung_atlas.h5ad"
    assert stew(catalog, "download", url, "--yes", "--species", "human") == 0
    d = Catalog.open(catalog).find("lung-atlas")
    assert d["source"] == url and d["format"] == "anndata" and d["organ"] == "lung"
    assert Path(d["path"]).read_bytes() == b"a" * 1000
    assert Path(d["path"]).parent == root / "data/public/url/lung-atlas-h5ad"
    assert not (Path(d["path"]).parent / ".stew-download").exists()

    server.requests.clear()
    capsys.readouterr()
    assert stew(catalog, "download", url, "--yes") == 1
    assert "already downloaded" in capsys.readouterr().err
    assert server.requests == []  # the check ran before anything was fetched
    assert stew(catalog, "ls", "--source", url, "--paths") == 0
    assert capsys.readouterr().out.strip() == d["path"]


def test_download_refuses_when_a_dataset_mentions_the_accession(env, capsys):
    _, catalog, root = env
    folder = root / "data/raw/GSE999_pancreas"
    touch(folder / "x.fastq.gz")
    assert stew(catalog, "add", str(folder), "--yes") == 0
    assert stew(catalog, "download", "GSE999", "--yes") == 1
    err = capsys.readouterr().err
    assert "mention GSE999" in err and "--force" in err


def test_download_refuses_while_someone_else_is_downloading(env, capsys):
    _, catalog, _ = env
    cat = Catalog.open(catalog)
    cat.conn.execute("INSERT INTO downloads (source_key, source, user, host, pid, dest, status, started_at) "
                     "VALUES ('zenodo:42', 'zenodo:42', 'alice', 'other-host', 1, '/x', 'running', 'now')")
    assert stew(catalog, "download", "zenodo:42", "--yes") == 1
    assert "alice has been downloading" in capsys.readouterr().err


def test_download_refuses_a_file_from_a_series_already_here(env, capsys):
    _, catalog, _ = env
    Catalog.open(catalog).insert({"name": "gse5", "path": "/d/gse5", "source": "GSE5", "source_key": "geo:GSE5",
                                  "modalities": []})
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE5/suppl/GSE5_counts.h5"
    assert stew(catalog, "download", url, "--yes") == 1
    assert "same as) GSE5" in capsys.readouterr().err


def test_fetch_checks_md5_and_resumes(tmp_path, env):
    server, _, _ = env
    data = b"a" * 1000
    from data_steward.fetch import FetchError, fetch
    f = sources.RemoteFile("lung_atlas.h5ad", f"{server.url}/lung_atlas.h5ad", 1000,
                           hashlib.md5(data).hexdigest())
    out = tmp_path / "out"
    out.mkdir()
    (out / "lung_atlas.h5ad.part").write_bytes(b"a" * 10)  # left by an interrupted run
    assert fetch(f, out).read_bytes() == data
    nested = sources.RemoteFile("sub/dir/copy.h5ad", f.url, 1000)  # Zenodo keys may contain '/'
    assert fetch(nested, out) == (out / "sub/dir/copy.h5ad").resolve()
    with pytest.raises(FetchError, match="outside"):
        fetch(sources.RemoteFile("../escape.h5ad", f.url, 1000), out)
    bad = sources.RemoteFile("lung_atlas.h5ad", f.url, 1000, "0" * 32)
    with pytest.raises(FetchError, match="md5"):
        fetch(bad, tmp_path)


def test_bare_stew_prints_help(capsys):
    assert cli.main([]) == 0
    out = capsys.readouterr().out
    assert "Find datasets:" in out and "download" in out


def test_json_listing(env, capsys):
    _, catalog, root = env
    touch(root / "data/x/mouse_kidney.h5ad")
    stew(catalog, "scan", str(root / "data/x"), "--yes")
    capsys.readouterr()
    assert stew(catalog, "ls", "--species", "mouse", "--json") == 0
    rows = json.loads(capsys.readouterr().out)
    assert [r["organ"] for r in rows] == ["kidney"]


def test_output_into_head_exits_cleanly(tmp_path):
    import subprocess
    import sys
    catalog = tmp_path / "c.db"
    Catalog.create(catalog)
    cat = Catalog.open(catalog)
    for i in range(3000):
        cat.insert({"name": f"d{i}", "path": f"/x/{i}", "modalities": []})
    src = str(Path(cli.__file__).parent.parent)
    import os
    r = subprocess.run(f"{sys.executable} -m data_steward --catalog {catalog} ls | head -2", shell=True,
                       capture_output=True, text=True, env={**os.environ, "PYTHONPATH": src})
    assert "Broken pipe" not in r.stderr and "Exception ignored" not in r.stderr


def test_only_stewards_change_the_catalogue(tmp_path, capsys, monkeypatch):
    catalog = tmp_path / "data/.data-steward/catalog.db"
    assert cli.main(["--catalog", str(catalog), "init", "--steward", "alice"]) == 0
    folder = tmp_path / "data/x"
    folder.mkdir()
    (folder / "a.h5ad").write_bytes(b"x")
    cat = Catalog.open(catalog)
    cat.insert({"name": "a", "path": str(folder / "a.h5ad"), "modalities": []})
    capsys.readouterr()
    # not a steward (the test runs as whoever we are, not alice)
    for args in (["scan", str(folder), "--yes"], ["set", "a", "--organ", "lung"], ["rm", "a", "--yes"],
                 ["download", "zenodo:1", "--yes"], ["stewards", "--add", "me"]):
        assert cli.main(["--catalog", str(catalog), *args]) == 1, args
        assert "stewards (alice)" in capsys.readouterr().err
    # reading and linking still work
    assert cli.main(["--catalog", str(catalog), "ls"]) == 0
    assert cli.main(["--catalog", str(catalog), "get", "a", "--to", str(tmp_path)]) == 0
    # a steward can do it
    import data_steward.util as u
    monkeypatch.setattr(u, "current_user", lambda: "alice")
    assert cli.main(["--catalog", str(catalog), "set", "a", "--organ", "lung"]) == 0
    assert cli.main(["--catalog", str(catalog), "stewards", "--add", "bob"]) == 0
    assert "alice, bob" in capsys.readouterr().out


def test_check_and_backup(tmp_path, capsys):
    catalog = tmp_path / "data/.data-steward/catalog.db"
    Catalog.create(catalog)
    f = tmp_path / "data/x.h5ad"
    f.write_bytes(b"abc")
    assert cli.main(["--catalog", str(catalog), "add", str(f), "--yes"]) == 0
    assert cli.main(["--catalog", str(catalog), "checksum", "--all"]) == 0
    capsys.readouterr()
    report = tmp_path / "r.txt"
    assert cli.main(["--catalog", str(catalog), "check", "--report", str(report)]) == 0
    assert "RESULT: OK" in report.read_text()
    f.write_bytes(b"abd")
    assert cli.main(["--catalog", str(catalog), "check"]) == 1
    assert "1 changed" in capsys.readouterr().out
    for _ in range(3):
        assert cli.main(["--catalog", str(catalog), "backup", "--keep", "2"]) == 0
        import time as _t
        _t.sleep(1.1)
    backups = sorted((tmp_path / "data/.data-steward/backups").glob("catalog-*.db"))
    assert len(backups) == 2
    assert Catalog.open(backups[-1]).find("x")["name"] == "x"


def test_storage_report(tmp_path, capsys):
    import os as _os
    catalog = tmp_path / "data/.data-steward/catalog.db"
    Catalog.create(catalog)
    for sub, size in (("projects/a", 3000), ("public/geo/GSE1", 1000)):
        d = tmp_path / "data" / sub
        d.mkdir(parents=True)
        (d / "x.h5ad").write_bytes(b"x" * size)
        _os.utime(d / "x.h5ad", (200000, 200000))  # 1970-01-03, in any time zone
        _os.utime(d, (200000, 200000))
        assert cli.main(["--catalog", str(catalog), "add", str(d), "--yes"]) == 0
    capsys.readouterr()
    assert cli.main(["--catalog", str(catalog), "report", "--top", "5"]) == 0
    out = capsys.readouterr().out
    assert "By area" in out and "projects" in out and "public/geo" in out
    assert "75.0%" in out  # projects holds 3000 of 4000 bytes
    assert out.count("1970-01-03") == 2  # both untouched since 1970, never linked
