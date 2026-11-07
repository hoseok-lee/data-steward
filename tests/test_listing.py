"""Tests for downloading folders: BioStudies/BioImage Archive studies whose files sit
in subfolders, and plain URLs that point at a web server's folder listing."""

from pathlib import Path

import pytest

from data_steward import cli, sources
from data_steward.db import Catalog

from test_stew import Server, touch

BASE = "https://ftp.example.org/fire/S-BIAD/1/S-BIAD1/Files/"


def apache(path: str, *entries: str) -> str:
    """An Apache "Index of" page, with the sort links and parent link it always has."""
    rows = "".join(f'<tr><td><a href="{e}">{e}</a></td><td align="right"> 1K</td></tr>\n' for e in entries)
    return (f"<html><head><title>Index of {path}</title></head><body><h1>Index of {path}</h1><table>"
            '<tr><th><a href="?C=N;O=D">Name</a></th><th><a href="?C=S;O=A">Size</a></th></tr>'
            f'<tr><td><a href="{path.rsplit("/", 2)[0]}/">Parent Directory</a></td></tr>\n{rows}'
            '<a href="https://elsewhere.org/x.ims">elsewhere</a><a href="/">home</a></table></body></html>')


LISTINGS = {
    BASE: apache("/fire/S-BIAD/1/S-BIAD1/Files", "Tumor/", "Tumor_metadata.json", "Tumor_metadata.tsv"),
    BASE + "Tumor/": apache("/fire/S-BIAD/1/S-BIAD1/Files/Tumor", "DMG/", "readme.txt"),
    BASE + "Tumor/DMG/": apache("/fire/S-BIAD/1/S-BIAD1/Files/Tumor/DMG", "a%201.ims", "b.ims"),
    BASE + "Tumor/extra/": apache("/fire/S-BIAD/1/S-BIAD1/Files/Tumor/extra", "c.ims"),
}


@pytest.fixture
def ebi(monkeypatch):
    """BioStudies' API and file server, answering from the dictionaries above."""
    api = {
        "studies/S-BIAD1/info": {"httpLink": BASE.removesuffix("Files/")},
        "studies/S-BIAD1": {"attributes": [{"name": "Title", "value": "Tumour profiler"}], "section": {
            "attributes": [{"name": "Study type", "value": "imaging"}],
            "subsections": [{"type": "Study Component",
                             "attributes": [{"name": "File List", "value": "Tumor_metadata.json"}]}]}},
        "files": [  # each file once per section that uses it, one folder, served two rows a page
            {"path": "Tumor/DMG/a 1.ims", "Size": 5_000_000_000, "isDirectory": "false", "Section": "A"},
            {"path": "Tumor/DMG/a 1.ims", "Size": 5_000_000_000, "isDirectory": "false", "Section": "B"},
            {"path": "Tumor/DMG/b.ims", "Size": 7_000_000_000, "isDirectory": "false", "Section": "A"},
            {"path": "Tumor/extra", "isDirectory": "true", "type": "directory", "Section": "A"},
        ],
    }
    calls = []

    def fake_json(url):
        calls.append(url)
        path = url.split("/api/v1/", 1)[1]
        if path.startswith("files/"):
            q = dict(p.split("=") for p in path.split("?", 1)[1].split("&"))
            start, length = int(q["start"]), int(q["length"])
            rows = api["files"]
            return {"recordsTotal": len(rows), "data": rows[start:start + length]}
        return api[path]

    def fake_text(url):
        calls.append(url)
        if url not in LISTINGS:
            raise sources.SourceError(f"{url} answered 404 Not Found")
        return LISTINGS[url]

    monkeypatch.setattr(sources, "BIOSTUDIES_PAGE", 2)
    monkeypatch.setattr(sources, "get_json", fake_json)
    monkeypatch.setattr(sources, "get_text", fake_text)
    monkeypatch.setattr(sources, "head", lambda url: (100, None))
    return api, calls


def test_biostudies_lists_files_in_subfolders(ebi):
    api, calls = ebi
    r = sources.resolve(sources.parse("S-BIAD1"))
    assert {f.name: f.size for f in r.files} == {
        "Tumor/DMG/a 1.ims": 5_000_000_000, "Tumor/DMG/b.ims": 7_000_000_000, "Tumor/extra/c.ims": 100,
        "Tumor_metadata.json": 100, "Tumor_metadata.tsv": 100}
    urls = {f.name: f.url for f in r.files}
    assert urls["Tumor/DMG/a 1.ims"] == BASE + "Tumor/DMG/a%201.ims"
    assert urls["Tumor/extra/c.ims"] == BASE + "Tumor/extra/c.ims"
    assert sum("/files/S-BIAD1?" in c for c in calls) == 2  # paged through the whole list
    assert BASE + "Tumor/" not in calls  # the API already listed that folder


def test_biostudies_crawls_the_folders_without_a_file_list(ebi):
    api, _ = ebi
    api["files"] = []
    r = sources.resolve(sources.parse("S-BIAD1"))
    assert sorted(f.name for f in r.files) == ["Tumor/DMG/a 1.ims", "Tumor/DMG/b.ims", "Tumor/readme.txt",
                                               "Tumor_metadata.json", "Tumor_metadata.tsv"]
    assert all(f.url.startswith(BASE) and f.size == 100 for f in r.files)


def test_folder_listing_stays_below_its_folder(ebi):
    files = sources.crawl_listing(BASE + "Tumor/")
    assert sorted(f.name for f in files) == ["DMG/a 1.ims", "DMG/b.ims", "readme.txt"]


# --- plain URLs of a folder listing, against a real local server ---------------------------------------

@pytest.fixture
def served(tmp_path):
    root = tmp_path / "served"
    touch(root / "run/a.txt", b"a" * 10)
    touch(root / "run/sub/b.txt", b"b" * 20)
    touch(root / "run/sub/deeper/c.txt", b"c" * 30)
    touch(root / "other.txt", b"x")  # beside the folder: must not be fetched
    (root / "page.html").write_text("<html><title>Hello</title><a href='other.txt'>x</a></html>")
    server = Server(root)
    yield server
    server.httpd.shutdown()


def test_url_of_a_folder_listing_is_a_folder(served):
    for url in (f"{served.url}/run/", f"{served.url}/run"):  # with or without the slash
        r = sources.resolve(sources.parse(url))
        assert {f.name: f.size for f in r.files} == {"a.txt": 10, "sub/b.txt": 20, "sub/deeper/c.txt": 30}
    r = sources.resolve(sources.parse(f"{served.url}/page.html"))  # an ordinary web page is still one file
    assert [f.name for f in r.files] == ["page.html"]
    r = sources.resolve(sources.parse(f"{served.url}/run/a.txt"))
    assert [(f.name, f.size) for f in r.files] == [("a.txt", 10)]


def test_download_url_of_a_folder_keeps_its_layout(served, tmp_path, capsys):
    catalog = tmp_path / "data/.data-steward/catalog.db"
    Catalog.create(catalog)
    url = f"{served.url}/run/"
    assert cli.main(["--catalog", str(catalog), "download", url, "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "sub/deeper/c.txt" in out and "3 file(s), 60 B" in out
    assert cli.main(["--catalog", str(catalog), "download", url, "--yes"]) == 0
    d = Catalog.open(catalog).find("run")
    root = Path(d["path"])
    assert root == tmp_path / "data/public/url/run"
    assert sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()) == [
        "a.txt", "sub/b.txt", "sub/deeper/c.txt"]
    assert (root / "sub/deeper/c.txt").read_bytes() == b"c" * 30
    assert not any(p.name.endswith(".html") for p in root.rglob("*"))
