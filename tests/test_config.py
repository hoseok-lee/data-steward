"""Finding the catalogue, and the per-folder config."""

from pathlib import Path

from data_steward import cli, config
from data_steward.db import Catalog


def stew(*args):
    return cli.main(list(args))


def test_init_creates_catalogue_and_config(tmp_path, capsys):
    root2 = tmp_path / "shared2"
    assert stew("init", str(root2)) == 0
    assert (root2 / ".data-steward/catalog.db").is_file() and (root2 / ".data-steward/config.toml").is_file()
    s = config.load(root2 / ".data-steward/catalog.db")
    assert s.move_on_add and s.layout == config.DEFAULT_LAYOUT and s.root == root2.resolve()
    assert stew("init", str(root2)) == 1  # already there


def test_discovery_order(tmp_path, monkeypatch, capsys):
    a, b = tmp_path / "a", tmp_path / "b"
    for r in (a, b):
        assert stew("init", str(r)) == 0
    sub = a / "projects/x"
    sub.mkdir(parents=True)
    monkeypatch.chdir(sub)
    assert config.locate()[0] == config.catalog_in(a.resolve())          # walked up from cwd
    monkeypatch.setenv("STEW_ROOT", str(b))
    assert config.locate()[0] == config.catalog_in(b)                     # env beats cwd
    assert config.locate(a)[0] == config.catalog_in(a)                    # --catalog beats env
    monkeypatch.delenv("STEW_ROOT")
    monkeypatch.chdir(tmp_path)
    assert config.locate() is None
    capsys.readouterr()
    assert stew("ls") == 1 and "no catalogue found" in capsys.readouterr().err
    assert stew("init", str(b / "c"), "--set-default") == 0               # user config
    assert config.locate()[0] == config.catalog_in((b / "c").resolve())
    capsys.readouterr()
    assert stew("config", "--get", "root") == 0
    assert capsys.readouterr().out.strip() == str((b / "c").resolve())


def test_system_config_and_site_config(tmp_path, monkeypatch):
    root = tmp_path / "srv"
    assert stew("init", str(root)) == 0
    monkeypatch.chdir(tmp_path)
    config.SYSTEM_CONFIG.write_text(f'root = "{root}"\n')
    assert config.locate() == (config.catalog_in(root), str(config.SYSTEM_CONFIG))
    config.SYSTEM_CONFIG.unlink()
    site = tmp_path / "site.toml"
    site.write_text(f'catalog = "{config.catalog_in(root)}"\n')
    monkeypatch.setattr(config, "site_config", lambda: site)
    assert config.locate()[0] == config.catalog_in(root)


def test_layout_and_move_on_add_from_config(tmp_path, monkeypatch):
    root = tmp_path / "data"
    assert stew("init", str(root)) == 0
    cfg = root / ".data-steward/config.toml"
    cfg.write_text(cfg.read_text().replace('project = "projects"', 'project = "lab/{name}-proj"')
                   .replace('public = "public/{kind}/{accession}"', 'public = "external/{accession}"'))
    monkeypatch.setenv("STEW_ROOT", str(root))
    mine = tmp_path / "me/pilot.h5ad"
    mine.parent.mkdir()
    mine.write_bytes(b"x")
    assert stew("add", str(mine), "--yes") == 0
    assert (root / "lab/pilot.h5ad-proj/pilot.h5ad").is_file() and mine.is_symlink()
    other = tmp_path / "me/GSE7"
    other.mkdir()
    (other / "a.h5ad").write_bytes(b"y")
    assert stew("add", str(other), "--whole", "--source", "GSE7", "--yes") == 0
    assert (root / "external/GSE7/a.h5ad").is_file()
    cli._SETTINGS.clear()
    cfg.write_text(cfg.read_text().replace("move_on_add = true", "move_on_add = false"))
    stay = tmp_path / "me/stay.h5ad"
    stay.write_bytes(b"z")
    assert stew("add", str(stay), "--yes") == 0
    assert not stay.is_symlink() and Catalog.open(config.catalog_in(root)).find("stay")["path"] == str(stay)


def test_init_in_place(tmp_path, monkeypatch):
    assert stew("init", str(tmp_path / "d"), "--in-place") == 0
    assert config.load(config.catalog_in(tmp_path / "d")).move_on_add is False


def test_bad_config_is_reported(tmp_path, monkeypatch, capsys):
    root = tmp_path / "d"
    assert stew("init", str(root)) == 0
    (root / ".data-steward/config.toml").write_text('[layout]\nmisc = "x"\n')
    cli._SETTINGS.clear()
    monkeypatch.setenv("STEW_ROOT", str(root))
    f = tmp_path / "x.h5ad"
    f.write_bytes(b"x")
    capsys.readouterr()
    assert stew("add", str(f), "--yes") == 1
    assert "unknown [layout] keys misc" in capsys.readouterr().err
    assert Path(f).is_file() and not f.is_symlink()
