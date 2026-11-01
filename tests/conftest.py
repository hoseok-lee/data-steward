import pytest

from data_steward import config


@pytest.fixture(autouse=True)
def isolated_config(tmp_path_factory, monkeypatch):
    """Keep tests away from the machine's real catalogues and configs."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path_factory.mktemp("xdg")))
    for var in ("STEW_ROOT", "DATA_STEWARD_CATALOG", "STEW_DOWNLOAD_DIR"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(config, "SYSTEM_CONFIG", tmp_path_factory.mktemp("etc") / "config.toml")
    monkeypatch.setattr(config, "site_config", lambda: None)


@pytest.fixture(autouse=True)
def fresh_settings():
    from data_steward import cli
    cli._SETTINGS.clear()
    yield
    cli._SETTINGS.clear()
