"""Where the catalogue is, and how its data folder is laid out.

A catalogue lives in `<root>/.data-steward/`: `catalog.db` and `config.toml`.
stew finds it, in order, from:

  1. --catalog (a catalogue file, or a data root folder)
  2. $DATA_STEWARD_CATALOG (a file) or $STEW_ROOT (a folder)
  3. the current folder or one above it holding .data-steward/  (like git)
  4. ~/.config/stew/config.toml          (`root = "/path"`)
  5. /etc/stew/config.toml               (same; for everyone on a server)
  6. a site.toml built into the package  (same; for site-specific builds)
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover
    tomllib = None

CATALOG_DIR = ".data-steward"
CATALOG_FILE = "catalog.db"
CONFIG_FILE = "config.toml"
SYSTEM_CONFIG = Path("/etc/stew/config.toml")

# Folders (relative to the root) that datasets are moved into; stew adds the
# dataset's own file or folder name.
DEFAULT_LAYOUT = {
    "public": "public/{kind}/{accession}",   # data with a source (GEO, Zenodo, ...)
    "reference": "references/{ref_kind}",    # genomes, indices, annotations (10x or annotations)
    "project": "projects",                   # everything else
}

CONFIG_TEMPLATE = '''\
# data-steward configuration for this data folder. All keys are optional.

# Move datasets added from elsewhere into this folder (leaving a symlink behind).
# false: register them where they are.
move_on_add = {move_on_add}

# Unix group for folders stew creates, and their mode. 3775 = group-writable,
# new files inherit the group, and only a dataset's owner may delete or rename it.
{group_line}
folder_mode = "3775"

# Where datasets go, relative to this folder. stew appends the dataset's name.
#   {{kind}} {{accession}}  e.g. geo / GSE123456, zenodo / 7216361
#   {{ref_kind}}          10x or annotations
[layout]
public = "{public}"
reference = "{reference}"
project = "{project}"
'''


class ConfigError(Exception):
    pass


@dataclass
class Settings:
    root: Path
    catalog: Path
    found: str = ""
    move_on_add: bool = True
    group: str | None = None
    folder_mode: int = 0o3775
    layout: dict = field(default_factory=lambda: dict(DEFAULT_LAYOUT))


def _read_toml(path: Path) -> dict:
    if tomllib is None:
        raise ConfigError("reading config files needs Python 3.11 or newer")
    try:
        with path.open("rb") as f:
            return tomllib.load(f)
    except (OSError, tomllib.TOMLDecodeError) as e:
        raise ConfigError(f"cannot read {path}: {e}")


def catalog_in(root: Path) -> Path:
    return root / CATALOG_DIR / CATALOG_FILE


def _as_catalog(p: Path) -> Path:
    p = p.expanduser()
    return catalog_in(p) if p.is_dir() and not p.name.endswith(".db") else p


def _pointer(path: Path) -> Path | None:
    """A catalogue named by a pointer config file (`root =` or `catalog =`)."""
    if not path.is_file():
        return None
    d = _read_toml(path)
    if d.get("catalog"):
        return Path(d["catalog"]).expanduser()
    if d.get("root"):
        return catalog_in(Path(d["root"]).expanduser())
    return None


def user_config_path() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or os.path.join(os.path.expanduser("~"), ".config")
    return Path(base) / "stew" / "config.toml"


def site_config() -> Path | None:
    p = Path(__file__).with_name("site.toml")
    return p if p.is_file() else None


def locate(explicit: Path | None = None, cwd: Path | None = None) -> tuple[Path, str] | None:
    """(catalogue file, how it was found), or None."""
    if explicit:
        return _as_catalog(explicit), "--catalog"
    if os.environ.get("DATA_STEWARD_CATALOG"):
        return Path(os.environ["DATA_STEWARD_CATALOG"]).expanduser(), "$DATA_STEWARD_CATALOG"
    if os.environ.get("STEW_ROOT"):
        return catalog_in(Path(os.environ["STEW_ROOT"]).expanduser()), "$STEW_ROOT"
    here = (cwd or Path.cwd()).resolve()
    for d in (here, *here.parents):
        if catalog_in(d).is_file():
            return catalog_in(d), f"the current folder ({d})"
    for path, label in ((user_config_path(), None), (SYSTEM_CONFIG, None), (site_config(), "built-in site config")):
        if path is not None:
            c = _pointer(path)
            if c:
                return c, label or str(path)
    return None


NOT_FOUND = ("no catalogue found. Create one with `stew init <data folder>`, run stew from inside a data "
             "folder, or point to one with --catalog, $STEW_ROOT, ~/.config/stew/config.toml or "
             "/etc/stew/config.toml (`root = \"/path/to/data\"`)")


def load(catalog: Path, found: str = "") -> Settings:
    """Settings for the data folder holding `catalog`."""
    catalog = catalog.expanduser()
    root = catalog.resolve().parent.parent if catalog.exists() else catalog.parent.parent
    s = Settings(root=root, catalog=catalog, found=found)
    path = catalog.parent / CONFIG_FILE
    if path.is_file():
        d = _read_toml(path)
        s.move_on_add = bool(d.get("move_on_add", True))
        s.group = d.get("group") or None
        try:
            s.folder_mode = int(str(d.get("folder_mode", "3775")), 8)
        except ValueError:
            raise ConfigError(f"{path}: folder_mode must be an octal string like \"3775\"")
        layout = d.get("layout") or {}
        unknown = set(layout) - set(DEFAULT_LAYOUT)
        if unknown:
            raise ConfigError(f"{path}: unknown [layout] keys {', '.join(sorted(unknown))} "
                              f"(known: {', '.join(DEFAULT_LAYOUT)})")
        s.layout.update({k: str(v).strip("/") for k, v in layout.items()})
    return s


def write_config(catalog: Path, group: str | None = None, move_on_add: bool = True,
                 layout: dict | None = None) -> Path:
    lay = {**DEFAULT_LAYOUT, **(layout or {})}
    path = catalog.parent / CONFIG_FILE
    path.write_text(CONFIG_TEMPLATE.format(
        move_on_add="true" if move_on_add else "false",
        group_line=f'group = "{group}"' if group else '# group = "lab"',
        **lay))
    return path


def write_user_config(root: Path) -> Path:
    path = user_config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f'# Default data-steward catalogue for this user\nroot = "{root}"\n')
    return path
