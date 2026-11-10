"""The `stew` command."""

from __future__ import annotations

import argparse
import csv
import re
import json
import shlex
import shutil
import subprocess
import time
from dataclasses import replace
import os
import sqlite3
import sys
from pathlib import Path

from . import __version__, checksum, config, contents, detect, provenance, sources, util, vocab
from .db import Catalog, CatalogError, Changes, Filter
from .detect import Candidate
from .fetch import FetchError, check_space, fetch_many, pending_bytes
from .util import Presence, human_size


class UserError(Exception):
    pass


def positive(text: str) -> int:
    try:
        n = int(text)
    except ValueError:
        n = 0
    if n < 1:
        raise argparse.ArgumentTypeError(f"must be a whole number of at least 1, not {text!r}")
    return n


# --- help text --------------------------------------------------------------------------

TOP_HELP = f"""\
stew {__version__}: a catalogue of your lab's datasets

Usage: stew [--catalog FILE] [-y] <command> [options]

Find datasets:
  ls          List datasets, filtered by organ, modality, species, source, ... [alias: list]
  search      Find datasets by words in their name, description, study, source or path
  saved       Your saved `ls` queries (save one with `stew ls ... --save NAME`)
  info        Everything known about one dataset, and whether you can read it
  summary     Counts by modality, organ, species and platform
  report      Storage report: space by area, owner, modality; largest and unused datasets
  vocab       The organs, modalities, species and platforms stew knows, with aliases
  serve       A read-only web page of the catalogue

Use datasets:
  get         Put a symlink to a dataset into one of your folders
  links       Links you have made with `get` (--all for everyone's)
  release     Remove a link made by `get`

Bring in data:
  download    Download from GEO, SRA, CELLxGENE, ArrayExpress, HCA, figshare, Zenodo or a URL,
              unless it is already here (--background to run it in tmux)
  downloads   Downloads in progress, finished and failed

Look after the catalogue:
  add         Register a file, a dataset folder, or every dataset inside a folder [alias: scan]
  set         Change a dataset's organ, species, modalities, source, tags, ...
  rm          Remove a dataset from the catalogue (files are never touched)
  inspect     Read cell/feature counts, organism, tissue and assay from inside the files
  relate      Record that one dataset comes from (or is part of) another; `unrelate` undoes it
  mv          Move a dataset or folder; updates the catalogue and the links handed out
  checksum    Record a checksum of every file; `verify` checks them later
  duplicates  Datasets and files stored more than once
  refresh     Re-check that datasets still exist; update sizes and owners
  stewards    Who may change the catalogue (everyone else can find and link)
  log         Recent changes, and who made them
  config      Which catalogue stew uses, how it found it, and its settings
  check       Health check: every dataset present, files unchanged since their checksum
  backup      Copy the catalogue to .data-steward/backups/ (keeps the newest 30)
  init        Set up a data folder: catalogue, config, optional default for you

Options:
  --catalog PATH  Catalogue file or data folder to use (default: $STEW_ROOT, the data folder
                  you are in, ~/.config/stew/config.toml, /etc/stew/config.toml; see `stew config`)
  -y, --yes       Don't ask questions; accept the suggested answers
  -h, --help      Show this help
  -V, --version   Show the version

Examples:
  stew ls                                   everything in the catalogue
  stew ls --organ lung --modality rna       lung datasets with RNA
  stew ls -m rna -m spatial --species human spatial transcriptomics from human tissue
  stew info hlca-core                       details for one dataset (name or id)
  stew get hlca-core --to ~/proj/data       symlink it into your folder
  stew search fibroblast pdac               datasets mentioning both words
  stew download GSE274673                   fetch a GEO series, unless it is already here
  stew add /data/me/project/data       register every dataset in a folder

A dataset can be named by its name, its id (from `stew ls`) or its path.
Run `stew help <command>` (or `stew <command> --help`) for a command's options.
"""

LS_EPILOG = """\
Repeating --organ, --species, --platform or --format matches ANY of the values;
repeating --modality or --tag requires ALL of them. Aliases are understood
(--organ lungs, --organ pdac, -m scrna, --species mm10); see `stew vocab`.

Examples:
  stew ls --organ lung --modality rna
  stew ls -m rna -m spatial                 spatial transcriptomics
  stew ls --organ lung --organ liver --species human
  stew ls --source GSE274673                datasets from one GEO series
  stew ls --platform xenium --long          include path, format, source and tags
  stew ls --owner alice --paths         just paths, e.g. for a script
  stew ls --under /data/shared/spatial --json
  stew ls -m spatial --species human --save human-spatial   remember this query
  stew ls -q human-spatial --csv > spatial.csv              run it again, as CSV"""

GET_EPILOG = """\
The link points at the shared copy; nothing is copied. You must already be able to
read the dataset; if you can't, stew tells you whom to ask. Links are recorded,
so `stew links` shows yours and `stew release <link>` removes one.

Examples:
  stew get hlca-core                        link into the current folder
  stew get 21 --to ~/proj/data --as lung1   by id, with a different link name"""

SET_EPILOG = """\
Examples:
  stew set hlca-core --species human --study HLCA --desc "Human Lung Cell Atlas core"
  stew set 7 --source https://cellxgene.cziscience.com/collections/6f6d381a
  stew set 7 --add-modality protein --tag atlas
  stew set 7 -m rna -m spatial              replace all modalities
  stew set 7 --organ ""                     clear a field"""

ADD_EPILOG = """\
What gets registered depends on what you point at:

  a file                                   that file, as one dataset
  a folder stew recognises as one dataset  the folder (Xenium, Visium, Cell Ranger, CosMx,
                                           MERSCOPE, .zarr, mtx folder, count tables,
                                           FASTQ folder, a folder with a .stew-dataset file)
  any other folder                         each dataset stew finds inside it, separately;
                                           loose scripts, logs and other files stay put.
                                           If it finds none, the folder as one dataset.

--whole registers a folder as one dataset whatever is inside; --each registers each folder
directly inside it as one dataset. Always fine to try --dry-run first.

Data lives in the shared data folder (the folder holding .data-steward/; see `stew config`).
Datasets added from anywhere else are moved there and a symlink is left at each old path,
so scripts and notebooks that use the old path keep working. Where they go:

  a GEO/SRA/Zenodo/... source        public/<kind>/<accession>/  (e.g. public/geo/GSE123456)
  a reference (genome, annotations)  references/10x/ or references/annotations/
  anything else                      projects/<name>/
(the default layout; change it in .data-steward/config.toml, or set move_on_add = false)

Pick the folder yourself with --to (e.g. --to projects/my_study), or keep the data where it
is with --in-place. The move is a rename, so it is instant, but only works on the same disk
as the data folder; data on another disk must be copied there by hand first.
(`stew scan` is another name for `stew add`.)

Examples:
  stew add /data/me/lung_atlas.h5ad --organ lung --species human --dry-run
  stew add /data/me/GSE123456 --source GSE123456
  stew add /data/me/project/data --dry-run      see how stew splits a folder
  stew add /data/me/pilot_xenium --to projects/pilot_2026
  stew add /data/me/analysis --whole             keep a folder together
  stew add /data/me/scratch.h5ad --in-place"""

INIT_EPILOG = """\
Creates <root>/.data-steward/ with catalog.db and config.toml (layout, whether datasets are
moved in, folder group and mode; edit it any time). stew then finds the catalogue when run
from inside <root>, or from anywhere once you either:

  stew init <root> --set-default            for yourself (~/.config/stew/config.toml)
  echo 'root = "<root>"' > /etc/stew/config.toml     for everyone on the machine (admin)
  export STEW_ROOT=<root>                   for a shell or a script

Examples:
  stew init /data/shared --group lab --set-default
  stew init ~/my-datasets --in-place        catalogue data where it is, never move it"""

DOWNLOAD_EPILOG = """\
SOURCE can be:
  GEO series or sample      GSE274673, GSM8454446 (add --raw for the reads on SRA)
  SRA / ENA reads           SRP…, SRX…, SRR…, PRJNA…, ERP…, PRJEB…
  ArrayExpress / BioStudies E-MTAB-5061, S-BIAD1759 (BioImage Archive)
  CELLxGENE collection      cellxgene:<id> or its cellxgene.cziscience.com URL
  Human Cell Atlas project  hca:<id> or its data.humancellatlas.org URL (--raw for FASTQs)
  figshare / Zenodo         figshare:123, zenodo:123, their DOIs or URLs
  anything else             an http(s)/ftp URL; a folder listing URL downloads the folder

Before downloading anything, stew checks that the data is not already here:
  1. a catalogued dataset already has this source (or the GEO series / Zenodo
     record that contains it);
  2. a catalogued dataset mentions the accession in its study, name or path
     (probably the same data, registered before sources were recorded);
  3. someone is downloading it right now;
  4. the destination folder already exists without being a stew download.
Check 1 and 3 always stop the download. Checks 2 and 4 can be overridden with --force.

To add files that a catalogued record has gained (or that were left out with
--files), run the same source again with --missing: stew fetches only the files
not yet in the dataset's folder, into that folder, and updates its size.

Large files (64 MiB and up) are fetched over several connections at once when the
server allows it (--connections, default 8), and up to 3 files at a time (--jobs),
at most 16 connections altogether; set the defaults in config.toml under
[download] (connections, jobs, max_connections). --connections 1 uses one
connection per file. A server that answers "busy" (429/503) is waited for, as long
as it asks, and given fewer connections.

Interrupted downloads resume when you re-run the same command, each part where it
stopped (`<file>.part` and `<file>.part.json`). Published checksums
(md5 or sha256) are verified. Finished downloads are registered with their source,
and the files are read for cell counts, organism and tissue. With --background the
download runs in a tmux session; `stew downloads` shows progress and stew tells
you when it has finished.

Examples:
  stew download GSE274673 --dry-run         see what would be fetched, and where
  stew download GSE274673 --files '*.h5'    only some of the files
  stew download GSE123456 --raw --background   the FASTQs from SRA, in the background
  stew download cellxgene:6f6d381a-7701-4781-935c-db10d30de293 --files '*core*'
  stew download zenodo:5979761 --organ lung --tag tracking
  stew download zenodo:5979761 --missing --files '*.zip'   add files left out before
  stew download https://example.org/atlas.h5ad --name my-atlas
  stew download S-BIAD1759 --jobs 2 --connections 6 --background"""


class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        sys.exit(f"error: {message}\nRun `{self.prog} --help` for help.")


class TopParser(Parser):
    def format_help(self):
        return TOP_HELP


def _meta_args(p, what="it", source=True):
    g = p.add_argument_group("metadata")
    g.add_argument("--organ", help="organ or tissue, e.g. lung")
    g.add_argument("--species", help="species, e.g. human")
    g.add_argument("-m", "--modality", dest="modalities", action="append", default=[],
                   help="modality, e.g. rna, atac, protein, spatial (repeatable)")
    g.add_argument("--platform", help="platform, e.g. xenium, visium, 10x-chromium")
    g.add_argument("--study", help="study or accession, e.g. GSE274673")
    if source:
        g.add_argument("--source", help=f"where {what} came from: URL or accession (see `stew help download`)")
    g.add_argument("--tag", dest="tags", action="append", default=[], help="free-form tag (repeatable)")
    g.add_argument("--desc", dest="description", help="one-line description")
    g.add_argument("--steward", help="person responsible for it (default: owner of the files)")
    g.add_argument("--derived-from", metavar="DATASET", help=f"the catalogued dataset {what} was made from")


def _home_args(p):
    g = p.add_argument_group("where the data goes")
    g.add_argument("--in-place", action="store_true",
                   help="leave the data where it is instead of moving it into the data folder")
    g.add_argument("--to", dest="home", metavar="DIR",
                   help="folder inside the data folder to move it into (default: chosen from its source)")


def _list_args(p):
    """Filters and output options shared by `ls` and `search`."""
    g = p.add_argument_group("filters")
    g.add_argument("--organ", action="append", default=[], help="organ (repeat to match any of several)")
    g.add_argument("--species", action="append", default=[], help="species (repeat to match any of several)")
    g.add_argument("-m", "--modality", action="append", default=[],
                   help="modality (repeat to require all: -m rna -m spatial)")
    g.add_argument("--platform", action="append", default=[], help="platform (repeat to match any of several)")
    g.add_argument("--format", action="append", default=[],
                   help="file format, e.g. anndata, xenium-bundle (repeat to match any)")
    g.add_argument("--tag", action="append", default=[], help="tag (repeat to require all)")
    g.add_argument("--owner", help="owned or stewarded by this user")
    g.add_argument("--study", help="study or accession")
    g.add_argument("--source", help="came from this URL or accession")
    g.add_argument("--name", help="name contains this text")
    g.add_argument("--under", type=Path, help="path lies under this folder")
    g.add_argument("--text", action="append", default=[], metavar="WORD",
                   help="name, description, study, source or path contains this word (repeat to require all)")
    g.add_argument("--kind", choices=("dataset", "reference", "all"),
                   help="datasets (default), references (genomes, indices, annotations), or all")
    g.add_argument("-q", "--saved", metavar="NAME", help="start from a saved query (see `stew saved`)")
    o = p.add_argument_group("output")
    o.add_argument("-l", "--long", action="store_true", help="also show cells, format, study, source, tags, path")
    x = o.add_mutually_exclusive_group()
    x.add_argument("--json", action="store_true", help="print JSON")
    x.add_argument("--csv", action="store_true", help="print CSV")
    x.add_argument("--paths", action="store_true", help="print only paths, one per line")
    o.add_argument("--save", metavar="NAME", help="save these filters as a query you can rerun with -q NAME")
    o.add_argument("--shared", action="store_true", help="with --save: let everyone use the saved query")


def _pick_args(p, what):
    p.add_argument("datasets", nargs="*", metavar="DATASET", help=f"datasets to {what} (name, id or path)")
    p.add_argument("--all", action="store_true", help=f"{what} every dataset")


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--catalog", type=Path, default=argparse.SUPPRESS, help="catalogue file to use")
    common.add_argument("-y", "--yes", action="store_true", default=argparse.SUPPRESS,
                        help="don't ask questions; accept the suggested answers")

    top = TopParser(prog="stew", add_help=False)
    top.add_argument("-h", "--help", action="help")
    top.add_argument("-V", "--version", action="version", version=f"stew {__version__}")
    top.add_argument("--catalog", type=Path, default=None)
    top.add_argument("-y", "--yes", action="store_true")
    sub = top.add_subparsers(dest="cmd", metavar="<command>", parser_class=Parser)

    def cmd(name, help, epilog=None, aliases=()):
        return sub.add_parser(name, help=help, description=help, epilog=epilog, aliases=list(aliases),
                              parents=[common], formatter_class=argparse.RawDescriptionHelpFormatter)

    p = cmd("init", "Set up a data folder: an empty catalogue and its config", INIT_EPILOG)
    p.add_argument("root", nargs="?", type=Path, help="the data folder (created if missing)")
    p.add_argument("--group", help="Unix group whose members may edit the catalogue (e.g. lab)")
    p.add_argument("--in-place", action="store_true",
                   help="never move datasets into the folder; `stew add` registers them where they are")
    p.add_argument("--set-default", action="store_true",
                   help="also make this your default catalogue (~/.config/stew/config.toml)")
    p.add_argument("--download-root", type=Path,
                   help="where `stew download` puts files (default: <data root>/downloads)")
    p.add_argument("--steward", action="append", default=[], metavar="USER",
                   help="restrict changes to these users (repeatable; default: anyone who can write the file)")

    p = cmd("add", "Register a file, a dataset folder, or every dataset inside a folder", ADD_EPILOG,
            aliases=["scan"])
    p.add_argument("path", type=Path, help="file or folder")
    p.add_argument("--name", help="name for a single dataset (default: from the file or folder name)")
    p.add_argument("--dry-run", action="store_true", help="show what would happen; move and register nothing")
    g = p.add_argument_group("folders")
    g.add_argument("--whole", action="store_true", help="register the folder as one dataset, whatever is inside")
    g.add_argument("--each", action="store_true",
                   help="register each folder directly inside it as one dataset, without looking inside")
    g.add_argument("--max-depth", type=int, default=8, help="how many folder levels to search (default 8)")
    _home_args(p)
    _meta_args(p)

    p = cmd("ls", "List datasets, optionally filtered", LS_EPILOG, aliases=["list"])
    _list_args(p)

    p = cmd("search", "Find datasets by words in their name, description, study, source or path",
            "Every word must appear somewhere (case doesn't matter); all `ls` filters work too.\n\n"
            "Examples:\n  stew search fibroblast\n  stew search caf pdac --species human\n"
            "  stew search GSE274673 --paths")
    p.add_argument("words", nargs="+", help="words to look for")
    _list_args(p)

    p = cmd("saved", "List your saved queries (and shared ones), or delete one")
    p.add_argument("--rm", metavar="NAME", help="delete this saved query")
    p.add_argument("--shared", action="store_true", help="with --rm: delete the shared query of that name")

    p = cmd("info", "Show everything known about one dataset")
    p.add_argument("dataset", help="dataset name, id or path")
    p.add_argument("--json", action="store_true", help="print JSON")

    p = cmd("set", "Change a dataset's metadata", SET_EPILOG)
    p.add_argument("dataset", help="dataset name, id or path")
    p.add_argument("--name", help="new name")
    for f in ("organ", "species", "platform", "study", "steward"):
        p.add_argument(f"--{f}", help=f"new {f} (empty string clears it)")
    p.add_argument("--desc", dest="description", help="new description (empty string clears it)")
    p.add_argument("--source", help="where it came from: URL, GEO accession or Zenodo record (empty clears)")
    p.add_argument("-m", "--modality", action="append", default=[], help="replace all modalities (repeatable)")
    p.add_argument("--add-modality", action="append", default=[], help="add a modality (repeatable)")
    p.add_argument("--rm-modality", action="append", default=[], help="remove a modality (repeatable)")
    p.add_argument("--tag", action="append", default=[], help="add a tag (repeatable)")
    p.add_argument("--untag", action="append", default=[], help="remove a tag (repeatable)")

    p = cmd("rm", "Remove a dataset from the catalogue (its files are never touched)")
    p.add_argument("dataset", help="dataset name, id or path")

    p = cmd("get", "Put a symlink to a dataset into one of your folders", GET_EPILOG)
    p.add_argument("dataset", help="dataset name, id or path")
    p.add_argument("--to", type=Path, default=Path("."), help="folder to put the link in (default: here)")
    p.add_argument("--as", dest="link_name", help="name of the link (default: the dataset's file or folder name)")

    p = cmd("links", "List links handed out by `get`")
    p.add_argument("dataset", nargs="?", help="only links to this dataset")
    p.add_argument("--all", action="store_true", help="everyone's links, not just yours")
    p.add_argument("--released", action="store_true", help="include links that were released")
    p.add_argument("--prune", action="store_true",
                   help="release links whose symlink was deleted by hand (yours; with --all, everyone's)")

    p = cmd("release", "Remove a link made by `get`")
    p.add_argument("link", type=Path, help="the link to remove")

    p = cmd("download", "Download a dataset unless it is already here", DOWNLOAD_EPILOG)
    p.add_argument("source", help="accession, record id or URL (see below)")
    p.add_argument("--raw", action="store_true", help="GEO/HCA: the raw reads instead of the processed files")
    p.add_argument("--background", action="store_true",
                   help="after the checks, download in a tmux session and return at once")
    p.add_argument("--to", type=Path,
                   help="folder to download into (default: <data root>/downloads/<kind>/<id>)")
    p.add_argument("--name", help="name for the new dataset (default: the accession or file name)")
    p.add_argument("--files", action="append", default=[], metavar="PATTERN",
                   help="only files matching this pattern, e.g. '*.h5ad' (repeatable)")
    p.add_argument("--dry-run", action="store_true", help="run the checks and show the plan; download nothing")
    p.add_argument("--force", action="store_true",
                   help="download even if a dataset merely mentions the accession, or the folder exists")
    p.add_argument("--missing", action="store_true",
                   help="the source is catalogued already: fetch only its files not yet in the dataset's "
                        "folder, into that folder")
    p.add_argument("--connections", type=positive, metavar="N",
                   help="connections per large file (default 8; 1 = one stream, as before 0.7)")
    p.add_argument("--jobs", type=positive, metavar="N",
                   help="files downloaded at once (default 3); at most 16 connections altogether")
    _meta_args(p, source=False)

    p = cmd("downloads", "Downloads in progress, finished and failed")
    p.add_argument("-n", "--limit", type=int, default=20, help="how many to show (default 20)")

    p = cmd("inspect", "Read cell/feature counts, organism, tissue and assay from inside the files",
            "stew reads only headers, indices and small annotation columns, never whole matrices.\n"
            "Fields that are already set are kept unless you pass --overwrite.\n\n"
            "Examples:\n  stew inspect hlca-core\n  stew inspect --all\n  stew inspect 7 --overwrite")
    _pick_args(p, "inspect")
    p.add_argument("--overwrite", action="store_true", help="replace organ/species/platform already set")

    p = cmd("relate", "Record that one dataset comes from, or is part of, another",
            "Examples:\n  stew relate hlca-core-filtered hlca-core                 derived from (default)\n"
            "  stew relate sample-3 gse274673 --kind part-of\n"
            "  stew relate copy-on-scratch hlca-core --kind same-as")
    p.add_argument("child", nargs="?", help="the dataset that came from the other")
    p.add_argument("parent", nargs="?", help="the dataset it came from or belongs to")
    p.add_argument("--kind", default="derived-from", choices=("derived-from", "part-of", "same-as"))
    p.add_argument("--suggest", action="store_true",
                   help="list relations the catalogue implies (accessions one dataset mentions that another "
                        "comes from, GEO samples of a catalogued series, 10x pipeline inputs)")
    p.add_argument("--apply", action="store_true", help="with --suggest: record the suggested relations")
    p = cmd("unrelate", "Remove a relation recorded with `relate`")
    p.add_argument("child")
    p.add_argument("parent")

    p = cmd("checksum", "Record a checksum of every file in datasets",
            "Full checksums read every byte; --quick reads the first and last MiB of each file (fast,\n"
            "good for finding copies). Run big ones in tmux. `stew verify` re-checks later.\n\n"
            "Examples:\n  stew checksum hlca-core\n  stew checksum --all --quick\n  stew checksum --all --jobs 8")
    _pick_args(p, "checksum")
    p.add_argument("--quick", action="store_true", help="hash only the first and last MiB of each file")
    p.add_argument("--jobs", type=int, default=4, help="files to hash at once (default 4)")
    p = cmd("verify", "Check datasets against their recorded checksums")
    _pick_args(p, "verify")
    p.add_argument("--jobs", type=int, default=4, help="files to hash at once (default 4)")
    cmd("duplicates", "Datasets and files that are stored more than once (needs `stew checksum`)")

    p = cmd("serve", "A read-only web page of the catalogue",
            "The page lists datasets with search, filters, details and CSV export. It only reads the\n"
            "catalogue, never the data. By default only this machine can reach it; from your own\n"
            "computer, tunnel through ssh:\n\n"
            "  ssh -N -L 8765:localhost:8765 <you>@<server>     then open http://localhost:8765/")
    p.add_argument("--host", default="127.0.0.1", help="address to listen on (default 127.0.0.1: this machine)")
    p.add_argument("--port", type=int, default=8765, help="port (default 8765)")

    p = cmd("mv", "Move a dataset (or a folder of datasets) and update the catalogue and links",
            "Moves on disk with a rename (same filesystem only, so it is instant), updates the paths of\n"
            "every catalogued dataset at or below the old path, and repoints the links handed out by\n"
            "`stew get`. --leave-link leaves a symlink at the old path so scripts that use it keep working.\n\n"
            "A destination ending in '/' is a folder to move into (created with --parents); without it,\n"
            "a destination that doesn't exist is the new name.\n\n"
            "Examples:\n  stew mv hlca-core /data/shared/public/hlca --dry-run\n"
            "  stew mv /data/shared/old_project /data/shared/projects/ --leave-link")
    p.add_argument("dataset", help="dataset name, id or path, or a folder holding catalogued datasets")
    p.add_argument("dest", help="new path; an existing folder, or one ending in '/', to move it into")
    p.add_argument("--dry-run", action="store_true", help="show what would happen")
    p.add_argument("--parents", action="store_true", help="create missing parent folders of the destination")
    p.add_argument("--leave-link", action="store_true", help="leave a symlink at the old path")

    p = cmd("refresh", "Re-check that datasets still exist and update their size and owner")
    p.add_argument("dataset", nargs="?", help="only this dataset (default: all)")

    cmd("summary", "Counts by modality, organ, species and platform")
    cmd("vocab", "Known modalities, organs, species and platforms, with accepted aliases")

    p = cmd("stewards", "Who may change the catalogue; stewards can add or remove stewards",
            "Everyone can find datasets, get and release links, and save queries. Only stewards can\n"
            "register, edit, remove, move, inspect, checksum or download datasets.\n\n"
            "Examples:\n  stew stewards\n  stew stewards --add gw\n  stew stewards --remove gw")
    p.add_argument("--add", action="append", default=[], metavar="USER")
    p.add_argument("--remove", action="append", default=[], metavar="USER")

    p = cmd("check", "Health check: datasets present, sizes and owners up to date, checksums still match",
            "Runs `refresh` and `verify` over the whole catalogue and prints only what needs attention,\n"
            "ending with a one-line verdict. Exit status 1 if anything is missing or changed. Meant to\n"
            "run nightly; --report also writes the result to a file.")
    p.add_argument("--report", type=Path, help="also write the result to this file")
    p.add_argument("--jobs", type=int, default=8, help="files to hash at once (default 8)")

    p = cmd("report", "Storage report: where the space goes, and what nobody uses",
            "Space by area of the data folder (projects/, public/geo/, references/, ...), by owner and by\n"
            "modality; the largest datasets; datasets nobody has linked with `stew get` whose files have not\n"
            "changed for --stale-days; and free space on the disk. Sizes are from the last refresh/check.")
    p.add_argument("-n", "--top", type=int, default=10, help="how many of the largest/unused to list (default 10)")
    p.add_argument("--stale-days", type=int, default=180, help="unused means untouched this long (default 180)")
    p.add_argument("--depth", type=int, default=2, help="folder levels for the area breakdown (default 2)")

    p = cmd("backup", "Copy the catalogue to .data-steward/backups/catalog-<date>.db")
    p.add_argument("--keep", type=int, default=30, help="how many backups to keep (default 30)")
    p.add_argument("--to", type=Path, help="backup folder (default: backups/ next to the catalogue)")

    p = cmd("config", "Which catalogue stew uses, how it found it, and its settings")
    p.add_argument("--get", choices=("root", "catalog", "layout", "move_on_add", "group"),
                   help="print just one value (for scripts)")

    p = cmd("log", "Recent changes to the catalogue")
    p.add_argument("-n", "--limit", type=int, default=20, help="how many changes to show (default 20)")

    p = sub.add_parser("help", help="Show help for a command")
    p.add_argument("topic", nargs="?")
    return top


# --- small helpers -------------------------------------------------------------------------

def opt(v) -> str:
    return "-" if v in (None, "", []) else ", ".join(v) if isinstance(v, list) else str(v)


def table(header: list[str], rows: list[list]) -> None:
    rows = [[opt(c) if not isinstance(c, str) else (c or "-") for c in r] for r in rows]
    w = [max([len(h)] + [len(r[i]) for r in rows]) for i, h in enumerate(header)]
    print("  ".join(h.ljust(w[i]) for i, h in enumerate(header)).rstrip())
    for r in rows:
        print("  ".join(c.ljust(w[i]) for i, c in enumerate(r)).rstrip())


def modalities(raw: list[str]) -> list[str]:
    out = set()
    for r in raw:
        for part in r.split(","):
            if part.strip():
                out.add(vocab.modality(part))
    return sorted(out)


def tags(raw: list[str]) -> list[str]:
    out = []
    for t in raw:
        t = t.strip().lower()
        if not t or any(c.isspace() for c in t) or "," in t:
            raise UserError(f"tags must be single words without commas: '{t}'")
        out.append(t)
    return out


def check_name(name: str) -> None:
    if not name or any(c.isspace() for c in name) or "/" in name:
        raise UserError(f"dataset names must be non-empty, without spaces or '/': '{name}'")
    if name.isdigit():
        raise UserError(f"dataset names can't be plain numbers (those are ids): '{name}'")


class Ui:
    def __init__(self, yes: bool):
        self.yes = yes
        self.tty = sys.stdin.isatty() and sys.stderr.isatty()

    @property
    def interactive(self) -> bool:
        return self.tty and not self.yes

    def ask(self, prompt: str, default: str | None = None, check=None) -> str | None:
        while True:
            suffix = f" [{default}]" if default else ""
            try:
                v = input(f"  {prompt}{suffix}: ").strip()
            except EOFError:
                v = ""
            v = v or default or ""
            if check and v:
                try:
                    check(v)
                except (UserError, ValueError) as e:
                    print(f"  {e}")
                    continue
            return v or None

    def confirm(self, prompt: str, default: bool) -> bool:
        if self.yes:
            return True
        if not self.tty:
            return False
        hint = "Y/n" if default else "y/N"
        try:
            v = input(f"{prompt} [{hint}] ").strip().lower()
        except EOFError:
            return default
        return default if not v else v.startswith("y")

    def choose(self, prompt: str, options: list[str]) -> int:
        print(prompt)
        for i, o in enumerate(options, 1):
            print(f"  {i}) {o}")
        while True:
            try:
                v = input(f"Choose 1-{len(options)} [1]: ").strip() or "1"
            except EOFError:
                return len(options) - 1
            if v.isdigit() and 1 <= int(v) <= len(options):
                return int(v) - 1


class Meta:
    """Normalised metadata given on the command line for scan, add and download."""

    def __init__(self, a):
        self.organ = vocab.organ(a.organ) if a.organ else None
        self.species = vocab.species(a.species) if a.species else None
        self.platform = vocab.platform(a.platform) if a.platform else None
        self.modalities = modalities(a.modalities)
        self.study = a.study
        self.tags = tags(a.tags)
        self.description = a.description
        self.steward = a.steward
        self.source = sources.parse(a.source) if getattr(a, "source", None) else None
        self.derived_from = a.derived_from

    def apply(self, c: Candidate) -> None:
        c.organ = self.organ or c.organ
        c.species = self.species or c.species
        if self.platform:
            c.platform = self.platform
            c.add_platform_modalities()
        c.modalities |= set(self.modalities)


def ask_metadata(ui: Ui, c: Candidate, everything: bool) -> None:
    """Asks about the fields a candidate lacks (or all of them)."""
    if everything:
        c.name = ui.ask("name", c.name, check_name) or c.name
    if everything or not c.organ:
        v = ui.ask("organ", c.organ)
        c.organ = vocab.organ(v) if v else None
    if everything or not c.species:
        v = ui.ask("species", c.species)
        c.species = vocab.species(v) if v else None
    if everything or not c.modalities:
        v = ui.ask("modalities (comma-separated)", ",".join(sorted(c.modalities)) or None,
                   lambda s: modalities([s]))
        c.modalities = set(modalities([v])) if v else set()


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def new_dataset(cat: Catalog, c: Candidate, meta: Meta, src: sources.Source | None = None) -> dict:
    owner = util.owner_of(c.path)
    size, n = util.disk_usage(c.path)
    src = src or meta.source
    mk = c.meta  # from a .stew-dataset file; command-line values win
    if not src and mk.get("source"):
        try:
            src = sources.parse(str(mk["source"]))
        except sources.SourceError:
            pass
    inspected = c.n_cells is not None or c.n_features is not None or bool(c.details)
    tag_list = list(meta.tags) + [t for t in (mk.get("tags") or []) if isinstance(t, str) and t not in meta.tags]
    return {
        "kind": c.kind,
        "n_cells": c.n_cells, "n_features": c.n_features,
        "details": json.dumps(c.details) if c.details else None, "inspected_at": now() if inspected else None,
        "name": cat.unique_name(c.name, detect.parent_name(c.path)),
        "path": str(c.path), "platform": c.platform, "format": c.format,
        "organ": c.organ, "species": c.species, "study": meta.study or mk.get("study"),
        "description": meta.description or mk.get("description"),
        "owner": owner, "steward": meta.steward or mk.get("steward") or owner, "size_bytes": size, "n_files": n,
        "source": src.display if src else None, "source_key": src.key if src else None,
        "modalities": sorted(c.modalities), "tags": tag_list,
    }


def register(cat: Catalog, c: Candidate, meta: Meta, src: sources.Source | None = None,
             parent: dict | None = None, kind: str = "derived-from") -> tuple[int, dict]:
    """Inserts a candidate, and relates it to `meta.derived_from` or `parent`."""
    d = new_dataset(cat, c, meta, src)
    did = cat.insert(d)
    up = parent or (cat.find(meta.derived_from) if meta.derived_from else None)
    if up:
        cat.relate(cat.by_id(did), up, kind if parent else "derived-from")
    return did, d


_SETTINGS: dict[str, config.Settings] = {}


def settings(cat: Catalog) -> config.Settings:
    key = str(cat.path)
    if key not in _SETTINGS:
        try:
            _SETTINGS[key] = config.load(cat.path)
        except config.ConfigError as e:
            raise UserError(str(e))
    return _SETTINGS[key]


def data_home(cat: Catalog) -> Path:
    return settings(cat).root


def layout_folder(cat: Catalog, key: str, **fields) -> Path:
    try:
        return data_home(cat) / settings(cat).layout[key].format(**fields)
    except (KeyError, IndexError) as e:
        raise UserError(f"[layout] {key} = {settings(cat).layout[key]!r} uses an unknown field {e}")


def plan_home(cat: Catalog, c: Candidate, src: sources.Source | None, to: str | None) -> Path:
    """Where in the data folder a dataset should live."""
    home = data_home(cat)
    name = c.path.name
    if to:
        base = Path(to).expanduser()
        base = base if base.is_absolute() else home / base
        return base / name
    if c.kind == "reference":
        ref_kind = "10x" if c.format == "10x-reference" else "annotations"
        return layout_folder(cat, "reference", ref_kind=ref_kind, name=name) / name
    if src and src.kind != "url" and src.accession:
        folder = layout_folder(cat, "public", kind=src.kind, accession=src.accession, name=name)
        return folder if util.slug(Path(name).stem) == util.slug(src.accession) else folder / name
    return layout_folder(cat, "project", name=name) / name


def move_problem(src: Path, dest: Path) -> str | None:
    """Why `src` can't be renamed to `dest`, or None."""
    if dest.exists() or dest.is_symlink():
        return f"{dest} already exists"
    if dest == src or dest.is_relative_to(src):
        return "can't move a folder into itself"
    probe = dest.parent
    while not probe.exists():
        probe = probe.parent
    if os.stat(probe).st_dev != os.lstat(src).st_dev:
        return (f"{src} is on a different disk from {dest.parent}, and stew only moves by renaming. Copy it "
                "into the data folder yourself (e.g. with rsync), check the copy, then `stew add` the copy.")
    need = {src.parent, probe}
    if src.is_dir() and not src.is_symlink():
        need.add(src)  # its ".." entry is rewritten
    blocked = sorted(str(p) for p in need if not os.access(p, os.W_OK))
    if blocked:
        return ("you don't have permission to move it there (no write access to "
                + ", ".join(f"{p}, owner {util.owner_of(Path(p))}" for p in blocked)
                + "). Ask a steward to add it, or use --in-place to register it where it is.")
    return None


def make_shared_dirs(folder: Path, mode: int = 0o3775, group: str | None = None) -> None:
    """Creates missing folders the way the shared data folder's are set up
    (by default group-writable, inheriting the group (setgid), and sticky, so
    everyone can add datasets but only a dataset's owner can delete or rename it)."""
    import grp
    gid = -1
    if group:
        try:
            gid = grp.getgrnam(group).gr_gid
        except KeyError:
            raise UserError(f"config: no such group {group!r}")
    missing = []
    while not folder.exists():
        missing.append(folder)
        folder = folder.parent
    for d in reversed(missing):
        d.mkdir()
        if gid != -1:
            try:
                os.chown(d, -1, gid)
            except OSError:
                pass
        for m in (mode, mode & 0o777):  # some sandboxes refuse the setgid/sticky bits
            try:
                os.chmod(d, m)
                break
            except OSError:
                continue


def shared_dirs(cat: Catalog, folder: Path) -> None:
    s = settings(cat)
    make_shared_dirs(folder, s.folder_mode, s.group)


def move_home(src: Path, dest: Path, cat: Catalog | None = None) -> None:
    """Renames src to dest and leaves a symlink at src."""
    if cat is not None:
        shared_dirs(cat, dest.parent)
    else:
        make_shared_dirs(dest.parent)
    os.rename(src, dest)
    try:
        src.symlink_to(dest)
    except OSError:
        os.rename(dest, src)
        raise


def homes_for(cat: Catalog, cs: list[Candidate], meta: "Meta", a) -> list[tuple[Candidate, Path | None]]:
    """(candidate, destination or None to stay) for each; raises if one can't move."""
    home = data_home(cat)
    out, problems = [], []
    for c in cs:
        if a.in_place or not settings(cat).move_on_add or c.path.is_relative_to(home):
            out.append((c, None))
            continue
        src = meta.source or (sources.parse(c.meta["source"]) if c.meta.get("source") else None)
        dest = plan_home(cat, c, src, a.home)
        why = move_problem(c.path, dest)
        if why:
            problems.append(f"{c.path}: {why}")
        out.append((c, dest))
    if problems:
        raise UserError("can't move into the data folder:\n  " + "\n  ".join(problems))
    return out


def link_provenance(cat: Catalog, d: dict) -> list[str]:
    """Relates a pipeline output to the catalogued FASTQs and reference its
    `_cmdline` names. Returns what was linked, for printing."""
    p = Path(d["path"])
    if not p.is_dir():
        return []
    prov = provenance.read(p)
    if prov.empty():
        return []
    done = []
    for paths, kind in ((prov.inputs, "derived-from"), (prov.references, "uses")):
        for q in paths:
            up = cat.containing(q)
            if up and up["id"] != d["id"] and not Path(up["path"]).is_relative_to(p):
                try:
                    cat.relate(d, up, kind)
                    done.append(f"{kind} {up['name']}")
                except CatalogError:
                    continue
                # Raw reads rarely say what they are; the output made from them does.
                fill = [(col, d[col]) for col in ("species", "organ")
                        if kind == "derived-from" and d[col] and not up[col]]
                if fill and up["kind"] == "dataset":
                    cat.update(up, Changes(fields=fill))
    details = dict(d.get("details") or {})
    if details.get("pipeline") != prov.as_details():
        details["pipeline"] = prov.as_details()
        cat.update(d, Changes(fields=[("details", json.dumps(details))]))
    return done


def link_all(cat: Catalog, ids: list[int]) -> None:
    lines = []
    for did in ids:
        d = cat.by_id(did)
        links = link_provenance(cat, d) if d else []
        if links:
            lines.append(f"  #{did} {d['name']}: " + ", ".join(links))
    if lines:
        print("Linked pipeline outputs to the inputs named in their _cmdline:")
        print("\n".join(lines))


def incomplete(d: dict) -> bool:
    return not (d["organ"] and d["species"] and d["modalities"])


def nag_missing(added: list[tuple[int, dict]]) -> None:
    missing = [(i, d) for i, d in added if incomplete(d) and d.get("kind", "dataset") == "dataset"]
    if missing:
        names = ", ".join(f"#{i} {d['name']}" for i, d in missing)
        print(f"{len(missing)} of them lack organ, species or modality: {names}\n"
              f"Fill these in with e.g. `stew set {missing[0][1]['name']} --organ lung --species human -m rna`.")


def warn_same_source(cat: Catalog, src: sources.Source | None, exclude: int | None = None) -> None:
    if not src:
        return
    others = [d for d in cat.by_source([src.key] + src.related) if d["id"] != exclude]
    if others:
        names = ", ".join(f"#{d['id']} {d['name']}" for d in others)
        print(f"note: {names} already come(s) from {src.display}")


# --- commands ------------------------------------------------------------------------------

def cmd_init(a, ui):
    if a.root and a.catalog:
        raise UserError("give either the data folder or --catalog, not both")
    target = a.root or a.catalog
    if target is None:
        raise UserError("which folder should hold the data? e.g. `stew init /data/shared`")
    target = target.expanduser()
    if a.root or (target.is_dir() and not target.name.endswith(".db")):
        target.mkdir(parents=True, exist_ok=True)
        catalog_path = config.catalog_in(target.resolve())
    else:
        catalog_path = target
    a.catalog = catalog_path
    cat = Catalog.create(catalog_path, a.group)
    cfg = config.write_config(catalog_path, a.group, move_on_add=not a.in_place)
    if a.download_root:
        cat.set_meta("download_root", str(a.download_root.expanduser().resolve()))
    if a.steward:
        cat.set_meta("stewards", ",".join(a.steward))
    print(f"Created an empty catalogue at {catalog_path}\nSettings: {cfg} (layout, moving, folder group/mode)")
    if a.set_default:
        print(f"Made it your default catalogue: {config.write_user_config(catalog_path.resolve().parent.parent)}")
    else:
        print(f"stew finds it from inside {catalog_path.resolve().parent.parent}; to use it from anywhere, "
              "re-run with --set-default, set $STEW_ROOT, or ask an admin for /etc/stew/config.toml.")
    if a.steward:
        print(f"Stewards (may register, edit, move and download): {', '.join(a.steward)}")
    if a.group:
        print(f"Members of group '{a.group}' can use it: find datasets, get links, save queries.")
    print("Next: stew add <folder>")


def print_candidates(root: Path, cs: list[Candidate]) -> None:
    rows = []
    for i, c in enumerate(cs, 1):
        try:
            rel = str(c.path.relative_to(root)) or "."
        except ValueError:
            rel = str(c.path)
        rows.append([str(i), rel, c.platform or "?", c.format, ",".join(sorted(c.modalities)) or "?",
                     c.organ or "?", c.species or "?", f"{c.n_cells:,}" if c.n_cells else "-"])
    table(["#", "PATH", "PLATFORM", "FORMAT", "MODALITIES", "ORGAN", "SPECIES", "CELLS"], rows)


def inside_catalogued(cat: Catalog) -> callable:
    """A test for whether a path lies inside a dataset that is already catalogued."""
    folders = {Path(d["path"]) for d in cat.list()}

    def inside(p: Path) -> Path | None:
        return next((q for q in p.parents if q in folders), None)
    return inside


def cmd_scan(a, ui, cat):
    """Registers every dataset found under a folder (reached through `add`)."""
    root = a.path.resolve()
    meta = Meta(a)
    if meta.derived_from:
        cat.find(meta.derived_from)
    if a.each:
        report = detect.ScanReport()
        for e in sorted(root.iterdir()):
            if e.is_dir() and not e.is_symlink() and not e.name.startswith("."):
                report.found.append(detect.detect_any(e))
    else:
        report = detect.scan(root, a.max_depth)
    for p in report.unreadable:
        print(f"warning: no permission to look inside {p}", file=sys.stderr)
    inside = inside_catalogued(cat)
    known = within = 0
    fresh = []
    for c in report.found:
        if cat.by_path(c.path):
            known += 1
        elif inside(c.path):
            within += 1
        else:
            c.read_contents()
            meta.apply(c)
            fresh.append(c)
    notes = [f"{known} already catalogued"] * bool(known) + [f"{within} inside catalogued datasets"] * bool(within)
    already = f" ({'; '.join(notes)})" if notes else ""
    if not fresh:
        print(f"No new datasets under {root}{already}.")
        return
    print(f"Found {len(fresh)} new dataset(s) under {root}{already}:\n")
    print_candidates(root, fresh)
    try:
        plan = dict((id(c), dest) for c, dest in homes_for(cat, fresh, meta, a))
    except UserError as e:
        if a.dry_run:
            print(f"\n{e}")
            plan = {}
        else:
            raise
    moving = [(c, plan[id(c)]) for c in fresh if plan.get(id(c))]
    if moving:
        print(f"\n{len(moving)} of them will be moved into {data_home(cat)}, leaving a link at the old path:")
        for c, dest in moving:
            print(f"  {c.path}  ->  {dest}")
    for c in fresh:
        if c.note:
            print(f"note: {c.path.name}: {c.note}", file=sys.stderr)
    print()
    if a.dry_run:
        print("Dry run: nothing registered.")
        return
    if ui.yes:
        chosen = fresh
    elif not ui.tty:
        print("Not running in a terminal, so not registering. Re-run with --yes to register all of these.")
        return
    else:
        choice = ui.choose("Register them?", ["yes, all as shown", "review one by one", "no"])
        if choice == 2:
            print("Nothing registered.")
            return
        chosen = fresh
        if choice == 1:
            chosen = []
            for c in fresh:
                print(f"\n{c.path}  [{c.platform or '?'} {c.format}]")
                if ui.confirm("Register this one?", True):
                    ask_metadata(ui, c, True)
                    chosen.append(c)
    added = []
    for c in chosen:
        dest = plan.get(id(c))
        if dest:
            move_home(c.path, dest, cat)
            c.path = dest
        added.append(register(cat, c, meta))
    if moving:
        print(f"Moved {sum(1 for c in chosen if plan.get(id(c)))} dataset(s) into {data_home(cat)}; "
              "the old paths are links to them.")
    refs = sum(1 for _, d in added if d["kind"] == "reference")
    print(f"Registered {len(added)} dataset(s)" + (f", {refs} of them references." if refs else "."))
    link_all(cat, [i for i, _ in added])
    if meta.source and len(added) > 1:
        print(f"All of them were given the source {meta.source.display}.")
    nag_missing(added)


def cmd_add(a, ui, cat):
    try:
        path = a.path.resolve(strict=True)
    except OSError:
        raise UserError(f"cannot open {a.path}")
    d = cat.by_path(path)
    if d:
        raise UserError(f"{path} is already catalogued as #{d['id']} {d['name']}")
    meta = Meta(a)
    if meta.derived_from:
        cat.find(meta.derived_from)
    if a.whole and a.each:
        raise UserError("--whole and --each contradict each other; pick one")
    c = detect.detect_any(path)
    # A folder that isn't itself a dataset: register what is inside instead.
    # With nothing recognisable inside, the folder itself is the dataset.
    if path.is_dir() and not a.whole and (a.each or c.format == "directory"):
        inside = [] if a.each else detect.scan(path, a.max_depth).found
        if a.each or inside:
            if a.name:
                raise UserError(f"{path} holds several datasets, so --name doesn't apply; name them afterwards "
                                "with `stew set <dataset> --name ...`, or add the folder with --whole")
            print(f"{path} isn't a dataset itself, so stew looks for datasets inside it "
                  "(--whole to add the folder as one).\n" if not a.each else "", end="")
            return cmd_scan(a, ui, cat)
    c.read_contents()
    if c.note:
        print(f"note: {c.note}", file=sys.stderr)
    meta.apply(c)
    if a.name:
        check_name(a.name)
        if cat.name_taken(a.name):
            raise UserError(f"the name '{a.name}' is already used")
        c.name = a.name
    (_, dest), = homes_for(cat, [c], meta, a)
    if a.dry_run:
        print(f"Would register '{c.name}': {c.format}, platform {opt(c.platform)}, modalities "
              f"{opt(sorted(c.modalities))}, organ {opt(c.organ)}, species {opt(c.species)}")
        print(f"and move {path}\n  -> {dest}\nleaving a link at the old path." if dest
              else f"in place, at {path}.")
        print("Dry run: nothing moved or registered.")
        return
    if ui.interactive:
        ask_metadata(ui, c, False)
    warn_same_source(cat, meta.source)
    if dest:
        move_home(path, dest, cat)
        c.path = dest
        print(f"Moved {path} -> {dest} (the old path is now a link to it)")
    did, d = register(cat, c, meta)
    print(f"Registered #{did} {d['name']}")
    link_all(cat, [did])
    nag_missing([(did, d)])
    print_info(cat, cat.by_id(did))


FILTER_KEYS = ("organ", "species", "modality", "platform", "format", "tag", "owner", "study", "source", "name",
               "under", "text", "kind")


def filter_values(a, cat: Catalog) -> dict:
    """The filters given on the command line, merged with a saved query (-q)."""
    vals = {k: getattr(a, k) for k in FILTER_KEYS}
    vals["under"] = str(a.under.resolve()) if a.under else None
    vals["text"] = list(a.text) + list(getattr(a, "words", []) or [])
    if a.saved:
        saved = cat.saved_query(a.saved)
        if saved is None:
            names = ", ".join(q["name"] for q in cat.saved_queries()) or "none yet"
            raise UserError(f"no saved query '{a.saved}' (saved: {names})")
        for k, v in saved.items():
            if isinstance(v, list):
                vals[k] = v + [x for x in vals.get(k) or [] if x not in v]
            elif not vals.get(k):
                vals[k] = v
    return vals


def run_filter(cat: Catalog, v: dict) -> tuple[list[dict], Filter]:
    src_key = src_text = None
    if v["source"]:
        try:
            src_key = sources.parse(v["source"]).key
        except sources.SourceError:
            src_text = v["source"]
    f = Filter(
        organs=[vocab.organ(s) for s in v["organ"]], species=[vocab.species(s) for s in v["species"]],
        modalities=modalities(v["modality"]), platforms=[vocab.platform(s) for s in v["platform"]],
        formats=[s.lower() for s in v["format"]], tags=tags(v["tag"]), owner=v["owner"], study=v["study"],
        name=v["name"], under=v["under"], source_key=src_key, source_text=src_text,
        kinds=[] if v.get("kind") == "all" else [v.get("kind") or "dataset"],
    )
    rows = cat.list(f)
    words = [w.lower() for w in v["text"]]
    if words:
        def hay(d):
            return " ".join(str(x) for x in (d["name"], d["description"], d["study"], d["source"], d["path"],
                                              d["organ"], d["species"], d["platform"], d["format"], d["tags"],
                                              d["modalities"], json.dumps(d.get("details") or {}))).lower()
        rows = [d for d in rows if all(w in hay(d) for w in words)]
    return rows, f


CSV_FIELDS = ("id", "name", "path", "platform", "format", "modalities", "organ", "species", "study", "source",
              "n_cells", "n_features", "size_bytes", "n_files", "owner", "steward", "tags", "description")


def print_rows(a, cat: Catalog, rows: list[dict], v: dict, f: Filter) -> None:
    if a.json:
        print(json.dumps(rows, indent=2))
        return
    if a.csv:
        w = csv.writer(sys.stdout)
        w.writerow(CSV_FIELDS)
        for d in rows:
            w.writerow([";".join(d[k]) if isinstance(d[k], list) else ("" if d[k] is None else d[k])
                        for k in CSV_FIELDS])
        return
    if a.paths:
        for d in rows:
            print(d["path"])
        return
    if not rows:
        if f.is_empty() and not v["text"]:
            print("The catalogue is empty. Add datasets with `stew add <folder>` or `stew download <source>`.")
            return
        print("No datasets match.")
        s = cat.summary()
        for label, key in (("organ", "organ"), ("modality", "modality"), ("species", "species"),
                           ("platform", "platform")):
            if v[key] and s[label]:
                print(f"  {label} values in the catalogue: " + ", ".join(f"{k} ({n})" for k, n in s[label]))
        return
    header = ["ID", "NAME", "PLATFORM", "MODALITIES", "ORGAN", "SPECIES", "SIZE", "OWNER"]
    if v.get("kind") in ("reference", "all"):
        header[2] = "PLATFORM/FORMAT"
    if a.long:
        header += ["CELLS", "FORMAT", "STUDY", "SOURCE", "TAGS", "PATH"]
    out = []
    for d in rows:
        r = [str(d["id"]), d["name"] + (" [ref]" if d["kind"] == "reference" else ""),
             d["platform"] or (d["format"] if d["kind"] == "reference" else None), ",".join(d["modalities"]),
             d["organ"], d["species"],
             human_size(d["size_bytes"]), d["owner"]]
        if a.long:
            r += [f"{d['n_cells']:,}" if d["n_cells"] else "-", d["format"], d["study"], d["source"],
                  ",".join(d["tags"]), d["path"]]
        out.append(r)
    table(header, out)
    total = sum(d["size_bytes"] or 0 for d in rows)
    hidden = ""
    if not v.get("kind"):
        n = len(cat.list(replace(f, kinds=["reference"])))
        hidden = f" ({n} matching reference(s) not shown; add --kind reference or --kind all)" if n else ""
    print(f"{len(rows)} dataset(s), {human_size(total)} in total{hidden}", file=sys.stderr)


def cmd_ls(a, ui, cat):
    v = filter_values(a, cat)
    rows, f = run_filter(cat, v)
    if a.save:
        check_name(a.save)
        saved = {k: x for k, x in v.items() if x}
        if not saved:
            raise UserError("nothing to save: give some filters")
        cat.save_query(a.save, saved, a.shared)
        print(f"Saved as '{a.save}'{' for everyone' if a.shared else ''}; rerun with `stew ls -q {a.save}`.",
              file=sys.stderr)
    print_rows(a, cat, rows, v, f)


def query_text(q: dict) -> str:
    parts = []
    for k, x in q.items():
        flag = "--modality" if k == "modality" else f"--{k}"
        for item in (x if isinstance(x, list) else [x]):
            parts.append(f"{flag} {shlex.quote(str(item))}")
    return " ".join(parts)


def cmd_saved(a, ui, cat):
    if a.rm:
        if not cat.delete_query(a.rm, a.shared):
            raise UserError(f"no {'shared ' if a.shared else ''}saved query '{a.rm}' of yours")
        print(f"Deleted '{a.rm}'.")
        return
    qs = cat.saved_queries()
    if not qs:
        print("No saved queries. Save one with e.g. `stew ls --organ lung -m rna --save lung-rna`.")
        return
    table(["NAME", "BY", "QUERY"], [[q["name"], "everyone" if q["user"] == "*" else q["user"],
                                     query_text(json.loads(q["args"]))] for q in qs])


def print_info(cat: Catalog, d: dict) -> None:
    path = Path(d["path"])
    access = {
        Presence.READABLE: "you can read it",
        Presence.NO_ACCESS: f"you can NOT read it; ask {opt(d['steward'])} for access",
        Presence.MISSING: f"MISSING: the path no longer exists (tell {opt(d['steward'])})",
    }[util.presence(path)]
    links = cat.links(dataset_id=d["id"])
    mine = sum(1 for lk in links if lk["user"] == cat.user)
    size = f"{human_size(d['size_bytes'])} in {d['n_files']} file(s)" if d["size_bytes"] is not None else "-"

    def num(n):
        return f"{n:,}" if n is not None else "-"

    up, down = cat.relations(d)
    print(f"#{d['id']}  {d['name']}")
    for k, v in (
        ("path", d["path"]), ("source", d["source"]), ("platform", d["platform"]), ("format", d["format"]),
        ("modalities", d["modalities"]), ("organ", d["organ"]), ("species", d["species"]),
        ("cells", num(d["n_cells"])), ("features", num(d["n_features"])),
        ("study", d["study"]), ("tags", d["tags"]), ("description", d["description"]), ("owner", d["owner"]),
        ("steward", d["steward"]), ("size", size), ("added", f"{d['added_at']} by {d['added_by']}"),
        ("updated", d["updated_at"]), ("checked", d["checked_at"]), ("access", access),
        ("links", f"{len(links)} active ({mine} yours)"),
        ("checksum", f"{d['fingerprint']} ({d['fingerprinted_at']})" if d["fingerprint"] else "-"),
    ):
        print(f"  {k:<12}{opt(v)}")
    for kind, x in up:
        print(f"  {'comes from':<12}{kind} #{x['id']} {x['name']}")
    for kind, x in down:
        print(f"  {'used by':<12}#{x['id']} {x['name']} ({kind})")
    details = {k: v for k, v in (d.get("details") or {}).items() if k != "obs_columns"}
    if details:
        print(f"\n  inside the files (read {d['inspected_at']}):")
        for k, v in details.items():
            text = json.dumps(v) if isinstance(v, (dict, list)) else str(v)
            print(f"    {k.replace('_', ' '):<14}{text if len(text) < 110 else text[:107] + '...'}")
    print(f"\nGet a link with: stew get {d['name']} --to <folder>")


def cmd_info(a, ui, cat):
    d = cat.find(a.dataset)
    if a.json:
        up, down = cat.relations(d)
        d["relations"] = {"from": [{"kind": k, "name": x["name"]} for k, x in up],
                          "to": [{"kind": k, "name": x["name"]} for k, x in down]}
        print(json.dumps(d, indent=2))
    else:
        print_info(cat, d)


def pick(a, cat: Catalog) -> list[dict]:
    if a.all:
        return cat.list()
    if not a.datasets:
        raise UserError("name some datasets, or use --all")
    return [cat.find(k) for k in a.datasets]


def cmd_inspect(a, ui, cat):
    rows = []
    for d in pick(a, cat):
        p = Path(d["path"])
        if util.presence(p) is not Presence.READABLE:
            rows.append([str(d["id"]), d["name"], "-", "-", "-", "-", "can't read it"])
            continue
        links = link_provenance(cat, d)
        d = cat.by_id(d["id"])
        c = contents.read(p, d["format"])
        if c.empty():
            rows.append([str(d["id"]), d["name"], "-", "-", "-", "-",
                         ", ".join(links) or c.error or f"nothing to read in a {d['format']}"])
            continue
        details = {**c.details, **({"pipeline": d["details"]["pipeline"]} if (d.get("details") or {}).get("pipeline")
                                   else {})}
        ch = Changes(fields=[("n_cells", c.n_cells), ("n_features", c.n_features),
                             ("details", json.dumps(details) if details else None), ("inspected_at", now())])
        for col in ("organ", "species", "platform"):
            val = getattr(c, col)
            if val and (a.overwrite or not d[col]) and val != d[col]:
                ch.fields.append((col, val))
        known = {t.name for t in vocab.MODALITIES}
        ch.add_modalities = sorted((c.modalities & known) - set(d["modalities"]))
        cat.update(d, ch)
        updated = [k for k, _ in ch.fields if k in ("organ", "species", "platform")]
        note = ", ".join([f"set {k}" for k in updated] + [f"+{m}" for m in ch.add_modalities] + links) or "ok"
        rows.append([str(d["id"]), d["name"], f"{c.n_cells:,}" if c.n_cells else "-",
                     f"{c.n_features:,}" if c.n_features else "-", c.species or "-", c.organ or "-", note])
    table(["ID", "NAME", "CELLS", "FEATURES", "SPECIES", "ORGAN", "RESULT"], rows)


ACCESSION = re.compile(r"\b(GSE\d+|GSM\d+|[SED]R[PXR]\d+|PRJ(?:NA|EB|DB)\d+|E-[A-Z]{4}-\d+|idr\d{4}|syn\d{6,}|"
                       r"(?:zenodo|figshare)[:.]\d+)\b", re.I)


def own_accessions(d: dict) -> set[str]:
    """Accessions naming where a dataset itself comes from."""
    out = set()
    if d.get("source_key") and d["source_key"].split(":")[0] != "url":
        out.add(d["source_key"].split(":")[1].upper())
    return out


def suggest_relations(cat: Catalog) -> list[tuple[dict, dict, str, str]]:
    """(child, parent, kind, reason) implied by what datasets say about each other."""
    ds = cat.list(Filter(kinds=["dataset", "reference"]))
    by_acc: dict[str, list[dict]] = {}
    for d in ds:
        for acc in own_accessions(d):
            by_acc.setdefault(acc, []).append(d)
    existing = {(r["child_id"], r["parent_id"]) for r in cat.all_relations()}
    existing |= {(p, c) for c, p in existing}
    out, seen = [], set()

    def add(child, parent, kind, why):
        key = (child["id"], parent["id"])
        if child["id"] != parent["id"] and key not in existing and key not in seen:
            seen.add(key)
            out.append((child, parent, kind, why))

    for d in ds:
        mine = own_accessions(d)
        # a GEO sample whose study is a catalogued series
        if (d.get("source_key") or "").startswith("geo:GSM") and d.get("study"):
            for parent in by_acc.get(d["study"].upper(), []):
                add(d, parent, "part-of", f"{next(iter(mine), 'it')} is a sample of {d['study']}")
        text = " ".join(str(d.get(k) or "") for k in ("study", "description", "name"))
        for acc in {m.group(1).upper().replace(".", ":") for m in ACCESSION.finditer(text)} - mine:
            for parent in by_acc.get(acc, []):
                add(d, parent, "derived-from", f"its study/description mentions {acc}")
    for d in ds:  # 10x pipeline outputs not linked yet
        prov = provenance.read(Path(d["path"])) if Path(d["path"]).is_dir() else None
        for q in (prov.inputs if prov else []):
            up = cat.containing(q)
            if up:
                add(d, up, "derived-from", "named as input in its _cmdline")
        for q in (prov.references if prov else []):
            up = cat.containing(q)
            if up:
                add(d, up, "uses", "named as reference in its _cmdline")
    return out


def cmd_relate(a, ui, cat):
    if a.suggest:
        found = suggest_relations(cat)
        if not found:
            print("No new relations suggested.")
            return
        table(["CHILD", "KIND", "PARENT", "WHY"],
              [[f"#{c['id']} {c['name']}", k, f"#{p['id']} {p['name']}", why] for c, p, k, why in found])
        if a.apply:
            for c, p, k, _ in found:
                cat.relate(c, p, k)
            print(f"Recorded {len(found)} relation(s).")
        else:
            print("Record them with `stew relate --suggest --apply`, or one by one with `stew relate`.")
        return
    if not (a.child and a.parent):
        raise UserError("name the child and the parent dataset, or use --suggest")
    child, parent = cat.find(a.child), cat.find(a.parent)
    cat.relate(child, parent, a.kind)
    print(f"#{child['id']} {child['name']} {a.kind} #{parent['id']} {parent['name']}")


def cmd_unrelate(a, ui, cat):
    child, parent = cat.find(a.child), cat.find(a.parent)
    if not cat.unrelate(child, parent):
        raise UserError(f"{child['name']} and {parent['name']} are not related")
    print(f"Removed the relation between {child['name']} and {parent['name']}.")


class HashProgress:
    def __init__(self, name: str, total: int):
        self.name, self.total, self.done, self.last = name, total, 0, 0.0
        self.tty = sys.stderr.isatty()

    def __call__(self, n: int) -> None:
        self.done += n
        t = time.monotonic()
        if self.tty and t - self.last > 0.5:
            self.last = t
            pct = 100 * self.done / self.total if self.total else 100
            sys.stderr.write(f"\r  {self.name}: {pct:5.1f}% of {human_size(self.total)}".ljust(70))
            sys.stderr.flush()

    def end(self) -> None:
        if self.tty:
            sys.stderr.write("\r" + " " * 70 + "\r")


def readable_or_row(d: dict, rows: list, width: int) -> bool:
    state = util.presence(Path(d["path"]))
    if state is Presence.READABLE:
        return True
    rows.append([str(d["id"]), d["name"]] + ["-"] * (width - 3) + [state.value])
    return False


def cmd_checksum(a, ui, cat):
    rows = []
    for d in pick(a, cat):
        if not readable_or_row(d, rows, 5):
            continue
        prog = HashProgress(d["name"], d["size_bytes"] or 0)
        entries = checksum.manifest(Path(d["path"]), a.quick, a.jobs, prog)
        prog.end()
        old = cat.files(d["id"])
        fp = checksum.fingerprint(entries)
        if old and any(f["hash"].startswith("quick:") for f in old.values()) == a.quick:
            change = checksum.compare(old, entries).summary()
        else:
            change = "first checksum" if not old else "recomputed (" + ("quick" if a.quick else "full") + ")"
        cat.store_files(d, entries, fp)
        rows.append([str(d["id"]), d["name"], str(len(entries)), fp[:22] + "…", change])
    table(["ID", "NAME", "FILES", "FINGERPRINT", "CHANGES"], rows)


def cmd_verify(a, ui, cat):
    rows, bad = [], 0
    for d in pick(a, cat):
        old = cat.files(d["id"])
        if not old:
            rows.append([str(d["id"]), d["name"], "-", "no checksum yet (run `stew checksum`)"])
            continue
        if not readable_or_row(d, rows, 4):
            bad += 1
            continue
        quick = any(f["hash"].startswith("quick:") for f in old.values())
        prog = HashProgress(d["name"], d["size_bytes"] or 0)
        diff = checksum.compare(old, checksum.manifest(Path(d["path"]), quick, a.jobs, prog))
        prog.end()
        bad += not diff.ok
        rows.append([str(d["id"]), d["name"], "quick" if quick else "full", "ok" if diff.ok else diff.summary()])
    table(["ID", "NAME", "MODE", "RESULT"], rows)
    if bad:
        print(f"{bad} dataset(s) differ from their checksums.", file=sys.stderr)
        return 1


def cmd_duplicates(a, ui, cat):
    ds = cat.list()
    groups: dict[str, list[dict]] = {}
    for d in ds:
        if d["fingerprint"] and d["n_files"]:  # folders of links only hold no data of their own
            groups.setdefault(d["fingerprint"], []).append(d)
    same = [g for g in groups.values() if len(g) > 1]
    missing = sum(1 for d in ds if not d["fingerprint"])
    if same:
        print("Whole datasets stored more than once:")
        for g in same:
            print(f"  {human_size(g[0]['size_bytes'])} x {len(g)}: " + ", ".join(f"#{d['id']} {d['name']}" for d in g))
    files = [f for f in cat.duplicate_files() if f["size"] > 0]
    if files:
        wasted = sum(f["size"] * (f["n"] - 1) for f in files)
        print(f"\nFiles found in more than one dataset ({human_size(wasted)} could be saved):")
        for f in files[:25]:
            places = f["places"].split("\x1f")
            print(f"  {human_size(f['size'])} x {f['n']}: " + ", ".join(p.rstrip(":") for p in places[:4])
                  + (" ..." if len(places) > 4 else ""))
    if not same and not files:
        print("No duplicates among checksummed datasets.")
    if missing:
        print(f"\n{missing} dataset(s) have no checksum yet, so they weren't compared "
              "(`stew checksum --all --quick` is fast).")


def cmd_serve(a, ui, cat):
    from .web import serve
    serve(cat.path.resolve(), a.host, a.port)


def cmd_set(a, ui, cat):
    d = cat.find(a.dataset)
    ch = Changes()
    if a.name:
        check_name(a.name)
        ch.fields.append(("name", a.name))
    for col, norm in (("organ", vocab.organ), ("species", vocab.species), ("platform", vocab.platform),
                      ("study", str.strip), ("description", str.strip), ("steward", str.strip)):
        v = getattr(a, col)
        if v is not None:
            ch.fields.append((col, norm(v) if v.strip() else None))
    src = None
    if a.source is not None:
        src = sources.parse(a.source) if a.source.strip() else None
        ch.fields += [("source", src.display if src else None), ("source_key", src.key if src else None)]
    if a.modality:
        ch.set_modalities = modalities(a.modality)
    ch.add_modalities = modalities(a.add_modality)
    ch.rm_modalities = modalities(a.rm_modality)
    ch.add_tags, ch.rm_tags = tags(a.tag), tags(a.untag)
    if ch.is_empty():
        raise UserError("nothing to change; see `stew set --help`")
    warn_same_source(cat, src, exclude=d["id"])
    cat.update(d, ch)
    print(f"Updated #{d['id']} {d['name']}: {ch.describe()}")


def cmd_rm(a, ui, cat):
    d = cat.find(a.dataset)
    links = cat.links(dataset_id=d["id"])
    if links:
        print(f"Note: {len(links)} link(s) to it are still handed out (see `stew links {d['name']} --all`).")
    if not ui.yes:
        if not ui.tty:
            raise UserError(f"refusing to remove #{d['id']} {d['name']} without --yes")
        if not ui.confirm(f"Remove #{d['id']} {d['name']} from the catalogue? Its files are not touched.", False):
            return
    cat.delete(d)
    print(f"Removed #{d['id']} {d['name']} from the catalogue ({d['path']} is untouched).")


def cmd_get(a, ui, cat):
    d = cat.find(a.dataset)
    src = Path(d["path"])
    state = util.presence(src)
    if state is Presence.NO_ACCESS:
        raise UserError(f"you ({cat.user}) can't read {d['path']}. It is owned by {opt(d['owner'])}; "
                        f"ask {opt(d['steward'])} to give you access.")
    if state is Presence.MISSING:
        raise UserError(f"#{d['id']} {d['name']} points to {d['path']}, which no longer exists. "
                        f"Tell its steward ({opt(d['steward'])}).")
    try:
        folder = a.to.resolve(strict=True)
    except OSError:
        raise UserError(f"cannot open folder {a.to}")
    if not folder.is_dir():
        raise UserError(f"{folder} is not a folder")
    link = folder / (a.link_name or src.name)
    if link.is_symlink() and Path(os.readlink(link)) == src:
        print(f"{link} already links to {d['path']}")
        return
    if link.exists() or link.is_symlink():
        raise UserError(f"{link} already exists; pick another name with --as")
    link.symlink_to(src)
    print(f"{link} -> {d['path']}")
    try:
        cat.add_link(d, link)
    except sqlite3.Error as e:
        print(f"warning: the link was made but could not be recorded in the catalogue: {e}", file=sys.stderr)


def cmd_links(a, ui, cat):
    did = cat.find(a.dataset)["id"] if a.dataset else None
    if a.prune:
        gone = [lk for lk in cat.links(None if a.all else cat.user, did) if not Path(lk["link_path"]).is_symlink()]
        for lk in gone:
            cat.release_link(lk)
            print(f"Released {lk['link_path']} ({lk['dataset_name']}, {lk['user']}): the link was already deleted")
        print(f"{len(gone)} link record(s) pruned." if gone else "Nothing to prune: every recorded link exists.")
        return
    links = cat.links(None if a.all else cat.user, did, a.released)
    if not links:
        print("No links." if a.all else "No links of yours (use --all for everyone's).")
        return
    header = ["DATASET", "USER", "LINK", "CREATED", "STATE"] + (["RELEASED"] if a.released else [])
    rows = []
    for lk in links:
        p = Path(lk["link_path"])
        state = "released" if lk["released_at"] else "ok" if p.is_symlink() else "gone"
        rows.append([lk["dataset_name"], lk["user"], lk["link_path"], lk["created_at"], state]
                    + ([lk["released_at"]] if a.released else []))
    table(header, rows)


def cmd_release(a, ui, cat):
    try:
        link = util.link_path(a.link)
    except OSError:
        raise UserError(f"cannot find {a.link}")
    lk = cat.active_link(link)
    if not lk:
        raise UserError(f"{link} is not a link handed out by stew (see `stew links`)")
    if lk["user"] != cat.user:
        raise UserError(f"that link belongs to {lk['user']}")
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise UserError(f"{link} is no longer a symlink; not touching it")
    else:
        print("(the link was already deleted)")
    cat.release_link(lk)
    print(f"Released {link} ({lk['dataset_name']})")


def cmd_refresh(a, ui, cat):
    ds = [cat.find(a.dataset)] if a.dataset else cat.list()
    rows, missing = [], 0
    for d in ds:
        p = Path(d["path"])
        state = util.presence(p)
        if state is not Presence.READABLE:
            missing += state is Presence.MISSING
            rows.append([str(d["id"]), d["name"], "MISSING" if state is Presence.MISSING else "no access", "-"])
            continue
        size, n = util.disk_usage(p)
        owner = util.owner_of(p)
        status = ("changed" if (d["size_bytes"], d["n_files"]) != (size, n)
                  else "new owner" if d["owner"] != owner else "ok")
        cat.update_stats(d["id"], owner, size, n)
        rows.append([str(d["id"]), d["name"], status, human_size(size)])
    table(["ID", "NAME", "STATUS", "SIZE"], rows)
    if missing:
        print(f"{missing} dataset(s) are missing. Remove them with `stew rm <name>`, or fix their path.")


def cmd_summary(a, ui, cat):
    s = cat.summary()
    print(f"{s['datasets']} datasets, {human_size(s['bytes'])} in total\n")
    for label in ("kind", "modality", "organ", "species", "platform"):
        print(f"  {label:<10}" + (", ".join(f"{k} {n}" for k, n in s[label]) or "-"))


def cmd_config(a, ui):
    found = config.locate(a.catalog)
    if not found:
        raise UserError(config.NOT_FOUND)
    path, how = found
    try:
        s = config.load(path, how)
    except config.ConfigError as e:
        raise UserError(str(e))
    if a.get:
        val = {"root": s.root, "catalog": path, "move_on_add": str(s.move_on_add).lower(),
               "group": s.group or "", "layout": " ".join(f"{k}={v}" for k, v in s.layout.items())}[a.get]
        print(val)
        return 0 if path.exists() else 1
    print(f"catalogue   {path}{'' if path.exists() else '   (MISSING: run `stew init`)'}")
    print(f"found via   {how}")
    print(f"data root   {s.root}")
    print(f"settings    {path.parent / config.CONFIG_FILE}"
          + ("" if (path.parent / config.CONFIG_FILE).is_file() else "   (none: defaults)"))
    print(f"move_on_add {str(s.move_on_add).lower()}")
    print(f"group       {s.group or '-'}   folder_mode {oct(s.folder_mode)[2:]}")
    print("download    " + "   ".join(f"{k} {v}" for k, v in s.download.items()))
    for k, v in s.layout.items():
        print(f"layout      {k:<10} {v}")
    return 0 if path.exists() else 1


def cmd_vocab(a, ui):
    def show(title, terms):
        print(title)
        for t in terms:
            line = f"  {t.name:<16}"
            if t.about:
                line += t.about + ". "
            if t.implies:
                line += f"[{', '.join(t.implies)}] "
            if t.aliases:
                line += "aka " + ", ".join(t.aliases)
            print(line.rstrip())
        print()
    show("Modalities (fixed list):", vocab.MODALITIES)
    show("Organs (others are accepted as typed):", vocab.ORGANS)
    show("Species (others are accepted as typed):", vocab.SPECIES)
    show("Platforms [modalities they imply] (others are accepted as typed):", vocab.PLATFORMS)


def cmd_check(a, ui, cat):
    started = now()
    problems, changed, checked = [], [], 0
    for d in cat.list():
        p = Path(d["path"])
        state = util.presence(p)
        if state is not Presence.READABLE:
            problems.append(f"#{d['id']} {d['name']}: {state.value} ({d['path']})")
            continue
        size, n = util.disk_usage(p)
        owner = util.owner_of(p)
        if (d["size_bytes"], d["n_files"]) != (size, n):
            changed.append(f"#{d['id']} {d['name']}: size {human_size(d['size_bytes'])} -> {human_size(size)}, "
                           f"files {d['n_files']} -> {n}")
        cat.update_stats(d["id"], owner, size, n)
        old = cat.files(d["id"])
        if old:
            quick = any(f["hash"].startswith("quick:") for f in old.values())
            diff = checksum.compare(old, checksum.manifest(p, quick, a.jobs))
            checked += 1
            if not diff.ok:
                problems.append(f"#{d['id']} {d['name']}: {diff.summary()}")
    lines = [f"stew check {started} -> {now()}: {len(cat.list())} entries, {checked} checksummed"]
    lines += ["PROBLEM " + x for x in problems] + ["size changed " + x for x in changed]
    verdict = "OK" if not problems else f"{len(problems)} PROBLEM(S)"
    lines.append(f"RESULT: {verdict}" + (f"; {len(changed)} size change(s) without checksum mismatch"
                                          if changed else ""))
    text = "\n".join(lines)
    print(text)
    if a.report:
        a.report.parent.mkdir(parents=True, exist_ok=True)
        a.report.write_text(text + "\n")
    cat.log("check", None, verdict)
    return 1 if problems else 0


def newest_mtime(p: Path, limit: int = 5000) -> float:
    """When anything under a path last changed (stops after `limit` entries)."""
    try:
        latest = p.lstat().st_mtime
    except OSError:
        return 0.0
    if not p.is_dir() or p.is_symlink():
        return latest
    n = 0
    for root, dirs, files in os.walk(p):
        for f in files + dirs:
            try:
                latest = max(latest, os.lstat(os.path.join(root, f)).st_mtime)
            except OSError:
                pass
            n += 1
            if n >= limit:
                return latest
    return latest


def cmd_report(a, ui, cat):
    ds = cat.list()
    root = data_home(cat)
    total = sum(d["size_bytes"] or 0 for d in ds)

    def share(n):
        return f"{100 * n / total:4.1f}%" if total else "-"

    def group(key):
        out: dict[str, list] = {}
        for d in ds:
            k = key(d) or "(unset)"
            out.setdefault(k, [0, 0])
            out[k][0] += d["size_bytes"] or 0
            out[k][1] += 1
        return sorted(out.items(), key=lambda kv: -kv[1][0])

    def area(d):
        try:
            parts = Path(d["path"]).relative_to(root).parts
        except ValueError:
            return "(outside " + str(root) + ")"
        return "/".join(parts[:min(a.depth, max(1, len(parts) - 1))]) or "."

    print(f"Storage report for {root}  ({now()})")
    print(f"{len(ds)} catalogued entries, {human_size(total)}")
    try:
        du = shutil.disk_usage(root)
        print(f"Disk: {human_size(du.used)} used of {human_size(du.total)}, {human_size(du.free)} free "
              f"({100 * du.free / du.total:.0f}%)")
    except OSError:
        pass
    for title, key in (("By area", area), ("By owner", lambda d: d["owner"]),
                       ("By kind", lambda d: d["kind"]),
                       ("By first modality", lambda d: (d["modalities"] or [None])[0])):
        print(f"\n{title}:")
        table(["", "SIZE", "SHARE", "ENTRIES"], [[k, human_size(s), share(s), str(n)] for k, (s, n) in group(key)])
    print(f"\nLargest {a.top}:")
    big = sorted(ds, key=lambda d: -(d["size_bytes"] or 0))[:a.top]
    table(["ID", "NAME", "SIZE", "OWNER", "PATH"],
          [[str(d["id"]), d["name"], human_size(d["size_bytes"]), d["owner"], d["path"]] for d in big])
    linked = {lk["dataset_id"] for lk in cat.links(released=True)}
    cutoff = time.time() - a.stale_days * 86400
    stale = []
    for d in sorted(ds, key=lambda d: -(d["size_bytes"] or 0)):
        if d["id"] in linked or not d["size_bytes"]:
            continue
        mt = newest_mtime(Path(d["path"]))
        if mt and mt < cutoff:
            stale.append((d, mt))
        if len(stale) >= a.top:
            break
    print(f"\nNever linked with `stew get` and unchanged for {a.stale_days}+ days (largest {a.top}):")
    if stale:
        table(["ID", "NAME", "SIZE", "LAST CHANGED", "OWNER"],
              [[str(d["id"]), d["name"], human_size(d["size_bytes"]), time.strftime("%Y-%m-%d", time.localtime(mt)),
                d["owner"]] for d, mt in stale])
    else:
        print("  none")
    dups = [f for f in cat.duplicate_files() if f["size"] > 0]
    if dups:
        print(f"\nFiles stored in more than one dataset: {human_size(sum(f['size'] * (f['n'] - 1) for f in dups))}"
              " could be saved (see `stew duplicates`).")


def cmd_backup(a, ui, cat):
    folder = (a.to or cat.path.resolve().parent / "backups")
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / f"catalog-{time.strftime('%Y-%m-%d_%H%M%S')}.db"
    out = sqlite3.connect(dest)
    with out:
        cat.conn.backup(out)  # a consistent copy even while others are writing
    out.close()
    ok = sqlite3.connect(dest).execute("PRAGMA integrity_check").fetchone()[0]
    if ok != "ok":
        dest.unlink()
        raise UserError(f"the backup failed its integrity check: {ok}")
    old = sorted(folder.glob("catalog-*.db"))[:-a.keep] if a.keep > 0 else []
    for f in old:
        f.unlink()
    print(f"Backed up to {dest} ({human_size(dest.stat().st_size)}); {len(old)} old backup(s) removed.")


def cmd_log(a, ui, cat):
    table(["WHEN", "WHO", "ACTION", "DATASET", "DETAIL"],
          [[e["at"], e["user"], e["action"], e["dataset"], e["detail"]] for e in cat.events(a.limit)])


def cmd_downloads(a, ui, cat):
    rows = cat.downloads(limit=a.limit)
    if not rows:
        print("No downloads yet.")
        return
    for d in rows:  # notice downloads whose process has died
        if d["status"] == "running":
            cat.active_downloads(d["source_key"])
    out = []
    for d in cat.downloads(limit=a.limit):
        status = d["status"]
        if status == "running" and d["total_bytes"]:
            dest = Path(d["dest"])
            done = util.disk_usage(dest)[0] - pending_bytes(dest) if dest.exists() else 0
            status = f"running {min(100, 100 * done / d['total_bytes']):.0f}%"
        out.append([str(d["id"]), d["source"], d["user"], status, d["started_at"], d["finished_at"], d["dest"],
                    d["log"] or ""])
    table(["ID", "SOURCE", "USER", "STATUS", "STARTED", "FINISHED", "DEST", "LOG"], out)


def cmd_mv(a, ui, cat):
    """Moves a dataset, or a folder holding catalogued datasets, by renaming
    it, then updates the catalogue and repoints the links handed out."""
    try:
        d = cat.find(a.dataset)
        src = Path(d["path"])
    except CatalogError:
        p = Path(a.dataset)
        if "/" not in a.dataset or not p.is_dir():
            raise
        src = p.resolve()
    if not (src.exists() or src.is_symlink()):
        raise UserError(f"{src} doesn't exist; fix the path with `stew rm` and `stew add`")
    into = a.dest.endswith("/")
    dest = Path(os.path.abspath(os.path.expanduser(a.dest)))
    if into or (dest.is_dir() and not dest.is_symlink()):
        dest = dest / src.name
    if dest.exists() or dest.is_symlink():
        raise UserError(f"{dest} already exists")
    if dest == src or dest.is_relative_to(src):
        raise UserError("can't move a folder into itself")
    if not dest.parent.is_dir() and not a.parents:
        raise UserError(f"{dest.parent} doesn't exist; create it, or pass --parents")
    probe = dest.parent
    while not probe.exists():
        probe = probe.parent
    if os.stat(probe).st_dev != os.lstat(src).st_dev:
        raise UserError(f"{dest} is on a different filesystem from {src}. stew only moves by renaming; copy "
                        "the data yourself (e.g. rsync), check it, then `stew rm` + `stew add` the new copy.")
    affected = cat.under(str(src))
    if not affected:
        raise UserError(f"nothing catalogued at or below {src}; move it with plain `mv`")
    # rename(2) needs write access to both parents, and to a moved folder itself
    # when its parent changes (its ".." entry is rewritten).
    need = [src.parent, probe if not dest.parent.exists() else dest.parent]
    if src.is_dir() and not src.is_symlink() and dest.parent != src.parent:
        need.append(src)
    if a.leave_link:
        need.append(src.parent)
    blocked = sorted({str(p) for p in need if not os.access(p, os.W_OK)})
    if blocked:
        owners = {p: util.owner_of(Path(p)) for p in blocked}
        raise UserError("no permission to move it; you need write access to: "
                        + ", ".join(f"{p} (owner {owners[p]})" for p in blocked))
    moves = [(x, str(dest) + x["path"][len(str(src)):]) for x in affected]
    links = [(lk, new) for x, new in moves for lk in cat.links(dataset_id=x["id"])]
    print(f"{src}\n  -> {dest}\n")
    table(["ID", "NAME", "NEW PATH"], [[str(x["id"]), x["name"], new] for x, new in moves])
    others = sorted({lk["user"] for lk, _ in links if lk["user"] != cat.user})
    print(f"\n{len(moves)} catalogued dataset(s) move; {len(links)} link(s) handed out by `stew get` will be "
          "repointed" + (f" (some belong to {', '.join(others)})" if others else "") + ".")
    if a.leave_link:
        print(f"A symlink will be left at {src}.")
    if a.dry_run:
        print("\nDry run: nothing moved.")
        return
    if not ui.yes:
        if not ui.tty:
            raise UserError("refusing to move without --yes when not in a terminal")
        if not ui.confirm("Move?", False):
            return
    if a.parents:
        dest.parent.mkdir(parents=True, exist_ok=True)
    os.rename(src, dest)
    try:
        cat.move_paths(moves, f"{src} -> {dest}")
    except Exception:
        os.rename(dest, src)  # keep disk and catalogue in step
        raise
    if a.leave_link:
        src.symlink_to(dest)
    old_of = {x["id"]: x["path"] for x, _ in moves}
    fixed, failed = 0, []
    for lk, new in links:
        p = Path(lk["link_path"])
        try:
            if not p.is_symlink() or os.readlink(p) != old_of[lk["dataset_id"]]:
                continue
            tmp = p.with_name(p.name + ".stew-tmp")
            tmp.symlink_to(new)
            os.replace(tmp, p)
            fixed += 1
        except OSError:
            failed.append(lk)
    print(f"Moved. Updated {len(moves)} catalogue entr{'y' if len(moves) == 1 else 'ies'}, repointed {fixed} link(s).")
    if failed:
        how = "they still work through the symlink left behind" if a.leave_link else \
            "they now point nowhere until their owners run `stew get` again"
        print(f"Couldn't repoint {len(failed)} link(s) (no permission in their folders); {how}:")
        for lk in failed:
            print(f"  {lk['user']}: {lk['link_path']}")


# --- download ------------------------------------------------------------------------------

def download_root(cat: Catalog) -> Path | None:
    """An explicit download folder ($STEW_DOWNLOAD_DIR or `init --download-root`), if any."""
    env = os.environ.get("STEW_DOWNLOAD_DIR") or cat.get_meta("download_root")
    return Path(env) if env else None


def download_dest(cat: Catalog, src: sources.Source) -> Path:
    base = download_root(cat)
    if base:
        return base / src.kind / src.folder
    return layout_folder(cat, "public", kind=src.kind, accession=src.folder, name=src.folder)


def logs_dir(cat: Catalog) -> Path:
    return cat.path.resolve().parent / "logs" / "downloads"


MARKER = ".stew-download"


def describe(ds: list[dict]) -> str:
    return "\n".join(f"  #{d['id']} {d['name']}  {d['path']}" for d in ds)


def blocking_checks(cat: Catalog, src: sources.Source, extra_keys: list[str]) -> None:
    """Checks that always stop a download: the source (or a collection that
    contains it) is catalogued, or someone is downloading it now."""
    exact = cat.by_source([src.key])
    if exact:
        raise UserError(f"{src.display} is already downloaded and catalogued:\n{describe(exact)}\n"
                        f"Use `stew get {exact[0]['name']}` to link it into your folder.")
    containing = cat.by_source([k for k in src.related + extra_keys if k != src.key])
    if containing:
        raise UserError(f"{src.display} is part of (or the same as) {containing[0]['source']}, which is already "
                        f"here:\n{describe(containing)}\nUse `stew get {containing[0]['name']}`.")
    for key in [src.key] + [k for k in extra_keys if k != src.key]:
        running = cat.active_downloads(key)
        if running:
            r = running[0]
            raise UserError(f"{r['user']} has been downloading {r['source']} into {r['dest']} since "
                            f"{r['started_at']}. Wait for it to finish (see `stew downloads`).")


def mention_check(cat: Catalog, words: list[str], what: str, force: bool, hint: str = "") -> None:
    mention = list({d["id"]: d for w in words for d in cat.mentioning(w)}.values())
    if mention and not force:
        raise UserError(f"these datasets mention {what} in their study, name or path, so it may already be "
                        f"here:\n{describe(mention)}\n{hint}Otherwise re-run with --force.")


def default_name(src: sources.Source, remote: sources.Remote, target: Path) -> str:
    if src.kind == "url":
        return detect.default_name(target)
    if src.kind in ("cellxgene", "hca") and remote.title:
        base = util.slug(remote.title)[:40].strip("-")
    elif src.kind in ("zenodo", "figshare"):
        base = f"{src.kind}-{src.accession}"
    else:
        base = src.accession.lower()
    return base + ("-raw" if src.raw else "")


def stew_command() -> list[str]:
    exe = Path(sys.argv[0]).resolve() if sys.argv and sys.argv[0] else None
    if exe and exe.name == "stew" and exe.exists():
        return [str(exe)]
    return [sys.executable, "-m", "data_steward"]


def start_background(a, cat: Catalog, src: sources.Source, c: Candidate, dest: Path) -> None:
    """Re-runs this download in a tmux session (or a detached process), after
    the checks and questions have been answered here."""
    logs = logs_dir(cat)
    logs.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    log = logs / f"{util.slug(src.folder)}-{stamp}.log"
    argv = [t for t in a._argv if t != "--background"]
    argv += ["--yes", "--catalog", str(cat.path.resolve())]
    # answers given interactively here
    for flag, val in (("--name", c.name), ("--organ", c.organ), ("--species", c.species)):
        if val and flag not in argv:
            argv += [flag, val]
    for m in sorted(c.modalities):
        argv += ["-m", m]
    cmd = stew_command() + argv
    env = {**os.environ, "STEW_BACKGROUND_LOG": str(log)}
    session = f"stew-dl-{util.slug(src.folder)[:30]}-{stamp[-6:]}"
    if shutil.which("tmux"):
        # A running tmux server doesn't pass on our environment, so set it in the command.
        shell = f"STEW_BACKGROUND_LOG={shlex.quote(str(log))} {shlex.join(cmd)} >> {shlex.quote(str(log))} 2>&1"
        subprocess.run(["tmux", "new-session", "-d", "-s", session, shell], check=True, env=env)
        where = f"tmux session {session} (watch it with `tmux attach -t {session}`)"
    else:
        with open(log, "ab") as out:
            subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                             start_new_session=True, env=env)
        where = "a background process"
    print(f"\nDownloading {src.display} in {where}.\nLog: {log}\n"
          "`stew downloads` shows progress; stew will tell you when it has finished.")


def download_missing(a, ui, cat, src: sources.Source):
    """`download --missing`: tops up the catalogued dataset holding this source
    (or the record containing it) with the remote files its folder lacks."""
    given = [f"--{f.replace('_', '-')}" for f in ("name", "to", "organ", "species", "platform", "study",
                                                  "description", "steward", "derived_from") if getattr(a, f)]
    if given or a.modalities or a.tags:
        raise UserError("--missing adds files to an existing dataset and keeps its folder and description; "
                        "leave out " + ", ".join(given + ["-m"] * bool(a.modalities) + ["--tag"] * bool(a.tags))
                        + " and change those with `stew set`.")
    held = cat.by_source([src.key]) or cat.by_source([k for k in src.related if k != src.key])
    if not held:
        raise UserError(f"{src.display} isn't catalogued yet; download it without --missing.")
    if len(held) > 1:
        raise UserError(f"{src.display} is catalogued more than once; --missing doesn't know which to top up:\n"
                        f"{describe(held)}")
    d = held[0]
    dest = Path(d["path"])
    if not dest.is_dir():
        raise UserError(f"#{d['id']} {d['name']} isn't a folder ({dest}), so there is nothing to add files to.")
    keys = list(dict.fromkeys([d["source_key"], src.key]))

    def not_running():
        for key in keys:
            for r in cat.active_downloads(key):
                raise UserError(f"{r['user']} has been downloading {r['source']} into {r['dest']} since "
                                f"{r['started_at']}. Wait for it to finish (see `stew downloads`).")
    not_running()

    print(f"Looking up {src.display} ...")
    try:
        remote = sources.resolve(src)
    except sources.SourceError as e:
        raise UserError(str(e))
    if remote.note and not remote.files:
        raise UserError(remote.note)
    files = remote.files
    if a.files:
        import fnmatch
        files = [f for f in files if any(fnmatch.fnmatch(f.name, p) for p in a.files)]
        if not files:
            raise UserError("no files match " + ", ".join(a.files) + ". Available: "
                            + ", ".join(f.name for f in remote.files[:20]) + (" ..." if len(remote.files) > 20 else ""))
    missing, differ = [], []
    for f in files:
        here = dest / f.name
        if not here.exists():
            missing.append(f)
        elif f.size is not None and here.is_file() and here.stat().st_size != f.size:
            differ.append(f"{f.name} ({human_size(here.stat().st_size)} here, {human_size(f.size)} remote)")
    if differ:
        raise UserError(f"these files in {dest} differ in size from {src.display}, so stew won't touch them:\n  "
                        + "\n  ".join(differ) + "\nMove them aside first if they should be replaced.")
    if not missing:
        print(f"Nothing to add: all {len(files)} file(s) of {src.display} are in #{d['id']} {d['name']} ({dest}).")
        return

    total = sum(f.size or 0 for f in missing)
    print(f"\n{src.display}: {len(files) - len(missing)} of {len(files)} file(s) already in #{d['id']} {d['name']}; "
          "missing:")
    table(["FILE", "SIZE"], [[f.name, human_size(f.size)] for f in missing[:30]])
    if len(missing) > 30:
        print(f"... and {len(missing) - 30} more")
    print(f"\n{len(missing)} file(s), {human_size(total)} -> {dest}")
    if a.dry_run:
        print("\nDry run: nothing downloaded.")
        return
    try:
        check_space(dest, total)
    except FetchError as e:
        raise UserError(str(e))
    if ui.interactive:
        if not ui.confirm(f"Add them to {d['name']}?", True):
            return
    elif not ui.yes:
        raise UserError("not running in a terminal; re-run with --yes to download")
    if a.background:
        start_background(a, cat, src, argparse.Namespace(name=None, organ=None, species=None, modalities=set()), dest)
        return

    dl = cat.start_download(d["source_key"], src.display, dest, not_running)
    log = os.environ.get("STEW_BACKGROUND_LOG")
    cat.set_download_meta(dl, bool(log), log, total or None)
    old_umask = os.umask(0o002)
    try:
        dl_conf = settings(cat).download
        fetch_many(missing, dest, jobs=a.jobs or dl_conf["jobs"], connections=a.connections or dl_conf["connections"],
                   max_connections=dl_conf["max_connections"])
    except (FetchError, OSError, KeyboardInterrupt) as e:
        cat.finish_download(dl, "failed", detail=str(e) or type(e).__name__)
        if isinstance(e, KeyboardInterrupt):
            raise UserError("interrupted; run the same command again to resume")
        raise UserError(f"download failed: {e}\nRun the same command again to resume.")
    finally:
        os.umask(old_umask)

    size, n = util.disk_usage(dest)
    with cat.tx():
        cat.update_stats(d["id"], util.owner_of(dest), size, n)
        cat.log("download", d["name"], f"added {len(missing)} missing file(s) of {src.display}")
    cat.finish_download(dl, "done", d["id"], f"{len(missing)} missing file(s) added")
    print(f"\nAdded {len(missing)} file(s) to #{d['id']} {d['name']}: now {human_size(size)} in {n} file(s) "
          f"at {dest}")
    if d["fingerprint"]:
        print(f"Its checksum predates these files: run `stew checksum {d['name']}` to update it.")


def cmd_download(a, ui, cat):
    try:
        src = sources.parse(a.source, raw=a.raw)
    except sources.SourceError as e:
        raise UserError(str(e))
    if a.missing:
        return download_missing(a, ui, cat, src)
    meta = Meta(a)
    if meta.derived_from:
        cat.find(meta.derived_from)

    # 1. Is it already here? Only the catalogue and the disk are consulted.
    blocking_checks(cat, src, [])
    mention_check(cat, src.mention_words(), src.accession or src.display, a.force,
                  f"If it really is the same data, record that with `stew set <dataset> --source {src.display}`. ")
    dest = a.to.expanduser().resolve() if a.to else download_dest(cat, src)
    marker = dest / MARKER
    resuming = marker.exists() and marker.read_text().strip() == src.key
    if dest.exists() and any(dest.iterdir()) and not resuming:
        d = cat.by_path(dest)
        if d:
            raise UserError(f"{dest} is already catalogued as #{d['id']} {d['name']}")
        if not a.force:
            raise UserError(f"{dest} already exists and wasn't made by `stew download`. It may hold this data "
                            "already: check it and register it with `stew add`, pick another folder with --to, "
                            "or re-run with --force.")

    # 2. Ask the remote side what there is.
    print(f"Looking up {src.display} ...")
    try:
        remote = sources.resolve(src)
    except sources.SourceError as e:
        raise UserError(str(e))
    if remote.parents:
        blocking_checks(cat, src, remote.parents)
        mention_check(cat, [p.split(":")[1] for p in remote.parents],
                      ", ".join(p.split(":")[1] for p in remote.parents), a.force)
    if remote.note and not remote.files:
        raise UserError(remote.note)
    files = remote.files
    if a.files:
        import fnmatch
        files = [f for f in files if any(fnmatch.fnmatch(f.name, p) for p in a.files)]
        if not files:
            shown = ", ".join(f.name for f in remote.files[:20]) + (" ..." if len(remote.files) > 20 else "")
            raise UserError("no files match " + ", ".join(a.files) + ". Available: " + shown)

    # 3. Plan: what goes where, described how.
    single = src.kind == "url" and len(files) == 1
    target = dest / files[0].name if single else dest
    c = Candidate.bare(target, None, "download")
    if a.name:
        check_name(a.name)
        if cat.name_taken(a.name):
            raise UserError(f"the name '{a.name}' is already used")
    c.name = a.name or default_name(src, remote, target)
    c.species, c.organ, c.platform = remote.species, remote.organ, remote.platform
    c.modalities = set(remote.modalities)
    words = " ".join([remote.title or "", str(target)] + [f.name for f in files[:50]])
    c.infer_from(words)
    if not remote.species and len(vocab.all_species(vocab.tokens(words))) > 1:
        c.species = None  # e.g. "... of the human and mouse pancreas": don't pick one
    meta.apply(c)
    meta.study = meta.study or remote.study
    meta.description = meta.description or remote.title

    total = sum(f.size or 0 for f in files)
    unknown = sum(f.size is None for f in files)
    print(f"\n{src.display}: {remote.title or files[0].name}")
    shown = files if len(files) <= 30 else files[:25]
    table(["FILE", "SIZE"], [[f.name, human_size(f.size)] for f in shown])
    if len(shown) < len(files):
        print(f"... and {len(files) - len(shown)} more")
    print(f"\n{len(files)} file(s), {human_size(total)}" + (f" (+{unknown} of unknown size)" if unknown else "")
          + f" -> {dest}" + (" (resuming)" if resuming else ""))
    print(f"Will register as '{c.name}': platform {opt(c.platform)}, modalities {opt(sorted(c.modalities))}, "
          f"organ {opt(c.organ)}, species {opt(c.species)}, study {opt(meta.study)}")
    if remote.note:
        print(f"note: {remote.note}")
    if a.dry_run:
        print("\nDry run: nothing downloaded.")
        return
    try:
        check_space(dest, total)
    except FetchError as e:
        raise UserError(str(e))
    if ui.interactive:
        ask_metadata(ui, c, False)
        if not ui.confirm("Download?", True):
            return
    elif not ui.yes:
        raise UserError("not running in a terminal; re-run with --yes to download")
    if a.background:
        start_background(a, cat, src, c, dest)
        return

    # 4. Download, holding a record so nobody else starts the same one.
    dl = cat.start_download(src.key, src.display, dest, lambda: blocking_checks(cat, src, remote.parents))
    log = os.environ.get("STEW_BACKGROUND_LOG")
    cat.set_download_meta(dl, bool(log), log, total or None)
    old_umask = os.umask(0o002)  # downloads are shared: keep them group-writable
    try:
        shared_dirs(cat, dest.parent)
        dest.mkdir(exist_ok=True)
        marker.write_text(src.key + "\n")
        dl_conf = settings(cat).download
        fetch_many(files, dest, jobs=a.jobs or dl_conf["jobs"], connections=a.connections or dl_conf["connections"],
                   max_connections=dl_conf["max_connections"])
    except (FetchError, OSError, KeyboardInterrupt) as e:
        cat.finish_download(dl, "failed", detail=str(e) or type(e).__name__)
        if isinstance(e, KeyboardInterrupt):
            raise UserError("interrupted; run the same command again to resume")
        raise UserError(f"download failed: {e}\nRun the same command again to resume.")
    finally:
        os.umask(old_umask)

    # 5. Register it, read what is inside, and register recognisable parts.
    marker.unlink(missing_ok=True)
    found = detect.detect_any(target)
    c.path = target
    c.format = found.format if found.format not in ("directory", "file") else f"{src.kind}-download"
    c.platform = c.platform or found.platform
    c.modalities |= found.modalities
    c.read_contents()
    meta.apply(c)
    did, d = register(cat, c, meta, src)
    cat.finish_download(dl, "done", did, f"{len(files)} file(s)")
    print(f"\nDownloaded and registered #{did} {d['name']} ({human_size(d['size_bytes'])}) at {target}")
    parts = [p for p in detect.scan(target).found if p.path != target] if target.is_dir() else []
    if len(parts) > 1 or (parts and found.format == "directory"):
        whole = cat.by_id(did)
        part_meta = Meta(argparse.Namespace(organ=None, species=None, platform=None, modalities=[], study=meta.study,
                                            tags=meta.tags, description=None, steward=meta.steward,
                                            derived_from=None))
        added = []
        for p in parts:
            p.read_contents()
            p.species = p.species or c.species
            p.organ = p.organ or c.organ
            added.append(register(cat, p, part_meta, None, parent=whole, kind="part-of"))
        print(f"Also registered {len(added)} dataset(s) found inside it, as parts of {d['name']}: "
              + ", ".join(f"#{i} {x['name']}" for i, x in added))
    nag_missing([(did, d)])


# --- entry point --------------------------------------------------------------------------

NO_CATALOG = {"init": cmd_init, "vocab": cmd_vocab, "config": cmd_config}
COMMANDS = {
    "scan": cmd_add, "add": cmd_add, "ls": cmd_ls, "list": cmd_ls, "search": cmd_ls, "saved": cmd_saved,
    "info": cmd_info, "set": cmd_set, "rm": cmd_rm, "get": cmd_get, "links": cmd_links, "release": cmd_release,
    "refresh": cmd_refresh, "summary": cmd_summary, "log": cmd_log, "download": cmd_download,
    "downloads": cmd_downloads, "inspect": cmd_inspect, "relate": cmd_relate, "unrelate": cmd_unrelate,
    "checksum": cmd_checksum, "verify": cmd_verify, "duplicates": cmd_duplicates, "serve": cmd_serve,
    "mv": cmd_mv,
}


# Commands that change datasets or the shared data area. Everyone else may
# read, search, get/release links and keep their own saved queries.
STEWARD_ONLY = {"scan", "add", "set", "rm", "mv", "relate", "unrelate", "inspect", "checksum", "refresh",
                "download", "check"}


def stewards(cat: Catalog) -> list[str] | None:
    """The stewards, or None for an old catalogue without any (then anyone may edit)."""
    v = cat.get_meta("stewards")
    return [s for s in v.split(",") if s] if v else None


def check_steward(cat: Catalog, cmd: str) -> None:
    allowed = stewards(cat)
    if cmd in STEWARD_ONLY and allowed is not None and cat.user not in allowed:
        what = "downloads into the shared data area" if cmd == "download" else "changes to the catalogue"
        raise UserError(f"only stewards ({', '.join(allowed)}) can make {what}. "
                        f"Ask {allowed[0]} to `stew {cmd}` it, or use `stew get` to link existing data.")


def cmd_stewards(a, ui, cat):
    current = stewards(cat)
    if a.add or a.remove:
        if current is not None and cat.user not in current:
            raise UserError(f"only stewards ({', '.join(current)}) can change who the stewards are")
        new = [s for s in (current or [cat.user]) if s not in a.remove]
        new += [s for s in a.add if s not in new]
        if not new:
            raise UserError("there must be at least one steward")
        with cat.tx():
            cat.set_meta("stewards", ",".join(new))
            cat.log("stewards", None, ",".join(new))
        current = new
    if current is None:
        print("No stewards set: anyone who can write the catalogue file may change it. "
              "Set some with `stew stewards --add <user>`.")
    else:
        print("Stewards: " + ", ".join(current))


COMMANDS["stewards"] = cmd_stewards
COMMANDS["check"] = cmd_check
COMMANDS["backup"] = cmd_backup
COMMANDS["report"] = cmd_report


def notify(cat: Catalog) -> None:
    for d in cat.unnotified():
        if d["status"] == "done":
            ds = cat.by_id(d["dataset_id"]) if d["dataset_id"] else None
            what = f" as #{ds['id']} {ds['name']}" if ds else ""
            print(f"stew: your background download of {d['source']} finished{what}.", file=sys.stderr)
        else:
            print(f"stew: your background download of {d['source']} {d['status']}: {d['detail'] or ''}"
                  f" (log: {d['log']}). Re-run it to resume.", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    argv = sys.argv[1:] if argv is None else argv
    sys.stdout.reconfigure(line_buffering=True)  # keep output in order with progress on stderr
    if not argv:
        print(TOP_HELP, end="")
        return 0
    a = parser.parse_args(argv)
    if a.cmd == "help":
        sub = parser._subparsers._group_actions[0].choices
        if a.topic in sub:
            sub[a.topic].print_help()
        else:
            print(TOP_HELP, end="")
        return 0
    ui = Ui(a.yes)
    a._argv = argv[argv.index(a.cmd):] if a.cmd in argv else argv
    try:
        if a.cmd in NO_CATALOG:
            return NO_CATALOG[a.cmd](a, ui) or 0
        found = config.locate(a.catalog)
        if not found:
            raise UserError(config.NOT_FOUND)
        a.catalog = found[0]
        cat = Catalog.open(a.catalog)
        if not os.environ.get("STEW_BACKGROUND_LOG"):
            notify(cat)
        check_steward(cat, a.cmd)
        return COMMANDS[a.cmd](a, ui, cat) or 0
    except (UserError, CatalogError, vocab.VocabError, sources.SourceError, config.ConfigError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    except sqlite3.OperationalError as e:
        print(f"error: {e}", file=sys.stderr)
        if "readonly" in str(e) or "unable to open" in str(e):
            print(f"hint: you can read the catalogue but not change it ({a.catalog}). "
                  "Ask its steward for write access.", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
    except BrokenPipeError:
        # The reader went away (e.g. `stew ls | head`): stop quietly, without
        # Python complaining when it flushes stdout at exit.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
