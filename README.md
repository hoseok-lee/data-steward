# data-steward (`stew`)

[![test](https://github.com/hoseok-lee/data-steward/actions/workflows/test.yml/badge.svg)](https://github.com/hoseok-lee/data-steward/actions/workflows/test.yml)

A command-line catalogue for a lab's datasets: single-cell, multi-omics, spatial
transcriptomics, raw reads and slide images, kept in one shared data folder. `stew`
records what each dataset is (organ, species, modality, platform, cell counts read
from the files), where it came from (GEO, SRA, Zenodo, ...) and who owns it; answers
filtered queries; downloads public data only when it isn't already there; and hands
out symlinks so everyone works from one copy instead of many.

```console
$ stew ls --organ lung -m rna -m spatial
ID  NAME                         PLATFORM  MODALITIES   ORGAN  SPECIES  SIZE      OWNER
21  human-lung-cancer-ffpe-rep1  xenium    rna,spatial  lung   human    10.5 GiB  alice
22  human-lung-cancer-ffpe-rep2  xenium    rna,spatial  lung   human    10.5 GiB  alice

$ stew download GSE274673
error: GSE274673 is already downloaded and catalogued:
  #86 gse274673  /data/shared/public/geo/GSE274673
Use `stew get gse274673` to link it into your folder.

$ stew get human-lung-cancer-ffpe-rep1 --to ~/project/data
```

## Install

```sh
pipx install "data-steward[hdf5]"     # or: pip install "data-steward[hdf5]"
nix profile install github:hoseok-lee/data-steward      # or with Nix
```

Python 3.11+. The `hdf5` extra (h5py) lets stew read cell counts, organism and tissue
from inside `.h5ad`, `.h5mu`, `.loom` and 10x `.h5` files; without it stew still works.

## Quick start

```sh
stew init /data/shared --group lab --set-default   # a data folder with an empty catalogue
stew add ~/old/project/data --dry-run              # how stew would split and place it
stew add ~/old/project/data                        # register it (and move it in, leaving links)
stew download GSE123456                            # fetch public data, unless it is already here
stew ls -m spatial --species human                 # find
stew get gse123456 --to ~/analysis/inputs          # link, don't copy
```

## Configuration

`stew init <root>` creates `<root>/.data-steward/` holding `catalog.db` and `config.toml`:

```toml
move_on_add = true          # move added datasets into <root>, leaving a symlink behind
group = "lab"               # group for folders stew creates
folder_mode = "3775"        # group-writable, inherited group, only owners delete

[layout]                    # where datasets go, relative to <root>
public = "public/{kind}/{accession}"   # with a source, e.g. public/geo/GSE123456
reference = "references/{ref_kind}"    # genomes and annotations: 10x or annotations
project = "projects"                   # everything else
```

stew finds the catalogue, in order, from: `--catalog`; `$STEW_ROOT` (or
`$DATA_STEWARD_CATALOG`); the current folder or one above it holding `.data-steward/`
(like git); `~/.config/stew/config.toml`; `/etc/stew/config.toml`; and a `site.toml`
built into the package. The config files hold one line, `root = "/data/shared"`.
For a shared server, an admin writes `/etc/stew/config.toml` once and everyone finds the
catalogue; or build the package with a `site.toml` next to `flake.nix`. `stew config`
shows which catalogue is in use and how it was found.

## Commands

`stew` on its own, `stew --help` or `stew help <command>` show help with examples.
A dataset can be named by its name, its numeric id, or its path.

| command | what it does |
|---|---|
| `ls` (`list`) | list datasets; filter with `--organ`, `--species`, `-m/--modality`, `--platform`, `--format`, `--tag`, `--owner`, `--study`, `--source`, `--name`, `--under`, `--text`; `-l` for cells, path, source and tags; `--paths`, `--json` or `--csv` for scripts; `--save NAME [--shared]` to keep the query, `-q NAME` to rerun it |
| `search <words>` | datasets whose name, description, study, source, path or file contents (cell types, tissues, ...) contain all the words; takes the `ls` filters too |
| `saved [--rm NAME]` | your saved queries and shared ones |
| `info <dataset>` | everything known about one dataset: contents, relations, checksum, whether you can read it |
| `serve [--port 8765]` | read-only web page of the catalogue (see below) |
| `get <dataset> [--to DIR] [--as NAME]` | make a symlink to the dataset in your folder and record it |
| `links [--all] [--released] [--prune]` | links handed out (yours by default); `--prune` releases links whose symlink was deleted by hand |
| `release <link>` | delete a link made by `get` |
| `download <source>` | download from GEO, SRA/ENA, CELLxGENE, ArrayExpress, HCA, figshare, Zenodo or a URL, unless it is already here; `--background` runs it in tmux (see below) |
| `downloads` | downloads in progress, finished and failed |
| `add <path>` (`scan`) | register a file, a folder stew recognises as one dataset, or every dataset found inside any other folder (scripts and loose files stay; nothing found = the folder as one); moved into the data folder with a symlink left behind. `--whole`, `--each`, `--to DIR`, `--in-place`, `--dry-run` |
| `set <dataset> ...` | edit metadata: `--organ`, `--species`, `-m` (replace), `--add-modality`, `--rm-modality`, `--tag`, `--untag`, `--study`, `--source`, `--desc`, `--steward`, `--name` |
| `rm <dataset>` | remove from the catalogue; files are never touched |
| `inspect <dataset…> \| --all` | read cell/feature counts, organism, tissue, assay and modalities from inside the files |
| `relate <child> <parent> [--kind]`, `relate --suggest [--apply]` | record that a dataset is `derived-from` (default), `part-of` or the `same-as` another; `unrelate` removes it |
| `checksum <dataset…> \| --all [--quick]` | record a hash of every file; `verify` re-checks them; `duplicates` lists datasets and files stored more than once |
| `mv <dataset\|folder> <dest>` | move on disk (a rename on the same filesystem), update the catalogue, repoint links; `--leave-link` keeps the old path working; `dest/` means into that folder |
| `refresh [dataset]` | re-check paths, sizes and owners; flags missing datasets |
| `report` | storage report: space by area, owner and modality; largest datasets; never-linked, unchanged ones |
| `check` / `backup` | nightly health check (present, checksums match) and catalogue backup |
| `summary` | counts by modality, organ, species and platform |
| `vocab` | the controlled vocabulary and accepted aliases |
| `log` | recent changes, with who made them |
| `init [--group lab]` | create a new, empty catalogue |

### Filters

Repeating a filter for organ, species, platform or format means *any of*
(`--organ lung --organ liver`); repeating a modality or tag means *all of*
(`-m rna -m spatial` finds spatial transcriptomics). Aliases are understood:
`--organ lungs`, `--organ pdac` (pancreas), `-m scrna` (rna), `--species mm10`
(mouse). See `stew vocab`.

### Sources

Every dataset can record where it came from: a GEO accession (`GSE…`, `GSM…`), a
Zenodo record (`zenodo:123456`) or any URL. Give it with `--source` to `add`, `scan`
or `set`; `download` records it automatically. Different spellings of the same
source are recognised as one: `GSE274673`, its GEO page URL and its FTP folder are
the same, as are `zenodo:123`, `10.5281/zenodo.123`, `https://doi.org/10.5281/zenodo.123`
and `https://zenodo.org/records/123`. `stew ls --source GSE274673` finds it.

### Downloading

```sh
stew download GSE274673 --dry-run        # checks, file list and plan; nothing fetched
stew download GSM8454446                 # one GEO sample's processed files
stew download GSE84133 --raw             # its FASTQs, via SRA/ENA
stew download SRP000001 --files 'SRR000003_*'
stew download E-MTAB-5061                # ArrayExpress / BioStudies
stew download cellxgene:6f6d381a-7701-4781-935c-db10d30de293 --files '*core*'
stew download hca:74b6d569-3b11-42ef-b6b1-a0454522b4a0        # add --raw for FASTQs
stew download figshare:5968960 --files metadata_droplet.csv
stew download zenodo:5979761 --organ lung --tag tracking --background
stew download https://example.org/atlas.h5ad --name my-atlas --derived-from hlca-core
```

| source | how stew gets the files |
|---|---|
| GEO `GSE…`/`GSM…` | supplementary (processed) files; `--raw` follows GEO's SRA/BioProject link to the FASTQs |
| SRA/ENA `SRP/SRX/SRR/PRJNA/ERP/PRJEB…` | FASTQs from ENA's mirror, md5-checked (submitted BAMs if there are no FASTQs); runs ENA doesn't mirror come as `.sra` files from NCBI's AWS open-data copy (convert with `fasterq-dump`) |
| ArrayExpress/BioStudies `E-MTAB-…` | the study's files (processed data, SDRF/IDF); raw reads via their ENA accession |
| CELLxGENE `cellxgene:<collection>` | one H5AD per dataset in the collection |
| HCA `hca:<project>` | processed files (matrices, loom, h5ad); `--raw` for FASTQs; sha256-checked |
| figshare, Zenodo | all files of the record, md5-checked |
| any URL | the one file |

Before fetching anything, `download` checks that the data isn't already here:

1. a catalogued dataset already has this source, or the GEO series / Zenodo record
   that contains it (a GSM's series, a file URL's record), or the SRA study that a
   GEO series' raw reads point to;
2. a catalogued dataset mentions the accession (or a containing series) in its
   study, name or path, which usually means it was registered before sources were
   recorded;
3. someone is downloading it right now;
4. the destination folder already exists without being a stew download.

Checks 1 and 3 always stop the download; 2 and 4 can be overridden with `--force`.
The first round of checks uses only the catalogue and the disk, so a repeat
request never reaches the network.

Files go to `<data root>/downloads/<kind>/<id>/` (the data root is the folder
holding `.data-steward/`; override with `--to` or `$STEW_DOWNLOAD_DIR`). Downloads
resume if interrupted (re-run the same command), published checksums are verified,
and free space is checked first. Finished downloads are registered with their
source, title and study, the files are read for cell counts, organism and tissue,
and recognisable datasets inside a download (e.g. each H5AD of a CELLxGENE
collection) are registered too, as `part-of` the download.

With `--background`, stew does the checks and asks its questions, then hands the
download to a tmux session (`stew-dl-…`) and returns. `stew downloads` shows
progress and the log file, and your next `stew` command tells you when it finished
or failed.

### References, project folders and provenance

Genome references are catalogued as `kind = reference`, separate from datasets:
10x reference packages (a folder with `reference.json` and `fasta/`) and folders
holding only annotation files (`.gtf`, `.bed`, `.chrom.sizes`, `.fa`, ...).
`stew ls` shows datasets only and says how many references matched;
`stew ls --kind reference` (or `--kind all`) lists them.

Some folders are one dataset as a whole even though stew would find pieces
inside, such as an analysis project of R objects and scripts. Either register
each subfolder whole with `stew add <parent> --each`, or put a `.stew-dataset`
file in the folder; it may hold JSON metadata (`name`, `organ`, `species`,
`modalities`, `study`, `source`, `description`, `tags`, `steward`, `kind`) that
stew uses when it registers the folder. Either way, later scans skip what is
inside a catalogued dataset.

10x pipeline outputs (cellranger, cellranger-arc, cellranger-atac, spaceranger)
record the command they ran in `<id>/_cmdline`. stew reads it, follows
`--fastqs`, `--libraries` CSVs and `--reference`, and records `derived-from`
relations to the catalogued FASTQ folders and a `uses` relation to the reference.
FASTQ folders then take the species and organ of what was made from them.

cellranger-atac outputs and folders of scATAC fragments (`fragments.tsv.gz` with
`barcodes.tsv`) are recognised, with cell counts; for barnyard references
(GRCh38-and-mm10) the species is set only when at least 90% of cells map to one genome.

### Moving data

`stew mv <dataset or folder> <destination>` renames on disk (same filesystem only,
so it is instant), updates the path of every catalogued dataset at or below the
old path, and repoints the symlinks handed out by `stew get`. `--leave-link`
leaves a symlink at the old path so scripts that use it keep working; `--dry-run`
shows the plan. A destination ending in `/` is a folder to move into (`--parents`
creates it). Moving a messy existing folder into a layout is a sequence of these,
with a link left at every old name so nothing that uses the old paths breaks.

### What is inside a dataset

`scan`, `add` and `download` read the files themselves, and `stew inspect` does it
again later. From `.h5ad`, `.h5mu`, `.loom`, 10x `.h5`, AnnData/SpatialData `.zarr`,
Xenium bundles (`experiment.xenium`) and Cell/Space Ranger outputs
(`metrics_summary.csv`) stew records the number of cells (or spots) and features,
and where present the organism, tissues, assay, diseases and cell types
(CELLxGENE-style `obs` columns), feature types (gene expression, antibody capture,
peaks), obsm keys and layers. Species falls back to the gene names (`ACTB` vs
`Actb`, `ENSG…` vs `ENSMUSG…`) and the reference genome. Only headers, indices and
small annotation columns are read, never expression matrices. File contents beat
guesses from folder names; values you set yourself are kept unless you run
`stew inspect --overwrite`.

### Relations, checksums, duplicates

`stew relate filtered raw` records that `filtered` was derived from `raw`; `--kind
part-of` and `--kind same-as` are the other kinds. `stew info` shows both
directions, and loops are refused.

`stew checksum` stores every file's size, mtime and BLAKE2b hash plus a fingerprint
of the whole dataset; `--quick` hashes only the first and last MiB of each file,
which is enough to find copies in seconds. `stew verify` re-hashes and reports
changed, missing and new files (exit status 1 if anything differs). `stew
duplicates` lists whole datasets with the same fingerprint and individual files
stored in more than one dataset, with the space that could be saved.

### Web page

`stew serve` serves a read-only page: a sortable table with search, organ/species/
modality/platform filters, a detail panel (contents, relations, checksum, the `stew
get` command) and CSV export. It opens the catalogue read-only and never serves data
files. By default it listens on the server's localhost only; people reach it through
an SSH tunnel from their own computer, which needs no root and no open port:

```sh
stew serve --port 8765                                 # on the server, e.g. in tmux
ssh -N -L 8765:localhost:8765 <you>@<server>           # on your computer; open http://localhost:8765/
```

To serve it on the network instead, use `--host 0.0.0.0` behind your firewall, or put
it behind a reverse proxy with a login.

### What `scan` recognises

| layout | platform / format |
|---|---|
| folder with `experiment.xenium` (organism and tissue read from it) | xenium / xenium-bundle |
| `outs/spatial/tissue_positions*` | visium / spaceranger |
| `outs/binned_outputs/` | visium-hd / spaceranger |
| `*exprMat_file.csv` | cosmx |
| `cell_by_gene.csv` + `detected_transcripts.csv` | merscope |
| `outs/filtered_feature_bc_matrix*` (+ `atac_fragments.tsv.gz`) | 10x-chromium / cellranger (10x-multiome / cellranger-arc) |
| `*.zarr` | spatialdata, anndata-zarr or zarr |
| `*.h5ad`, `*.h5mu`, `*.loom`, `*.rds`, `*.h5seurat`, `*feature_bc_matrix.h5` | analysis objects |
| folders of `*.fastq.gz` / `*.bam` | fastq / bam (one dataset per folder) |
| `*.svs`, `*.ndpi`, `*.mrxs`, `*.ome.tif` | slide images |

Organ, species, platform and extra modalities are guessed from words in the path
(`human_lung_cancer`, `mouse_brain`, `pdac`, `visium_hd`, ...). Anything not
guessed is shown as `?` before registering, and can be filled in during review or
later with `set`.

## The shared catalogue

The catalogue is a single SQLite file in `<root>/.data-steward/` (see
[Configuration](#configuration) for how stew finds it).

`init --download-root DIR` sets a separate folder for downloads (by default they follow the `public` layout).
`init --group lab` makes the catalogue folder group-owned and group-writable, so lab
members can register datasets, download, and record their links; everyone else can
still read it. Every change is logged with the user who made it (`log`). Newer
versions of stew upgrade an older catalogue in place the first time they open it.

Optionally, stewards (set with `init --steward` or `stew stewards --add/--remove`) are
the only users who may register, edit, remove, move, inspect, checksum or download
datasets; everyone else can list, search, `get`/`release` links and keep their own
saved queries. Without stewards (the default) anyone who can write the catalogue may change it.

"Access" here means tracking, not granting: `get` only links datasets you can
already read. If you can't, it tells you who owns the data and whom to ask.

## Nightly checks

`stew backup` copies the catalogue safely (even while in use) to
`.data-steward/backups/`, keeping the newest 30; `stew check` confirms every dataset is
present and unchanged since its checksum and ends with `RESULT: OK` or the problems.
`scripts/nightly.sh` runs both and logs to `.data-steward/logs/`. Schedule it with cron
(`0 3 * * * /path/to/nightly.sh`), a systemd timer, or, where neither is available,
`scripts/scheduler.sh` in tmux. `stew report` summarises storage: space by area, owner
and modality, the largest datasets, and the ones nobody links to.

For lab members, [docs/USER_GUIDE.md](docs/USER_GUIDE.md) is a one-page guide you can
adapt (replace `/data/shared` and the steward's name).

## Development

```sh
nix develop                     # python 3.12 with pytest, h5py, numpy; ruff
python -m pytest                # tests (also run by `nix build` and CI)
ruff check src tests scripts
python -m data_steward --help   # run from the checkout
nix build                       # result/bin/stew
python scripts/make-dummy-data.py /tmp/dummy   # fake tree: real small HDF5 files + sparse big ones,
                                               # plus mirror/: a typical lab data folder in miniature
```

Releases: bump `version` in `pyproject.toml` (and `__init__.py`), add a CHANGELOG entry, then
create a GitHub release tagged `v<version>`; `.github/workflows/publish.yml` tests, builds and
uploads it to PyPI.

The flake is pinned to nixpkgs 25.05 so it also evaluates on old Nix releases (2.11).
