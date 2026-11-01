# Finding and using lab data with `stew`

<!-- A one-page guide for lab members. Replace /data/shared, <server> and <steward>. -->

The lab's shared datasets live in **`/data/shared`**, and `stew` is the catalogue
of them: what each dataset is (organ, species, modality, platform, cell counts),
where it came from (GEO, Zenodo, ...), who looks after it, and where it is on disk.
Use it to find data and to link it into your project instead of copying it.

## Install (once)

```sh
pipx install "data-steward[hdf5]"     # or however your lab installs it
stew            # shows the help
```

Run the same command again to update. It works from any folder.

## Find data

```sh
stew ls                                   # everything
stew ls --organ lung -m rna               # lung, with RNA
stew ls -m rna -m spatial --species human # human spatial transcriptomics
stew ls --platform xenium -l              # -l adds cells, study, source and path
stew search fibroblast pdac               # words in names, descriptions, cell types, ...
stew ls --source GSE274673                # where a GEO series (or Zenodo, URL) ended up
stew info gse274673                       # everything about one dataset
```

Aliases work (`--organ lungs`, `--organ pdac`, `-m scrna`, `--species mm10`); see
`stew vocab`. Repeating `--organ` means *any of*; repeating `-m` means *all of*.
`stew ls --kind reference` lists genome references and annotations (10x refdata,
gene sets, chromosome sizes).

## Use data

```sh
cd ~/my-project/data
stew get gse274673                        # a symlink to the shared copy, here
stew get 21 --as lung1 --to ../inputs     # by id, with your own name and folder
stew links                                # the links you have made
stew release lung1                        # remove one
stew links --prune                        # forget links you deleted by hand
```

**Please don't copy datasets** into your own folders; link them. Treat
`/data/shared` as read-only: write your results in your own project folder.

## Add data you already have

All datasets live in `/data/shared`. `stew add` **moves your data there and leaves
a symlink at the old path**, so your scripts and notebooks keep working.

```sh
stew add /data/me/lung_atlas.h5ad --organ lung --species human --dry-run  # preview
stew add /data/me/lung_atlas.h5ad --organ lung --species human -m rna \
  --desc "what it is, in one line"
stew add /data/me/GSE123456 --source GSE123456      # -> public/geo/GSE123456
stew add /data/me/project/data --dry-run            # a folder: see how it is split
stew set lung-atlas --tag pilot                          # fix or add details later
```

- **What gets registered** depends on what you point at: a **file** is one dataset; a
  folder stew recognises as **one dataset** stays together (Xenium, Visium, Cell Ranger
  incl. ARC/ATAC, CosMx, MERSCOPE, `.zarr`, `matrix.mtx` + barcodes + features, a folder
  of count/expression tables, fragments, FASTQ/BAM); **any other folder** is searched
  and each dataset inside (including each `.h5ad`, `.h5mu`, `.loom`, `.rds`, slide
  image) is registered separately, while scripts, logs and other files stay where they
  are (a folder with nothing recognisable inside is added as one dataset). Run `--dry-run` first. `--whole` keeps a folder as one dataset; `--each` makes
  each subfolder one dataset; a `.stew-dataset` file in a folder also keeps it whole.
- **Where it goes:** `public/<geo|zenodo|...>/<accession>/` if you give a `--source`,
  `references/` for genomes and annotations, otherwise `projects/<name>/`. Choose
  yourself with `--to projects/my_study`; `--in-place` registers without moving.
- stew reads cell counts, organism and tissue from the files; add what it couldn't
  work out with `--organ`, `--species`, `-m`, `--source`. Moves only work on
  the same disk as `/data/shared`. You stay the owner of what you add, and nobody
  else can delete or rename it.

## Get new data

Public data is downloaded once, into `/data/shared/public/`, by a steward
(<steward>). Before asking, check it isn't already here:

```sh
stew download GSE123456 --dry-run         # says "already here", or lists what it would fetch
```

GEO (`--raw` for the reads), SRA/ENA, ArrayExpress, CELLxGENE, the Human Cell Atlas,
figshare and Zenodo accessions or plain URLs all work. Then ask the steward to run it.

## Good to know

- **Web page.** For browsing without the terminal:
  `ssh -N -L 8765:localhost:8765 <you>@<server>`, then open http://localhost:8765/.
- **Old paths still work.** Datasets moved into `/data/shared` leave a symlink at
  their old path.
- **Layout.** `projects/` (lab and collaborator data), `public/<source>/<accession>/`
  (published data), `references/` (genomes, annotations).
- **"You can NOT read it".** `stew info` and `stew get` say who to ask for access.
- **Something missing or wrong?** Tell the steward. `stew log` shows recent changes,
  and the catalogue is backed up and checked every night.
