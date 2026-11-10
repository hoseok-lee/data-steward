# Changelog

## 0.7.1
- `stew download <source> --missing` tops up a catalogued dataset with the remote files
  its folder lacks (e.g. ones left out with `--files`), into the same folder, and
  updates the dataset's size and file count. It works with `--files`, `--dry-run` and
  `--background`, and leaves files that are there with a different size alone (and
  says so). Without `--missing`, a catalogued source is still refused.

## 0.7.0
- `stew download` fetches large files (64 MiB and up) over several connections at once
  when the server accepts ranged requests: 8 per file by default (`--connections N`;
  `--connections 1` is the old single stream). Connections that finish early take over
  half of the largest part still running. On Zenodo this is several times faster.
- Interrupted downloads resume each part where it stopped (`<file>.part.json` beside the
  `<file>.part` lists what is missing). A `.part` left by an earlier stew is resumed, not
  restarted: its bytes count as done and only the rest is split.
- A dataset's files are downloaded 3 at a time (`--jobs N`), with at most 16
  connections altogether. Progress: one line on a terminal; in logs, a line per file
  every 30 s. One failed file no longer stops the others.
- Servers without ranged requests, files of unknown size and signed URLs (HCA) use one
  connection. A "busy" answer (429/503) is waited out, honouring Retry-After, and the
  file continues over fewer connections instead of failing.
- New `[download]` settings in `config.toml`: `connections`, `jobs`, `max_connections`.
- `stew downloads` no longer counts the not-yet-written gaps of a `.part` as done.

## 0.6.1
- `stew download` of a BioStudies / BioImage Archive study (e.g. `S-BIAD1759`) now gets
  every file, including those in subfolders, from the BioStudies file list (paths and
  sizes kept); before, only the top-level metadata files were found. Studies without a
  file list are read folder by folder. BioImage Archive titles are picked up.
- `stew download <url>` of a web server's folder listing ("Index of …") downloads the
  folder's files, subfolders included, instead of saving the listing page. Links to the
  parent folder, sort orders and other sites are not followed.

## 0.6.0
- First public release.
- Configuration: `stew init <root>` writes `<root>/.data-steward/config.toml` (folder
  layout, whether added data is moved in, folder group and mode). stew finds the
  catalogue from `--catalog`, `$STEW_ROOT`, the current folder (like git),
  `~/.config/stew/config.toml`, `/etc/stew/config.toml`, or a built-in `site.toml`.
- `stew config` shows which catalogue is in use and its settings.
- Downloads follow the `public` layout.
- MIT licence; PyPI packaging (`pip install data-steward`, h5py optional).
