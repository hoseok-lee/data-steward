# Changelog

## 0.6.1
- `stew download` of a BioStudies / BioImage Archive study (e.g. `S-BIAD1759`) now gets
  every file, including those in subfolders, from the BioStudies file list (paths and
  sizes kept); before, only the top-level metadata files were found. Studies without a
  file list are read folder by folder.
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
