# Changelog

## 0.6.0
- First public release.
- Configuration: `stew init <root>` writes `<root>/.data-steward/config.toml` (folder
  layout, whether added data is moved in, folder group and mode). stew finds the
  catalogue from `--catalog`, `$STEW_ROOT`, the current folder (like git),
  `~/.config/stew/config.toml`, `/etc/stew/config.toml`, or a built-in `site.toml`.
- `stew config` shows which catalogue is in use and its settings.
- Downloads follow the `public` layout.
- MIT licence; PyPI packaging (`pip install data-steward`, h5py optional).
