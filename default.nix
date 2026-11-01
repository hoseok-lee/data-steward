# Lets non-flake tools use this repo, e.g. `nix-env -f . -i`.
(builtins.getFlake (toString ./.)).packages.${builtins.currentSystem}.default
