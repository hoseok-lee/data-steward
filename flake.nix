{
  description = "data-steward (stew): a command-line catalogue for a lab's datasets";

  # nixos-25.05 also evaluates on old Nix releases (2.11), which some lab servers still run.
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAll = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
      pyproject = builtins.fromTOML (builtins.readFile ./pyproject.toml);
    in
    {
      packages = forAll (pkgs: rec {
        data-steward = pkgs.python3Packages.buildPythonApplication {
          pname = "data-steward";
          version = pyproject.project.version;
          pyproject = true;
          src = pkgs.lib.cleanSourceWith {
            src = ./.;
            filter = path: type:
              let base = baseNameOf path; in
              base != "__pycache__" && (builtins.elem base [ "src" "tests" "scripts" "pyproject.toml" "README.md" "LICENSE" "site.toml" ]
                || pkgs.lib.hasInfix "/src/" path || pkgs.lib.hasInfix "/tests/" path
                || pkgs.lib.hasInfix "/scripts/" path);
          };
          build-system = [ pkgs.python3Packages.setuptools ];
          # h5py lets stew read inside .h5ad/.h5mu/.loom/10x .h5 files.
          dependencies = [ pkgs.python3Packages.h5py ];
          nativeCheckInputs = [ pkgs.python3Packages.pytestCheckHook pkgs.python3Packages.numpy ];
          # A site.toml next to this file (`root = "/path/to/data"`) becomes the built-in default
          # catalogue, so everyone who installs from that checkout finds it without any setup.
          postPatch = ''
            if [ -f site.toml ]; then cp site.toml src/data_steward/site.toml; fi
          '';
          meta = {
            description = pyproject.project.description;
            homepage = "https://github.com/hoseok-lee/data-steward";
            license = pkgs.lib.licenses.mit;
            mainProgram = "stew";
          };
        };
        default = data-steward;
      });

      devShells = forAll (pkgs: {
        default = pkgs.mkShell {
          packages = [ (pkgs.python3.withPackages (ps: [ ps.pytest ps.h5py ps.numpy ])) pkgs.ruff ];
          shellHook = ''export PYTHONPATH="$PWD/src''${PYTHONPATH:+:$PYTHONPATH}"'';
        };
      });
    };
}
