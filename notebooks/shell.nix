{ pkgs ? import <nixpkgs> {}}:
let
  fhs =  pkgs.buildFHSEnv rec{
    name = "Statistics_in_structural_bioinformatics";

    targetPkgs = _: [
      pkgs.jupyter
      pkgs.micromamba
    ];

    profile = ''
      set -e
      export MAMBA_ROOT_PREFIX=${builtins.getEnv "HOME"}/.mamba
      # Check if the environment already exists
      if ! micromamba env list | grep -q "$MAMBA_ROOT_PREFIX/envs/${name}"; then
        echo "Creating environment '${name}'..."
        micromamba create -q -y -n ${name}
      else
        echo "Environment '${name}' already exists. Skipping creation."
      fi

      echo "$MAMBA_ROOT_PREFIX/envs/${name}"
      chmod u+w "$MAMBA_ROOT_PREFIX/envs/${name}"

      # Install packages only if env.yml has been updated
      if [ ! -f $MAMBA_ROOT_PREFIX/envs/${name}/.env-hash ] || ! sha256sum -c $MAMBA_ROOT_PREFIX/envs/${name}/.env-hash; then
        echo "Installing/updating packages from env.yml..."
        micromamba install -y -n ${name} -f env.yml -c conda-forge
        sha256sum env.yml > $MAMBA_ROOT_PREFIX/envs/${name}/.env-hash
      else
        echo "Packages already installed. Skipping."
      fi
      
      echo activating the environment ...
      # Activate the environment
      eval "$(micromamba shell hook --shell=posix)"
      micromamba activate ${name}

      exec jupyter lab
      
      set +e
    '';
  };
in fhs.env


 
