{ pkgs ? import <nixpkgs> {}}:
let
  PGDATA = "${toString ./.}/.pg";
  PORT = "5432";
  CRSRA = "SRR14243987";
  postgresConf = pkgs.writeText "postgresql.conf"
      ''
        # Add Custom Settings
        log_min_messages = warning
        log_min_error_statement = error
        log_min_duration_statement = 100  # ms
        log_connections = on
        log_disconnections = on
        log_duration = on
        #log_line_prefix = '[] '
        log_timezone = 'UTC'
        log_statement = 'all'
        log_directory = 'pg_log'
        log_filename = 'postgresql-%Y-%m-%d_%H%M%S.log'
        logging_collector = on
        log_min_error_statement = error
      '';
   fhs =  pkgs.buildFHSEnv rec{
     name = "immune-nets";

    targetPkgs = _: [
      pkgs.micromamba
	    pkgs.postgresql_15
    ];

    profile = ''
      set -e
      export MAMBA_ROOT_PREFIX=${builtins.getEnv "PWD"}/.mamba
      # Check if the environment already exists
      if micromamba env list | grep -q "^${name}"; then
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

      echo "Using ${pkgs.postgresql_15.name}."

      # Setup: other env variables
      export SRA=${CRSRA}
      export PGPORT="${PORT}"
      export PGHOST="${PGDATA}"
      export PGDATA="${PGDATA}"
      [ ! -d ${PGDATA} ] && pg_ctl initdb -o "-U postgres --no-locale" && cat "${postgresConf}" >> ${PGDATA}/postgresql.conf
      pg_ctl -o "-p $PGPORT -k ${PGDATA}" start
      micromamba run scripts/envrc
      echo all done


      set +e
    '';
  };
in fhs.env


 
