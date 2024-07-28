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
  fhs = pkgs.buildFHSUserEnv {
    name = "my-fhs-environment";

    targetPkgs = _: [
      pkgs.micromamba
	  pkgs.postgresql_15
    ];

    profile = ''
      set -e
      eval "$(micromamba shell hook --shell=posix)"
      export MAMBA_ROOT_PREFIX=${builtins.getEnv "PWD"}/.mamba
	 	  
      micromamba create -q -n immune-nets
      micromamba activate immune-nets
	  micromamba install --yes -f env.yml -c conda-forge

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


 
