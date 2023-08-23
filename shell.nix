with import <nixpkgs> {};

stdenv.mkDerivation {
  name = "postgres-env";
  buildInputs = [];

  nativeBuildInputs = with pkgs.buildPackages; [
    python38
    sratoolkit
    vim
    postgresql_15
  ];

  postgresConf =
    writeText "postgresql.conf"
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

  PGDATA = "${toString ./.}/.pg";
  PORT = "5432";
  CRSRA = "SRR14243987";
  shellHook = ''
    echo "Using ${postgresql_15.name}."

    # Setup: other env variables
    export SRA=$CRSRA
    export PGPORT="$PORT"
    export PGHOST="$PGDATA"
    [ ! -d $PGDATA ] && pg_ctl initdb -o "-U postgres" && cat "$postgresConf" >> $PGDATA/postgresql.conf
    pg_ctl -o "-p $PGPORT -k $PGDATA" start
    alias fin="pg_ctl stop && exit"
    alias pg="psql -p $PGPORT -U postgres"
    alias dcrds="sh ./scripts/clonotype-processor.sh -s $SRA -n $(date -I) -d ../datasets"
  '';
}