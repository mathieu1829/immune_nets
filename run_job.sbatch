#!/bin/bash
#SBATCH --job-name=immune_nets_job
#SBATCH --output=logs/immune_nets_%j.out
#SBATCH --error=logs/immune_nets_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1

cd /home/teuthid/Documents/python/immune_nets

NIX_SHELL_COMMAND="python -m src.analysis.clusterJob" nix-shell
