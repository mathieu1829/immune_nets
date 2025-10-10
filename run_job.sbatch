#!/bin/bash
#SBATCH --job-name=immune_nets_job
#SBATCH --output=logs/immune_nets_%j.out
#SBATCH --error=logs/immune_nets_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

cd /home/teuthid/Documents/python/immune-nets

nix-shell --run "python -m src.analysis.clusterJob"
