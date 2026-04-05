#!/bin/bash
#SBATCH -p obl
#SBATCH -w obl1
#SBATCH --ntasks=9
#SBATCH -A inz_tcr
#SBATCH --gres=gpu:0
#SBATCH --time=02-00:00:00

ulimit -l 4194304 || true
eval "$(/home/inf155621/miniconda3/bin/conda shell.bash hook)"
conda activate immune-nets
mpirun -np 9 python -m immune_nets.cluster.optimizerClusterJob --group-paths "covid:/home/inf155621/immune-nets-datasets/covid,healthy:/home/inf155621/immune-nets-datasets/healthy" --run-id 1
