#!/bin/bash
#SBATCH -p obl
#SBATCH -w obl1
#SBATCH --ntasks=9
#SBATCH -A inz_tcr
#SBATCH --gres=gpu:0
#SBATCH --time=05-01:05:10

ulimit -l 4194304 || true
eval "$(/home/inf155621/miniconda3/bin/conda shell.bash hook)"
conda activate immune-nets
mpirun -np 9 python -m src.cluster.optimizerClusterJob --group-paths "covid:/home/inf155621/immune-nets-datasets/covid,healthy:/home/inf155621/immune-nets-datasets/healthy" --run-id 1
