#!/bin/bash
#SBATCH -p obl
#SBATCH -w obl1
#SBATCH -n1
#SBATCH -A inz_tcr
#SBATCH --gres=gpu:0
#SBATCH --time=05-01:05:10

/home/inf155621/miniconda3/bin/conda run -n immune-nets python -m tests.test_clusterjob

