#!/bin/bash
#SBATCH --job-name=midas_pipeline
#SBATCH --account=xuek5
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=snakemake_%j.out
#SBATCH --error=snakemake_%j.err

source /dfs7/xuek5-lab/aespelet/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

cd /dfs7/xuek5-lab/aespelet/myrepos/abx-response-invitro

snakemake --cores 4 --rerun-incomplete
