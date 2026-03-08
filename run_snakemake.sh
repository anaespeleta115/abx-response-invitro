#!/bin/bash
#SBATCH --job-name=midas_pipeline
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=snakemake_%j.out
#SBATCH --error=snakemake_%j.err

source /dfs7/xuek5-lab/aespelet/miniforge3/etc/profile.d/conda.sh
conda activate snakemake

cd /dfs7/xuek5-lab/aespelet/abx-response-invitro

snakemake --cores 16
