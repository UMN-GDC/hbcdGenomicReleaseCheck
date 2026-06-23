#!/bin/bash
#SBATCH --job-name=33report
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/33report.out
#SBATCH --error=logs/33report.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline
quarto render 33-cnv-qc-report.qmd
