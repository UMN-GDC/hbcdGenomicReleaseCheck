#!/bin/bash
#SBATCH --job-name=32report
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/32report.out
#SBATCH --error=logs/32report.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline
quarto render 32-cnv-qc-report.qmd
