#!/bin/bash
#SBATCH --job-name=31report
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/31report.out
#SBATCH --error=logs/31report.err

# Get the download link hash from the TopMed Server job results summary

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline
quarto render 31-cnv-qc-report.qmd
