#!/bin/bash
#SBATCH --job-name=17report
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/17report.out
#SBATCH --error=logs/17report.err

# step 17: CNV genomic profile report (rolling median per chr × race × type)

source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline
quarto render "$(dirname "$0")/17-cnv-qc-report.qmd"
