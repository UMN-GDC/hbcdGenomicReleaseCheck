#!/bin/bash -l
#SBATCH --job-name=installpcaAIR
#SBATCH --output=logs/pcair.out
#SBATCH --error=logs/pcair.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH -p agsmall

# step 6: compute GCTA GRM from harmonized PLINK files
# Output written to DATA_DIR.
# The GRM is computed in chromosome chunks (--make-grm-part) then merged.

module load R/4.4.2-openblas-rocky8
Rscript install-bioc-packages.R
