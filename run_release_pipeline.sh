#!/bin/bash
#SBATCH --job-name=hbcd_release_pipeline
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=msismall
#SBATCH --output=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_%x_%j.out
#SBATCH --error=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_%x_%j.err

# run_release_pipeline.sh — Create the full HBCD release directory
# Runs all pipeline steps in order: 25 → 26 → 27 → 28 → 29 → 30.
#
# Usage (as SLURM batch job):
#   export HBCD_RELEASE=br_21p3
#   sbatch run_release_pipeline.sh
#
# Or with explicit release tag:
#   HBCD_RELEASE=br_21p3 sbatch run_release_pipeline.sh
#
# Env vars:
#   HBCD_RELEASE         — release tag (default: br_21p3)
#   HBCD_DATA_DIR        — source data dir (default: shared/data)
#   HBCD_IMPUTATION_DIR  — imputation data root (default: shared/hbcdSandboxData)
#   RELEASE_DIR          — override release directory path

set -euo pipefail

# ── Initialize environment for SLURM batch ───────────────────────────────
# Source profile for modules and conda
source /etc/profile.d/modules.sh
source ~/miniconda3/etc/profile.d/conda.sh

# Activate environments
source /projects/standard/basu_hbcd/shared/.venv/bin/activate
conda activate /projects/standard/gdc/public/envs/gdcPipeline

# ── defaults ────────────────────────────────────────────────────
: "${HBCD_RELEASE:=br_21p3}"
: "${HBCD_DATA_DIR:=/projects/standard/basu_hbcd/shared/data}"
: "${HBCD_IMPUTATION_DIR:=/projects/standard/basu_hbcd/shared/hbcdSandboxData}"
: "${RELEASE_DIR:=$(dirname "$HBCD_DATA_DIR")/HBCD_genomics_release_${HBCD_RELEASE}/genotype_microarray}"

RELEASE_BASE="$(dirname "$RELEASE_DIR")"

# Exported paths for downstream scripts / manual steps
export HBCD_IDENTIFIERS_FILE="${HBCD_DATA_DIR}/release_identifiers_20260628.csv"
export HBCD_PAR_VISIT_FILE="${HBCD_DATA_DIR}/par_visit_data_br21_1.tsv"

echo "=========================================="
echo "  HBCD Release Pipeline"
echo "  Release tag : $HBCD_RELEASE"
echo "  Release dir : $RELEASE_DIR"
echo "  Data dir    : $HBCD_DATA_DIR"
echo "  Imputation  : $HBCD_IMPUTATION_DIR"
echo "  SLURM Job   : ${SLURM_JOB_ID:-local}"
echo "=========================================="

# ── Step 25: De-ID + Filter genotypes ───────────────────────────
echo ""
echo "=== Step 25: De-ID + Filter genotypes ==="
python 25-filterGenotypeFiles.py

# ── Step 26: PLINK2 --keep → GDA/merged_chroms ─────────────────
echo ""
echo "=== Step 26: PLINK2 --keep ==="
module load plink/2.00-alpha-091019
bash 26-run_plink_filter.sh
module unload plink/2.00-alpha-091019 2>/dev/null || true

# ── Step 27: Filter imputed VCFs to release subjects ────────────
# This is a SLURM array job (24 tasks). Submit and poll until complete.
echo ""
echo "=== Step 27: Filter imputed VCFs (SLURM array) ==="
JOB_ID=$(sbatch --parsable 27-filter_imputed_vcf.SLURM)
echo "  Submitted job: $JOB_ID"
while squeue -j "$JOB_ID" --noheader 2>/dev/null | grep -q .; do
    echo "  Waiting for job $JOB_ID..."
    sleep 30
done
echo "  Job $JOB_ID complete."

# ── Step 28: De-identify CNV calls ──────────────────────────────
echo ""
echo "=== Step 28: De-identify CNV ==="
python 28-cnv-deid.py

# ── Step 29: Filter all handoff derivatives to release subjects ─
echo ""
echo "=== Step 29: Filter handoff derivatives ==="
python 29-filter_release_outputs.py

# ── Step 30: Validate exclusion integrity ───────────────────────
echo ""
echo "=== Step 30: Validate exclusions ==="
python 30-validateExclusions.py

echo ""
echo "=========================================="
echo "  Release pipeline complete!"
echo "  Output: $RELEASE_DIR"
echo "=========================================="
echo "  Run tests:  python -m pytest tests/ -v"
echo "=========================================="
