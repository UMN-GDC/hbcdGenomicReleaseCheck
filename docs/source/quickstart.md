# Quick Start Guide

This guide walks you through running the complete Phase C pipeline (de-identification, filtering, and validation) for an HBCD genomic data release.

## Prerequisites

- Access to HBCD shared data directories
- MSI/Slurm cluster access (for SLURM jobs)
- `gdcPipeline` conda environment
- Required modules loaded (plink2, bcftools, R, gcta, king)

## 1. Environment Setup

```bash
# Load required modules (MSI)
module load plink/2.00-alpha-091019
module load bcftools
module load R/4.4.2-openblas-rocky8
module load gcta/1.94.1
module load king/2.2.5

# Set release tag (br31p2 for current target)
export HBCD_RELEASE=br31p2

# Set data directories (adjust for your environment)
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData

# Optional: Override crosswalk/par_visit files for this release
# (run_release_pipeline.sh exports these automatically)
# export HBCD_IDENTIFIERS_FILE=$HBCD_DATA_DIR/release_identifiers_20261201.csv
# export HBCD_PAR_VISIT_FILE=$HBCD_DATA_DIR/par_visit_data_br31_2.tsv
```

> **Note**: The wrapper script `run_release_pipeline.sh` uses `conda run -n gdcPipeline` (full path to env) to avoid activate/deactivate bugs. No manual `conda activate` needed.

## 2. Verify Phase A & B Complete

Before running Phase C, ensure Phase A (QC & derivatives) and Phase B (de-identification) are complete:

```bash
# Check Phase A outputs exist
ls -la $HBCD_DATA_DIR/onlyQc.{bed,bim,fam}
ls -la $HBCD_DATA_DIR/data_handoff/hbcd_pcrelate_grm.grm.*
ls -la $HBCD_DATA_DIR/data_handoff/hbcd_pcair_32PCs_clean.tsv

# Check Phase B de-identification complete
# (onlyQc.fam should have release_candid-based IIDs like 1234567890C)
head -5 $HBCD_DATA_DIR/onlyQc.fam
```

## 3. Run Phase C Pipeline

### Option A: Automated Wrapper (Recommended)

```bash
# Runs all Phase C steps in order with SLURM job monitoring
./run_release_pipeline.sh
```

The wrapper script:
1. Runs step 25 (de-ID + filter genotypes)
2. Runs step 26 (PLINK2 --keep)
3. Submits step 27 (SLURM array for VCF filtering) and waits
4. Runs step 28 (CNV de-identification)
5. Runs step 29 (filter all derivatives)
6. Runs step 30 (validate exclusions)

### Option B: Manual Step-by-Step

```bash
# Step 25: Build keep list & temp.fam (de-ID + filter)
python 25-filterGenotypeFiles.py

# Step 26: PLINK2 --keep → release PLINK files
./26-run_plink_filter.sh

# Step 27: Filter imputed VCFs (SLURM array, 24 tasks)
sbatch 27-filter_imputed_vcf.SLURM
# Wait for completion:
squeue -u $USER --name=filter_imputed

# Step 28: CNV de-identification (pscid → release_candid)
python 28-cnv-deid.py

# Step 29: Filter all derivatives to release subjects
python 29-filter_release_outputs.py

# Step 30: Validate no excluded subjects in outputs
python 30-validateExclusions.py

# Run full test suite
python -m pytest tests/ -v
```

## 4. Verify Release Outputs

```bash
# Check release directory structure
RELEASE_DIR=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}/genotype_microarray
tree $RELEASE_DIR

# Expected structure:
# genotype_microarray/
# ├── GDA/
# │   ├── merged_chroms.{bed,bim,fam}
# │   ├── batch.info
# │   └── removed_individuals.txt
# ├── genesis/
# │   ├── pcair_weights.tsv
# │   ├── pcrelate_relatedness.grm.{id,bin,N.bin}
# │   ├── pcrelate_relatedness.grm.gz
# │   └── pcrelate_relatedness.tsv
# ├── cnv/
# │   ├── CNV_slim.txt
# │   └── CNV_bookmarks.csv
# ├── imputed/
# │   ├── chr*_dose.vcf.gz + .tbi
# │   ├── chr*.info.gz
# │   ├── batch*-snps-typed-only.txt
# │   ├── batch*-snps-excluded.txt
# │   ├── batch*-chunks-excluded.txt
# │   └── batch*-quality-control.html
# ├── keep_list.txt
# ├── temp.fam
# └── staging/
#     └── cnv/
#         ├── CNV_slim_deid.txt
#         └── CNV_bookmarks_deid.csv
```

## 5. Run Validation Tests

```bash
# Full test suite
python -m pytest tests/ -v

# Specific test categories
python -m pytest tests/test_release_data.py -v      # De-ID & data integrity
python -m pytest tests/test_exclusions.py -v        # Exclusion compliance
```

## 6. Check Logs

Each step produces detailed logs:

```bash
# Phase C logs in release base directory
RELEASE_BASE=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}
ls -la $RELEASE_BASE/log_*

# Key log files:
# log_25_filterGenotypeFiles_*.txt
# log_26_run_plink_filter_*.txt
# log_29_filter_release_outputs_*.txt
# log_30_validateExclusions_*.txt
```

## Expected Runtime (MSI)

| Step | Type | Est. Time | Resources |
|------|------|-----------|-----------|
| 25 | Python | ~5 min | 1 CPU, 8 GB |
| 26 | PLINK2 | ~10 min | 1 CPU, 16 GB |
| 27 | SLURM array (24) | ~30 min | 1 CPU, 8 GB each |
| 28 | Python | ~2 min | 1 CPU, 4 GB |
| 29 | Python | ~5 min | 1 CPU, 8 GB |
| 30 | Python | ~3 min | 1 CPU, 4 GB |
| Tests | pytest | ~2 min | 1 CPU, 4 GB |

**Total**: ~1 hour (mostly waiting for SLURM queue)

## Common Issues

| Problem | Solution |
|---------|----------|
| SLURM jobs stuck in queue | Check partition/account: `sbatch --partition=... --account=...` |
| `keep_list.txt` missing | Re-run step 25; check `par_visit` and exclusion files |
| VCF filtering fails | Verify `$HBCD_IMPUTATION_DIR/imputed/c{1,2,3}/chr*.dose.vcf.gz` exist |
| CNV de-ID fails | Check `$HBCD_DATA_DIR/data_handoff/CNV_slim_clean.txt` exists |
| Test failures | Run `python 30-validateExclusions.py` first for detailed output |

## Next Steps

- Review [Outputs](outputs.md) for detailed file descriptions
- See [Configuration](configuration.md) for environment variables
- Check [Troubleshooting](troubleshooting.md) for common problems