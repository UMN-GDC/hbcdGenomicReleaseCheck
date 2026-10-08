# Running the Pipeline

## Overview

The pipeline has three phases. Phase A and B are typically run once to produce the `onlyQc` dataset and de-identified derivatives. Phase C is run for each release to filter to the release subject subset.

```mermaid
graph LR
    A[Phase A: QC & Derivatives] --> B[Phase B: De-ID]
    B --> C[Phase C: Release Filter]
    C --> D[Release Outputs]
```

## Phase A: QC & Derivative Computation (Scripts 01-17)

Run these **once** to produce the `onlyQc` dataset and all downstream derivatives. These scripts operate on raw `pscid` IDs and produce **non-de-identified** outputs.

### Prerequisites

- Raw HBCD genotype data in `$HBCD_DATA_DIR/HBCD.{bed,bim,fam}`
- QC exclusion Excel file
- par_visit participation data
- MSI/Slurm access for GRM, PC-AiR, PC-Relate, imputation steps

### Step-by-Step

```bash
# 0. Environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate gdcPipeline
module load plink/2.00-alpha-091019
module load bcftools
module load R/4.4.2-openblas-rocky8
module load gcta/1.94.1
module load king/2.2.5

# 1. QC Removal → onlyQc
python 01-onlyQCremoved.py          # Maps HBCD.fam IIDs → release_candid; produces QC_removed.txt
bash 02-plink_qc_remove.sh          # PLINK2 --remove QC + controls → onlyQc.{bed,bim,fam}

# 2. GRM Computation
sbatch 03-grm_gcta.SLURM            # GCTA GRM (per-chr parts, merged) → hbcd_gcta_grm.*
sbatch 04-grm_plink.SLURM           # PLINK2 GRM (--make-rel) → hbcd_plink_grm.rel.*

# 3. PC-AiR
sbatch 05-pcair.SLURM               # Runs 06-pca_ir_pipeline.R (KING → PC-AiR)
# Outputs: hbcd_rsid_harmonized_pc_scores.txt, PC-AiR RDS, unrelated/related IDs

# 4. PC-Relate
sbatch 07-pcrelate.sh               # Runs 08-pca_relate_pipeline.R
Rscript 09-export_pcrelate.R        # Export kinship/IBD to CSV + TSV

# 5. Imputation: Prep → TOPMed → Download → Unzip
sbatch --array=0-22 10a-prepare_imputation.SLURM    # Autosomes → per-chr VCF
sbatch 10b-prepare_imputation_X.SLURM               # ChrX VCF
bash 11-submit_topmed.sh            # Upload to TOPMed (needs ~/topmedKey)
sbatch 12a-download_chunk1.SLURM    # Download chunk 1
sbatch 12b-download_chunk2.SLURM    # Download chunk 2
sbatch 12c-download_chunk3.SLURM    # Download chunk 3
sbatch 13a-unzip.SLURM              # Unzip chunk 1 (array 1-7)
sbatch 13b-unzip.SLURM              # Unzip chunk 2
sbatch 13c-unzip.SLURM              # Unzip chunk 3
sbatch 13d-unzip.SLURM              # Unzip chrX

# 6. Imputation QC & Reports
bash futureScripts/14-submit_gp_prep.sh              # Submits GP precompute array
sbatch futureScripts/15-aggregate_gp.SLURM           # Aggregate GP results
quarto render futureScripts/16-imputation-quality-report.qmd
sbatch futureScripts/17-cnv-qc-report.SLURM          # CNV genomic profile report
```

### Phase A Outputs (All use `pscid` IDs)

```
$HBCD_DATA_DIR/
├── onlyQc.{bed,bim,fam}                    # QC-passing genotypes
├── hbcd_gcta_grm.*                         # GCTA GRM
├── hbcd_plink_grm.rel.*                    # PLINK GRM
├── hbcd_rsid_harmonized_pc_scores.txt      # PC-AiR scores
├── pcair_*.rds                             # PC-AiR RDS object
├── pcrelate_*                              # PC-Relate outputs
├── data_handoff/
│   ├── hbcd_pcrelate_grm.grm.*             # For Phase C
│   ├── hbcd_pcair_32PCs_clean.tsv
│   ├── hbcd_pcrelate_grm_pairwise.tsv
│   ├── CNV_slim_clean.txt                  # pscid IIDs
│   └── HBCD_CNV_bookmark_metrics_clean.csv # Already de-IDed
└── imputed/                                # TOPMed results
    ├── c1/chr*.dose.vcf.gz
    ├── c2/chr*.dose.vcf.gz
    ├── c3/chr*.dose.vcf.gz
    └── chrX.dose.vcf.gz
```

## Phase B: De-identification (Intermediate, Run Once)

After Phase A, map all `pscid` IDs to anonymous `release_candid` integers using the identifiers crosswalk.

```bash
# This is a batch rename step - typically run once after Phase A completes
# Maps pscid → release_candid in all Phase A outputs
# Produces de-identified versions of every file in place

# Example (custom script - not in repo):
python deidentify_phase_a_outputs.py \
  --identifiers $HBCD_DATA_DIR/release_identifiers_20260628.csv \
  --input-dir $HBCD_DATA_DIR \
  --output-dir $HBCD_DATA_DIR_deid
```

**After Phase B:**
- Every PLINK IID becomes `{release_candid}{C|M}` (e.g., `1234567890C`)
- Every derivative file (GRM, PC-AiR, PC-Relate) uses de-identified IDs
- Original `pscid` IDs no longer present in any `data/` file

## Phase C: Release Filtering & Validation (Scripts 25-30)

Run these **for each release** to filter de-identified data to the release subject subset.

### Automated (Recommended)

```bash
export HBCD_RELEASE=br31p2
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
./run_release_pipeline.sh
```

### Manual Steps

#### Step 25: Genotype De-ID + Filter (`25-filterGenotypeFiles.py`)

```bash
python 25-filterGenotypeFiles.py
```

**What it does:**
1. Loads identifiers, batch info, par_visit, exclusion lists
2. Two-stage filter:
   - Stage 1: QC + control removal → `onlyQc` set
   - Stage 2: Exclusion list removal (HBCDexclusions.csv + Excel)
3. Merges .fam with identifiers (pscid → release_candid) and batch info
4. Relationship filter: keeps only matching C/M suffix
5. Outputs:
   - `temp.fam` — ALL subjects with de-ID FID/IID (preserves .bed row count)
   - `keep_list.txt` — Release whitelist (valid subjects with metadata)
   - `GDA/batch.info` — IID, visit, plate_number for release subjects
   - `GDA/removed_individuals.txt` — Excluded IIDs for documentation

#### Step 26: PLINK2 Filter (`26-run_plink_filter.sh`)

```bash
./26-run_plink_filter.sh
```

**What it does:**
```bash
plink2 --bfile onlyQc \
       --allow-extra-chr \
       --fam temp.fam \
       --keep keep_list.txt \
       --make-bed --out GDA/merged_chroms
```
- `--bfile` = de-identified `onlyQc` data
- `--fam temp.fam` = supplies de-ID FID/IID + original row count
- `--keep` = restricts to release whitelist
- Output: `GDA/merged_chroms.{bed,bim,fam}`

#### Step 27: Filter Imputed VCFs (`27-filter_imputed_vcf.SLURM`)

```bash
sbatch 27-filter_imputed_vcf.SLURM
# Wait for completion:
watch "squeue -u $USER --name=filter_imputed"
```

**What it does (SLURM array, 24 tasks = chr 1-22 + X):**
1. Finds imputed VCF in `$HBCD_IMPUTATION_DIR/imputed/c{1,2,3}/`
2. VCF samples already de-identified (`{release_candid}{C|M}`)
3. Intersects VCF samples with `keep_list.txt`
4. `bcftools view --samples-file` → filtered `chr{chr}_dose.vcf.gz`
5. Copies `.info.gz` + QC statistics files (batch*-prefixed)

#### Step 28: CNV De-identification (`28-cnv-deid.py`)

```bash
python 28-cnv-deid.py
```

**What it does:**
- **CNV slim clean**: Reads `data_handoff/CNV_slim_clean.txt` (pscid IIDs like `GSM0000000_Grn_12345C`)
  - Extracts pscid + relationship (C/M)
  - Maps pscid → release_candid via crosswalk
  - Writes `staging/cnv/CNV_slim_deid.txt` with de-ID IIDs
- **CNV bookmarks**: Copies `data_handoff/HBCD_CNV_bookmark_metrics_clean.csv` (already de-IDed) to `staging/cnv/CNV_bookmarks_deid.csv`

#### Step 29: Filter All Derivatives (`29-filter_release_outputs.py`)

```bash
python 29-filter_release_outputs.py
```

**What it does:** Reads each de-IDed derivative from `data_handoff/`, keeps only rows with IIDs in `keep_list.txt`, writes to release directory.

| Output | Source | Filter Column | Method |
|--------|--------|---------------|--------|
| `genesis/pcrelate_relatedness.grm.*` | `hbcd_pcrelate_grm.grm.*` | IID | Binary GRM subset |
| `genesis/pcrelate_relatedness.grm.gz` | `hbcd_pcrelate_grm.grm.gz` | IID1, IID2 | Text GRM index mapping |
| `genesis/pcrelate_relatedness.tsv` | `hbcd_pcrelate_grm_pairwise.tsv` | ID1, ID2 | CSV filter |
| `genesis/pcair_weights.tsv` | `hbcd_pcair_32PCs_clean.tsv` | subject_id | TSV filter |
| `cnv/CNV_slim.txt` | `staging/cnv/CNV_slim_deid.txt` | sample_id | TSV filter |
| `cnv/CNV_bookmarks.csv` | `staging/cnv/CNV_bookmarks_deid.csv` | sample_id | CSV filter |

#### Step 30: Validate Exclusions (`30-validateExclusions.py`)

```bash
python 30-validateExclusions.py
```

**What it checks:**
- `GDA/merged_chroms.fam` — FID and IID level
- `HBCDexclusions.csv` — Per-column, mapped to release_candid
- All derivative files — IID columns vs excluded-IID set
- `cnv/CNV_slim.txt` — sample_id column
- `cnv/CNV_bookmarks.csv` — sample_id column
- Cross-check — warns of unexpected IIDs (neither release nor excluded)

### Run Tests

```bash
python -m pytest tests/ -v
```

## Re-running After Exclusion List Changes

If `HBCDexclusions.csv`, Excel QC sheet, or other exclusion sources change:

```bash
# Phase A & B do NOT need re-run (derivatives still valid)

# Re-run Phase C only:
python 25-filterGenotypeFiles.py      # Rebuild keep list
./26-run_plink_filter.sh              # Re-filter PLINK
sbatch 27-filter_imputed_vcf.SLURM    # Re-filter VCFs
python 28-cnv-deid.py                 # Re-de-ID CNV (source always pscid)
python 29-filter_release_outputs.py   # Re-filter all derivatives
python 30-validateExclusions.py       # Re-validate
python -m pytest tests/ -v            # Re-test
```

## Monitoring SLURM Jobs

```bash
# Check job status
squeue -u $USER

# Check specific job
squeue -j <JOB_ID>

# View job output
cat slurm-<JOB_ID>.out

# Cancel job
scancel <JOB_ID>

# For array jobs (step 27)
squeue -u $USER --name=filter_imputed
```

## Log Files

Each step produces logs in the release base directory:

```
$RELEASE_BASE/
├── log_25_filterGenotypeFiles_YYYYMMDD_HHMMSS.txt
├── log_26_run_plink_filter_YYYYMMDD_HHMMSS.txt
├── log_29_filter_release_outputs_YYYYMMDD_HHMMSS.txt
├── log_30_validateExclusions_YYYYMMDD_HHMMSS.txt
└── slurm-<JOB_ID>.out               # SLURM job outputs
```

## Parallel Execution Notes

- Steps 27, 28, 29 are **independent** after step 26 completes
- Can run 27 (SLURM), 28, 29 in parallel
- Step 30 and tests must run last (depend on all outputs)