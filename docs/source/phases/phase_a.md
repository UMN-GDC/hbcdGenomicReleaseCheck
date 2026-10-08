# Phase A: QC & Derivative Computation

## Overview

Phase A runs **once** to transform raw HBCD genotype data into the `onlyQc` dataset and compute all downstream derivatives (GRM, PC-AiR, PC-Relate, imputation). All outputs use original `pscid` IDs and are **not de-identified**.

## Inputs

| File | Location | Description |
|------|----------|-------------|
| `HBCD.{bed,bim,fam}` | `$HBCD_DATA_DIR` | Raw genotype data |
| `HBCD.fam` | `$HBCD_DATA_DIR` | Original PLINK .fam with pscid IIDs |
| `batch.info` | `$HBCD_DATA_DIR` | Batch metadata (IID={pscid}{C\|M}) |
| QC Excel | `$HBCD_DATA_DIR` | `HBCD_genetics_QC1_missing_race_LORIS.xlsx` |
| par_visit | `$HBCD_DATA_DIR` | `par_visit_data_br21_1.tsv` |

## Steps

### Step 01: `01-onlyQCremoved.py`

**Purpose**: Map HBCD.fam IIDs → `release_candid`; produce `QC_removed.txt`

```python
# Key operations:
# 1. Load HBCD.fam (pscid IIDs like GSM0000000_Grn_12345C)
# 2. Load identifiers (pscid → release_candid)
# 3. Load QC exclusion Excel (Study_ID column)
# 4. Identify QC failures + control subjects
# 5. Output: QC_removed.txt (pscids to remove)
# 6. Remaining subjects = onlyQc set
```

**Outputs**:
- `QC_removed.txt` — PSCIDs removed by QC

### Step 02: `02-plink_qc_remove.sh`

**Purpose**: PLINK2 `--remove` to produce `onlyQc.{bed,bim,fam}`

```bash
plink2 --bfile HBCD \
       --remove QC_removed.txt \
       --make-bed --out onlyQc
```

**Outputs**:
- `onlyQc.{bed,bim,fam}` — QC-passing genotypes (pscid IIDs)

### Step 03: `03-grm_gcta.SLURM`

**Purpose**: GCTA GRM computation (per-chromosome parts, merged)

```bash
# SLURM array 1-22 (autosomes)
# Per chromosome:
gcta64 --bfile onlyQc --chr {chr} --make-grm-part {chr} --out hbcd_gcta_grm_chr{chr}
# Merge:
gcta64 --mgrm hbcd_gcta_grm_chr*.grm.id --make-grm --out hbcd_gcta_grm
```

**Outputs**:
- `hbcd_gcta_grm.grm.{id,bin,N.bin}` — Merged GCTA GRM

### Step 04: `04-grm_plink.SLURM`

**Purpose**: PLINK2 GRM (`--make-rel`)

```bash
plink2 --bfile onlyQc \
       --make-rel square \
       --out hbcd_plink_grm_rel
```

**Outputs**:
- `hbcd_plink_grm_rel.rel.{id,bin}` — PLINK GRM

### Step 05: `05-pcair.SLURM`

**Purpose**: PC-AiR via KING kinship → ancestry-adjusted PCs

```bash
# Runs 06-pca_ir_pipeline.R
# 1. KING kinship estimation on onlyQc
# 2. PC-AiR on unrelated subset
# 3. Project PCs to all subjects
```

**Script**: `06-pca_ir_pipeline.R`

**Outputs**:
- `hbcd_rsid_harmonized_pc_scores.txt` — PC scores (all subjects)
- `pcair_*.rds` — PC-AiR RDS object
- `unrelated_ids.txt`, `related_ids.txt` — Subject partitions

### Step 06: `06-pca_ir_pipeline.R`

**Purpose**: R implementation of PC-AiR pipeline

**Dependencies**: `GENESIS`, `SNPRelate`, `GWASTools`, `GDSfmt`, `tidyverse`

### Step 07: `07-pcrelate.sh`

**Purpose**: Wrapper for PC-Relate pipeline

```bash
# Runs 08-pca_relate_pipeline.R
```

### Step 08: `08-pca_relate_pipeline.R`

**Purpose**: PC-Relate pairwise relatedness estimation using GENESIS

```r
# Key steps:
# 1. Load PC-AiR results
# 2. PC-Relate on all subject pairs
# 3. Output: IBD probabilities, kinship matrices
```

**Outputs**:
- `pcrelate_pairs` — Pairwise IBD
- `pcrelate_ibd` — IBD estimates
- `pcrelate_self` — Self-relatedness
- `kinmat_wide` — Kinship matrix (wide)
- `kinmat_long` — Kinship matrix (long)

### Step 09: `09-export_pcrelate.R`

**Purpose**: Export kinship/IBD to CSV + TSV for handoff

**Outputs** (to `data_handoff/`):
- `hbcd_pcrelate_grm.grm.{id,bin,N.bin}` — Binary GRM
- `hbcd_pcrelate_grm.grm.gz` — Text GRM (gzipped)
- `hbcd_pcrelate_grm_pairwise.tsv` — Pairwise relatedness
- `hbcd_pcair_32PCs_clean.tsv` — PC-AiR weights (32 PCs)

### Steps 10-13: Imputation Pipeline

#### Step 10a: `10a-prepare_imputation.SLURM` (array 0-22)

**Purpose**: Prepare autosomes 1-22 for TOPMed imputation

```bash
plink2 --bfile onlyQc --chr {chr} \
       --snps-only just-acgt \
       --exclude palindromic_snps.txt \
       --export vcf-4.2 bgz \
       --output-chr chrM \
       --out onlyQc_chr{chr}
```

#### Step 10b: `10b-prepare_imputation_X.SLURM`

**Purpose**: Prepare chrX with special handling

```bash
plink2 --bfile onlyQc --chr X \
       --set-invalid-haploid-missing \
       --export vcf-4.2 bgz \
       --output-chr chrM \
       --out onlyQc_chrX
```

#### Step 11: `11-submit_topmed.sh`

**Purpose**: Upload per-chr VCFs to TOPMed Imputation Server

```bash
# Requires ~/topmedKey (API token)
# Uses TOPMed REST API
# Submits each chromosome as separate job
```

#### Steps 12a-c: Download Chunks

```bash
sbatch 12a-download_chunk1.SLURM  # c1
sbatch 12b-download_chunk2.SLURM  # c2
sbatch 12c-download_chunk3.SLURM  # c3
```

#### Steps 13a-d: Unzip Chunks

```bash
sbatch 13a-unzip.SLURM   # array 1-7
sbatch 13b-unzip.SLURM
sbatch 13c-unzip.SLURM
sbatch 13d-unzip.SLURM   # chrX
```

**Outputs** (to `$HBCD_IMPUTATION_DIR/imputed/`):
- `c1/chr{1-22}.dose.vcf.gz` + `.info.gz`
- `c2/chr{1-22}.dose.vcf.gz` + `.info.gz`
- `c3/chr{1-22}.dose.vcf.gz` + `.info.gz`
- `chrX.dose.vcf.gz` + `.info.gz`

### Steps 14-17: Imputation QC & Reports (Future Scripts)

#### Step 14: `futureScripts/14a-precompute_gp.SLURM`

**Purpose**: GP precomputation per chromosome × ancestry group

```bash
# Ancestry groups: WHT, BLK, AIAN, ASN, HPI, 2PLUS, OTH, UNK
# Low-MAF variants only
bash futureScripts/14-submit_gp_prep.sh  # Submits array
```

#### Step 15: `futureScripts/15-aggregate_gp.SLURM`

**Purpose**: Aggregate per-chr GP results → final CSVs

#### Step 16: `futureScripts/16-imputation-quality-report.qmd`

**Purpose**: Quarto report with imputation quality visualizations

#### Step 17: `futureScripts/17-cnv-qc-report.qmd` / `17-cnv-qc-report.SLURM`

**Purpose**: CNV genomic profile report (rolling median LRR, BAF, CNV per chr)

## Phase A Outputs Summary

```
$HBCD_DATA_DIR/
├── onlyQc.{bed,bim,fam}                           # Main QC-passing genotypes
├── QC_removed.txt                                 # QC exclusions
├── hbcd_gcta_grm.grm.{id,bin,N.bin}              # GCTA GRM
├── hbcd_plink_grm_rel.rel.{id,bin}               # PLINK GRM
├── hbcd_rsid_harmonized_pc_scores.txt            # PC-AiR scores
├── pcair_*.rds                                    # PC-AiR RDS
├── unrelated_ids.txt / related_ids.txt           # Subject partitions
├── pcrelate_*                                     # PC-Relate outputs
├── data_handoff/
│   ├── hbcd_pcrelate_grm.grm.{id,bin,N.bin}     # For Phase C
│   ├── hbcd_pcrelate_grm.grm.gz
│   ├── hbcd_pcrelate_grm_pairwise.tsv
│   ├── hbcd_pcair_32PCs_clean.tsv
│   ├── CNV_slim_clean.txt                        # pscid IIDs
│   └── HBCD_CNV_bookmark_metrics_clean.csv       # Already de-IDed
└── $HBCD_IMPUTATION_DIR/imputed/
    ├── c1/chr*.dose.vcf.gz + .info.gz
    ├── c2/chr*.dose.vcf.gz + .info.gz
    ├── c3/chr*.dose.vcf.gz + .info.gz
    └── chrX.dose.vcf.gz + .info.gz
```

## Notes

- All Phase A outputs use **pscid-based IIDs** (e.g., `GSM0000000_Grn_12345C`)
- Phase B will de-identify these to `release_candid` format
- Imputation steps (10-13) currently reference sandbox paths — adjust `DATA_DIR`/`OUT` in scripts for production
- Steps 14-17 are in `futureScripts/` — not yet integrated into main pipeline