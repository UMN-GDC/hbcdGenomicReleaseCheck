# Phase A Scripts Detailed Reference

## 01-onlyQCremoved.py

**Purpose**: Map HBCD.fam IIDs → `release_candid`; produce `QC_removed.txt`

**Usage**:
```bash
python 01-onlyQCremoved.py
```

**Inputs**:
- `$HBCD_DATA_DIR/HBCD.fam` — Original PLINK .fam (pscid IIDs)
- `$HBCD_DATA_DIR/release_identifiers_20260628.csv` — PSCID → release_candid
- `$HBCD_DATA_DIR/HBCD_genetics_QC1_missing_race_LORIS.xlsx` — QC exclusions (Exclude_Summary sheet, Study_ID column)

**Logic**:
1. Load HBCD.fam, extract pscid from IID (`{array}_{channel}_{pscid}{C|M}`)
2. Load identifiers crosswalk (pscid → release_candid)
3. Load QC exclusion Excel, extract Study_ID → pscid
4. Identify QC failures + control subjects
5. Write `QC_removed.txt` (pscids to remove)
6. Print summary statistics

**Outputs**:
- `$HBCD_DATA_DIR/QC_removed.txt` — One pscid per line

**Key Code Patterns**:
```python
# IID format: GSM0000000_Grn_12345C
pscid = iid.split('_')[-1][:-1]  # Extract pscid
relationship = iid[-1]            # C or M
```

---

## 02-plink_qc_remove.sh

**Purpose**: PLINK2 `--remove` to produce `onlyQc.{bed,bim,fam}`

**Usage**:
```bash
bash 02-plink_qc_remove.sh
```

**Command**:
```bash
plink2 --bfile HBCD \
       --remove QC_removed.txt \
       --make-bed --out onlyQc
```

**Inputs**:
- `HBCD.{bed,bim,fam}` — Raw genotype data
- `QC_removed.txt` — From step 01

**Outputs**:
- `onlyQc.{bed,bim,fam}` — QC-passing genotypes (pscid IIDs)

---

## 03-grm_gcta.SLURM

**Purpose**: GCTA GRM computation (per-chromosome parts, merged)

**Usage**:
```bash
sbatch 03-grm_gcta.SLURM
```

**SLURM Configuration**:
```bash
#SBATCH --array=1-22          # Autosomes
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --partition=msismall
```

**Per-Task Logic** (chromosome `chr`):
```bash
gcta64 --bfile onlyQc --chr {chr} \
       --make-grm-part {chr} \
       --out hbcd_gcta_grm_chr{chr}
```

**Merge Step** (after array completes):
```bash
# Create merge list
ls hbcd_gcta_grm_chr*.grm.id > merge_list.txt
gcta64 --mgrm merge_list.txt --make-grm --out hbcd_gcta_grm
```

**Outputs**:
- `hbcd_gcta_grm_chr{chr}.grm.{id,bin,N.bin}` — Per-chr parts
- `hbcd_gcta_grm.grm.{id,bin,N.bin}` — Merged GRM

---

## 04-grm_plink.SLURM

**Purpose**: PLINK2 GRM (`--make-rel`)

**Usage**:
```bash
sbatch 04-grm_plink.SLURM
```

**SLURM Configuration**:
```bash
#SBATCH --array=1-22
#SBATCH --mem=8G
#SBATCH --time=2:00:00
```

**Command**:
```bash
plink2 --bfile onlyQc --chr {chr} \
       --make-rel square \
       --out hbcd_plink_grm_chr{chr}
```

**Merge** (manual or script):
```bash
# Combine per-chr .rel files
```

**Outputs**:
- `hbcd_plink_grm_chr{chr}.rel.{id,bin}`
- `hbcd_plink_grm_rel.rel.{id,bin}` — Merged

---

## 05-pcair.SLURM

**Purpose**: PC-AiR wrapper — submits R script

**Usage**:
```bash
sbatch 05-pcair.SLURM
```

**SLURM Configuration**:
```bash
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=msilarge
```

**Command**:
```bash
Rscript 06-pca_ir_pipeline.R
```

---

## 06-pca_ir_pipeline.R

**Purpose**: PC-AiR ancestry-adjusted PCs via KING kinship

**Dependencies**: `GENESIS`, `SNPRelate`, `GWASTools`, `GDSfmt`, `tidyverse`, `SessionInfo`

**Key Steps**:
```r
# 1. Convert onlyQc PLINK → GDS
snpgdsBED2GDS("onlyQc.bed", "onlyQc.fam", "onlyQc.bim", "onlyQc.gds")

# 2. KING kinship estimation
king <- snpgdsIBDKING(gdsobj, snp.id=..., num.thread=...)

# 3. Identify unrelated pairs (kinship < 2^(-11/2))
unrelated <- king$kinship < 0.0442

# 4. PC-AiR on unrelated subset
pcair <- pcair(gdsobj, kinobj=king, divobj=..., 
               snp.include=..., num.thread=...)

# 5. Project PCs to all samples
pcs_all <- pcair_project(pcair, gdsobj, ...)

# 6. Save outputs
saveRDS(pcair, "pcair_object.rds")
write.table(pcs_all, "hbcd_rsid_harmonized_pc_scores.txt", ...)
write.table(unrelated_ids, "unrelated_ids.txt", ...)
write.table(related_ids, "related_ids.txt", ...)
```

**Outputs**:
- `hbcd_rsid_harmonized_pc_scores.txt` — PC scores (all subjects)
- `pcair_object.rds` — PC-AiR RDS object
- `unrelated_ids.txt`, `related_ids.txt` — Subject partitions
- `hbcd_pcair_32PCs_clean.tsv` — 32 PCs for handoff (via step 09)

---

## 07-pcrelate.sh

**Purpose**: PC-Relate wrapper

**Usage**:
```bash
bash 07-pcrelate.sh
```

**Command**:
```bash
Rscript 08-pca_relate_pipeline.R
```

---

## 08-pca_relate_pipeline.R

**Purpose**: PC-Relate pairwise relatedness using GENESIS

**Dependencies**: `GENESIS`, `GWASTools`, `GDSfmt`, `SessionInfo`

**Key Steps**:
```r
# 1. Load PC-AiR results
pcair <- readRDS("pcair_object.rds")

# 2. PC-Relate
pcrelate <- pcrelate(gdsobj, pcobj=pcair, 
                     training.set=unrelated_ids,
                     snp.include=..., num.thread=...)

# 3. Extract results
pairs <- pcrelate$pairs
ibd <- pcrelate$ibd
self <- pcrelate$self
kinmat_wide <- pcrelate$kinmat_wide
kinmat_long <- pcrelate$kinmat_long
```

**Outputs** (R objects, exported by step 09):
- `pcrelate_pairs` — Pairwise IBD probabilities
- `pcrelate_ibd` — IBD estimates
- `pcrelate_self` — Self-relatedness
- `kinmat_wide` — Kinship matrix (wide)
- `kinmat_long` — Kinship matrix (long)

---

## 09-export_pcrelate.R

**Purpose**: Export PC-Relate results to handoff directory

**Usage**:
```bash
Rscript 09-export_pcrelate.R
```

**Outputs** (to `$HBCD_DATA_DIR/data_handoff/`):
```r
# Binary GRM (GCTA format)
write_grm_bin(kinmat, "hbcd_pcrelate_grm.grm")

# Text GRM (gzipped)
write_grm_gz(kinmat, "hbcd_pcrelate_grm.grm.gz")

# Pairwise TSV
write.table(pairs, "hbcd_pcrelate_grm_pairwise.tsv", sep="\t", ...)

# PC-AiR 32 PCs clean
write.table(pcair_weights, "hbcd_pcair_32PCs_clean.tsv", sep="\t", ...)
```

**Files produced**:
- `hbcd_pcrelate_grm.grm.{id,bin,N.bin}`
- `hbcd_pcrelate_grm.grm.gz`
- `hbcd_pcrelate_grm_pairwise.tsv`
- `hbcd_pcair_32PCs_clean.tsv`

---

## 10a-prepare_imputation.SLURM

**Purpose**: Prepare autosomes 1-22 for TOPMed imputation

**Usage**:
```bash
sbatch --array=0-22 10a-prepare_imputation.SLURM
```

**SLURM Configuration**:
```bash
#SBATCH --array=0-22
#SBATCH --mem=8G
#SBATCH --time=2:00:00
```

**Per-Task Logic** (task 0-22 → chr 1-22):
```bash
CHR=$((SLURM_ARRAY_TASK_ID + 1))

plink2 --bfile onlyQc --chr ${CHR} \
       --snps-only just-acgt \
       --exclude palindromic_snps.txt \
       --export vcf-4.2 bgz \
       --output-chr chrM \
       --out ${OUT}/onlyQc_chr${CHR}
```

**Hardcoded Paths** (need adjustment for production):
```bash
DATA_DIR=~/hbcdData  # Sandbox path
OUT=${DATA_DIR}/imputation_prep
```

**Outputs** (to `$OUT/`):
- `onlyQc_chr{1-22}.vcf.gz` — Per-chr VCFv4.2

---

## 10b-prepare_imputation_X.SLURM

**Purpose**: Prepare chrX with special handling

**Usage**:
```bash
sbatch 10b-prepare_imputation_X.SLURM
```

**Command**:
```bash
plink2 --bfile onlyQc --chr X \
       --set-invalid-haploid-missing \
       --snps-only just-acgt \
       --exclude palindromic_snps.txt \
       --export vcf-4.2 bgz \
       --output-chr chrM \
       --out ${OUT}/onlyQc_chrX
```

---

## 11-submit_topmed.sh

**Purpose**: Submit per-chr VCFs to TOPMed Imputation Server

**Usage**:
```bash
bash 11-submit_topmed.sh
```

**Requirements**:
- `~/topmedKey` — API token (chmod 600)
- VCFs in `$OUT/` from steps 10a/10b

**Logic**:
```bash
for chr in {1..22} X; do
    curl -X POST "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v1/jobs" \
         -H "Authorization: Bearer $(cat ~/topmedKey)" \
         -F "vcf=@${OUT}/onlyQc_chr${chr}.vcf.gz" \
         -F "build=GRCh38" \
         -F "population=auto" \
         -F "email=your@email.edu"
done
```

**Output**: Job IDs for tracking

---

## 12a-download_chunk1.SLURM / 12b-download_chunk2.SLURM / 12c-download_chunk3.SLURM

**Purpose**: Download imputed results from TOPMed

**Usage**:
```bash
sbatch 12a-download_chunk1.SLURM
sbatch 12b-download_chunk2.SLURM
sbatch 12c-download_chunk3.SLURM
```

**Logic**: Poll TOPMed API for job completion, download ZIP archives

**Outputs** (to `$HBCD_IMPUTATION_DIR/imputed/`):
- `c1/chr{1-22}.dose.vcf.gz` + `.info.gz`
- `c2/chr{1-22}.dose.vcf.gz` + `.info.gz`
- `c3/chr{1-22}.dose.vcf.gz` + `.info.gz`

---

## 13a-unzip.SLURM / 13b-unzip.SLURM / 13c-unzip.SLURM / 13d-unzip.SLURM

**Purpose**: Extract TOPMed ZIP archives with password

**Usage**:
```bash
sbatch 13a-unzip.SLURM   # array 1-7
sbatch 13b-unzip.SLURM
sbatch 13c-unzip.SLURM
sbatch 13d-unzip.SLURM   # chrX
```

**Logic**:
```bash
# TOPMed ZIPs are encrypted
unzip -P "$PASSWORD" "chr{chr}.zip" -d ${OUT}/
```

---

## Future Scripts

### 14a-precompute_gp.SLURM

**Purpose**: GP precomputation per chromosome × ancestry group

**Submitted by**: `futureScripts/14-submit_gp_prep.sh`

**Ancestry groups**: WHT, BLK, AIAN, ASN, HPI, 2PLUS, OTH, UNK

**Logic**: Extract genotype probabilities for low-MAF variants from imputed VCFs

### 15-aggregate_gp.SLURM

**Purpose**: Combine per-chr GP results → final CSVs

### 16-imputation-quality-report.qmd

**Purpose**: Quarto report with imputation quality visualizations

**Contents**:
- Per-chunk SNP metrics
- Excluded SNP breakdowns
- Typed SNP chromosome counts
- Interactive QC plots

### 17-cnv-qc-report.qmd / 17-cnv-qc-report.SLURM

**Purpose**: CNV genomic profile report

**Contents**:
- Rolling median (k=51) Log R Ratio per chromosome
- B Allele Frequency
- CNV value
- Faceted by chromosome with `ggh4x::facet_wrap2`
- Colored by self-reported race
- Styled by subject type (C/M)

**Render**:
```bash
sbatch futureScripts/17-cnv-qc-report.SLURM
# Or directly:
quarto render futureScripts/17-cnv-qc-report.qmd
```

---

## 00-filterConcept.qmd

**Purpose**: Quarto document with filter concept visualization

**Contents**:
- 3-panel Venn diagrams: QC-removed ∩ Excluded ∩ onlyQc
- Pipeline flowchart
- Subject count summaries

**Render**:
```bash
quarto render 00-filterConcept.qmd
```

**Output**: `00-filterConcept.html`