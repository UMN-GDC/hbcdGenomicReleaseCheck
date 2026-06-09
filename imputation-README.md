# HBCD Imputation Pipeline (steps 13-27)

Pre-imputation preparation and TOPMed Imputation Server submission pipeline for the HBCD study.
Input is `onlyQc.{bed,bim,fam}` (output of step 2) with release_candid-based IIDs.

---

## Pipeline Overview

### 13-14. Prepare: PLINK → Per-Chromosome VCF
Converts the onlyQc PLINK fileset into sorted, per-chromosome, bgzipped VCFv4.2 files:
- Splits by chromosome (1–22, X)
- Outputs VCFv4.2 with `--export vcf-4.2 bgz`
- Chromosome IDs encoded with `chr` prefix via `--output-chr chrM`
- Excludes palindromic SNPs
- Includes only SNPs
- Realigns to GRCh38 reference via bcftools
- X chromosome uses `--set-invalid-haploid-missing`

### 15-16. Submit: Upload to TOPMed
Uploads per-chromosome VCF files to the TOPMed Imputation Server via REST API.

### 17-20. Download
Downloads imputed ZIP archives from TOPMed.

### 21-24. Unzip
Extracts ZIP archives with password.

### 25. Precompute GP
Processes genotype probabilities for rare imputed variants.

### 26. Aggregate GP
Combines per-chromosome summaries into final CSVs.

### 27. Imputation Quality Report
Quarto document with imputation quality visualizations.

---

## Security

- All processing on safebox server (NIST SP 800-53 Moderate baseline)
- Standard TOPMed ZIP encryption
- API token stored locally in `~/topmedKey`
