# Output Files

## Release Directory Structure

```
HBCD_genomics_release_${HBCD_RELEASE}/
├── genotype_microarray/                    # Main release deliverables
│   ├── GDA/                                # Genotype Data Archive
│   │   ├── merged_chroms.bed               # PLINK binary genotype data
│   │   ├── merged_chroms.bim               # Variant information
│   │   ├── merged_chroms.fam               # Sample information (de-IDed)
│   │   ├── batch.info                      # Batch metadata (visit, plate)
│   │   └── removed_individuals.txt         # Excluded subject IIDs
│   ├── genesis/                            # Genetic relationship derivatives
│   │   ├── pcair_weights.tsv               # PC-AiR 32 PC weights
│   │   ├── pcrelate_relatedness.grm.id     # GRM sample IDs
│   │   ├── pcrelate_relatedness.grm.bin    # GRM binary (lower triangle)
│   │   ├── pcrelate_relatedness.grm.N.bin  # GRM N binary (sample counts)
│   │   ├── pcrelate_relatedness.grm.gz     # GRM text (gzipped, indexed)
│   │   └── pcrelate_relatedness.tsv        # Pairwise relatedness
│   ├── cnv/                                # Copy number variants
│   │   ├── CNV_slim.txt                    # CNV calls (de-IDed, filtered)
│   │   └── CNV_bookmarks.csv               # CNV bookmark metrics
│   ├── imputed/                            # TOPMed imputed data
│   │   ├── chr{1-22}_dose.vcf.gz           # Dosage VCFs (filtered)
│   │   ├── chr{1-22}_dose.vcf.gz.tbi       # Tabix indices
│   │   ├── chrX_dose.vcf.gz                # Chr X dosage VCF
│   │   ├── chrX_dose.vcf.gz.tbi            # Chr X tabix index
│   │   ├── chr{1-22}.info.gz               # Imputation INFO scores
│   │   ├── chrX.info.gz                    # Chr X INFO scores
│   │   ├── batch1-snps-typed-only.txt      # Typed SNPs (batch 1)
│   │   ├── batch1-snps-excluded.txt        # Excluded SNPs (batch 1)
│   │   ├── batch1-chunks-excluded.txt      # Excluded chunks (batch 1)
│   │   ├── batch1-quality-control.html     # QC report (batch 1)
│   │   ├── batch2-snps-typed-only.txt      # Typed SNPs (batch 2)
│   │   ├── batch2-snps-excluded.txt        # Excluded SNPs (batch 2)
│   │   ├── batch2-chunks-excluded.txt      # Excluded chunks (batch 2)
│   │   ├── batch2-quality-control.html     # QC report (batch 2)
│   │   ├── batch3-snps-typed-only.txt      # Typed SNPs (batch 3)
│   │   ├── batch3-snps-excluded.txt        # Excluded SNPs (batch 3)
│   │   ├── batch3-chunks-excluded.txt      # Excluded chunks (batch 3)
│   │   └── batch3-quality-control.html     # QC report (batch 3)
│   ├── keep_list.txt                       # Release subject whitelist
│   ├── temp.fam                            # De-ID fam (all QC subjects)
│   └── staging/                            # Intermediate files (not release)
│       └── cnv/
│           ├── CNV_slim_deid.txt           # CNV de-IDed (pre-filter)
│           └── CNV_bookmarks_deid.csv      # Bookmarks de-IDed (pre-filter)
└── logs/                                   # Execution logs
    ├── log_25_filterGenotypeFiles_*.txt
    ├── log_26_run_plink_filter_*.txt
    ├── log_29_filter_release_outputs_*.txt
    └── log_30_validateExclusions_*.txt
```

---

## GDA/ — Genotype Data Archive

### `merged_chroms.{bed,bim,fam}`

**Description**: Release genotype data in PLINK 1.9 binary format, filtered to release subjects with de-identified IDs.

**Format**: Standard PLINK binary fileset
- `.bed` — Genotype data (bit-packed, 2 bits per genotype)
- `.bim` — Variant info: CHR, SNP, GD, BP, A1, A2
- `.fam` — Sample info: FID, IID, PAT, MAT, SEX, PHENO

**De-identification**:
- `FID` = 10-digit `release_candid` (e.g., `1234567890`)
- `IID` = `FID` + relationship suffix (`C` or `M`) (e.g., `1234567890C`)
- `PAT`/`MAT` = 0 (not used)
- `SEX` = 1/2/0 (male/female/unknown)
- `PHENO` = "NONE" (not used)

**Validation**:
```bash
# Check sample count
wc -l GDA/merged_chroms.fam

# Check variant count
wc -l GDA/merged_chroms.bim

# Verify de-ID pattern
awk '{print $2}' GDA/merged_chroms.fam | grep -E '^\d{10}[CM]$' | wc -l
```

### `batch.info`

**Description**: Batch metadata for release subjects.

**Format**: TSV with header

| Column | Type | Description |
|--------|------|-------------|
| `IID` | string | De-identified IID (`{release_candid}{C\|M}`) |
| `visit` | string | Visit identifier (e.g., `V02`, `V03`) |
| `plate_number` | string | Genotyping plate number |

**Example**:
```tsv
IID	visit	plate_number
1234567890C	V02	PLATE001
1234567890M	V02	PLATE001
2345678901C	V03	PLATE002
```

**Validation**: 1:1 match with `merged_chroms.fam` IIDs.

### `removed_individuals.txt`

**Description**: List of IIDs excluded from the release (for documentation).

**Format**: Space-separated, no header

| Column | Type | Description |
|--------|------|-------------|
| `IID` | string | De-identified IID (`{release_candid}{C\|M}`) |

**Sources**: Excel `Exclude_Summary` sheet + `HBCDexclusions.csv` mapped via identifiers.

---

## Genesis/ — Genetic Relationship Derivatives

### `pcair_weights.tsv`

**Description**: PC-AiR ancestry-adjusted principal component weights (32 PCs) for release subjects.

**Format**: TSV with header

| Column | Type | Description |
|--------|------|-------------|
| `subject_id` | string | De-identified IID (`{release_candid}{C\|M}`) |
| `PC1`–`PC32` | float | Principal component scores |

**Source**: `data_handoff/hbcd_pcair_32PCs_clean.tsv` filtered to release IIDs.

**Validation**: All `subject_id` values must be in `keep_list.txt`.

### `pcrelate_relatedness.grm.id`

**Description**: Sample ID list for binary GRM files.

**Format**: Space-separated, no header

| Column | Type | Description |
|--------|------|-------------|
| `FID` | string | 10-digit `release_candid` |
| `IID` | string | De-identified IID (`{release_candid}{C\|M}`) |

**Row count**: = number of release subjects = `merged_chroms.fam` rows.

### `pcrelate_relatedness.grm.bin`

**Description**: Genetic Relationship Matrix (GRM) in binary format (lower triangle, row-major).

**Format**: Binary float32 array, lower triangle packed
- Size: `N*(N+1)/2` float32 values where N = number of subjects
- Diagonal elements = 0.5 * (1 + F) where F = inbreeding coefficient
- Off-diagonal = genetic relatedness estimate

**Read with Python**:
```python
import numpy as np
N = 1000  # number of subjects
grm_flat = np.fromfile('pcrelate_relatedness.grm.bin', dtype=np.float32)
grm = np.zeros((N, N), dtype=np.float32)
grm[np.tril_indices(N)] = grm_flat
grm += grm.T
np.fill_diagonal(grm, np.diag(grm) / 2)
```

### `pcrelate_relatedness.grm.N.bin`

**Description**: Sample count matrix for GRM (number of variants contributing to each pairwise estimate).

**Format**: Same as `.grm.bin` but with integer counts (float32).

### `pcrelate_relatedness.grm.gz`

**Description**: GRM in text format (gzipped), with 1-based indices referencing `.grm.id`.

**Format**: Gzipped space-separated, no header

| Column | Type | Description |
|--------|------|-------------|
| `IID1` | int | 1-based index into `.grm.id` |
| `IID2` | int | 1-based index into `.grm.id` |
| `N` | int | Number of variants |
| `GRM` | float | Relatedness estimate |

**Mapping to IIDs**:
```python
import pandas as pd
grm_id = pd.read_csv('pcrelate_relatedness.grm.id', sep='\s+', header=None, names=['FID', 'IID'])
iid_list = grm_id['IID'].tolist()  # 0-based, so index i-1
```

### `pcrelate_relatedness.tsv`

**Description**: Pairwise relatedness estimates with explicit IIDs.

**Format**: TSV with header

| Column | Type | Description |
|--------|------|-------------|
| `ID1` | string | De-identified IID 1 |
| `ID2` | string | De-identified IID 2 |
| `kinship` | float | Kinship coefficient |
| `ibd0` | float | IBD0 probability |
| `ibd1` | float | IBD1 probability |
| `ibd2` | float | IBD2 probability |
| ... | | Additional PC-Relate columns |

---

## CNV/ — Copy Number Variants

### `CNV_slim.txt`

**Description**: Filtered CNV calls for release subjects.

**Format**: TSV with header (columns from source `CNV_slim_clean.txt`)

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | string | De-identified IID (`{release_candid}{C\|M}`) |
| `chr` | string | Chromosome |
| `start` | int | Start position (bp) |
| `end` | int | End position (bp) |
| `num_probes` | int | Number of probes in CNV |
| `length` | int | CNV length (bp) |
| `cnv` | float | CNV value (log R ratio) |
| `baf` | float | B Allele Frequency |
| ... | | Additional columns |

**Source**: `data_handoff/CNV_slim_clean.txt` → de-IDed by step 28 → filtered by step 29.

### `CNV_bookmarks.csv`

**Description**: CNV bookmark metrics for release subjects.

**Format**: CSV with header (columns from source)

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | string | De-identified IID (`{release_candid}{C\|M}`) |
| ... | | Bookmark metrics columns |

**Source**: `data_handoff/HBCD_CNV_bookmark_metrics_clean.csv` (already de-IDed) → filtered by step 29.

---

## Imputed/ — TOPMed Imputed Data

### `chr{1-22}_dose.vcf.gz` / `chrX_dose.vcf.gz`

**Description**: Imputed genotype dosages filtered to release subjects.

**Format**: VCF 4.2, bgzipped, tabix-indexed

**Key header fields**:
```
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Imputed dosage (0-2)">
##FORMAT=<ID=GP,Number=3,Type=Float,Description="Genotype probabilities (AA,AB,BB)">
```

**Sample IDs**: De-identified `{release_candid}{C|M}` (same as `merged_chroms.fam`)

**Validation**:
```bash
# Check sample count per chromosome
for f in imputed/chr*_dose.vcf.gz; do
    echo "$f: $(bcftools query -l $f | wc -l) samples"
done

# Verify all chromosomes have same sample set
bcftools query -l imputed/chr1_dose.vcf.gz | sort > /tmp/chr1.samples
bcftools query -l imputed/chr22_dose.vcf.gz | sort > /tmp/chr22.samples
diff /tmp/chr1.samples /tmp/chr22.samples
```

### `chr{1-22}.info.gz` / `chrX.info.gz`

**Description**: Imputation INFO scores per variant.

**Format**: Gzipped TSV with header (TOPMed standard)

| Column | Type | Description |
|--------|------|-------------|
| `SNP` | string | Variant ID |
| `CHR` | string | Chromosome |
| `POS` | int | Position (GRCh38) |
| `REF` | string | Reference allele |
| `ALT` | string | Alternate allele |
| `INFO` | float | Imputation INFO score (0-1) |
| `MAF` | float | Minor allele frequency |
| ... | | Additional columns |

### Batch QC Statistics Files

| File | Description |
|------|-------------|
| `batch{N}-snps-typed-only.txt` | SNPs used for haplotype estimation but not in reference panel |
| `batch{N}-snps-excluded.txt` | SNPs excluded from imputation |
| `batch{N}-chunks-excluded.txt` | Genomic regions excluded from imputation |
| `batch{N}-quality-control.html` | Per-chromosome-batch QC report (HTML) |

Where N = 1, 2, 3 (imputation chunks).

---

## Base Files

### `keep_list.txt`

**Description**: Release subject whitelist (FID + IID).

**Format**: Space-separated, no header

| Column | Type | Description |
|--------|------|-------------|
| `FID` | string | 10-digit `release_candid` |
| `IID` | string | De-identified IID (`{release_candid}{C\|M}`) |

**Use**: Input to PLINK2 `--keep`, bcftools `--samples-file`.

### `temp.fam`

**Description**: De-identified .fam for ALL QC-passing subjects (not just release). Preserves original .bed row count for PLINK compatibility.

**Format**: Same as `merged_chroms.fam` but includes non-release subjects.

**Critical**: Row count must exactly match `onlyQc.fam` (source .bed).

| Column | Type | Description |
|--------|------|-------------|
| `FID` | int | `release_candid` (0 if no match) |
| `IID` | string | `{release_candid}{C\|M}` or `0_unknown` |
| `PAT` | int | 0 |
| `MAT` | int | 0 |
| `SEX` | int | 1/2/0 |
| `PHENO` | string | "NONE" |

---

## Staging/ (Intermediate, Not Release Deliverables)

### `staging/cnv/CNV_slim_deid.txt`

CNV slim file after de-identification (step 28), before release filtering (step 29). Contains all QC-passing subjects.

### `staging/cnv/CNV_bookmarks_deid.csv`

CNV bookmarks after de-identification (step 28), before release filtering (step 29).

---

## Log Files

### `log_25_filterGenotypeFiles_*.txt`

Step 25 execution log: identifier merges, filter counts, output sizes.

### `log_26_run_plink_filter_*.txt`

Step 26 PLINK2 output: variant/sample counts, warnings.

### `log_29_filter_release_outputs_*.txt`

Step 29 filtering log: before/after row counts per file, cross-form consistency.

### `log_30_validateExclusions_*.txt`

Step 30 validation report: exclusion overlap checks per file, summary statistics.

---

## File Integrity Checks

```bash
# 1. All VCFs have tabix indices
for f in imputed/chr*_dose.vcf.gz; do
    test -f "${f}.tbi" || echo "MISSING: ${f}.tbi"
done

# 2. Sample counts match across VCFs
python -c "
import subprocess
counts = {}
for chr in range(1, 23):
    f = f'imputed/chr{chr}_dose.vcf.gz'
    n = len(subprocess.run(['bcftools', 'query', '-l', f], capture_output=True, text=True).stdout.strip().split())
    counts[chr] = n
print('Sample counts:', counts)
assert len(set(counts.values())) == 1, 'Mismatch!'
"

# 3. GRM dimensions match sample count
N_fam=$(wc -l < GDA/merged_chroms.fam)
N_grm=$(wc -l < genesis/pcrelate_relatedness.grm.id)
test $N_fam -eq $N_grm && echo "GRM dims OK" || echo "GRM dims MISMATCH"

# 4. All IIDs follow de-ID pattern
for f in GDA/merged_chroms.fam GDA/batch.info genesis/pcair_weights.tsv; do
    awk -F'\t' '{print $1}' $f | grep -E '^\d{10}[CM]$' > /dev/null || echo "BAD IIDs in $f"
done
```