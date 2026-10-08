# Phase C: Release Filtering & Validation

## Overview

Phase C runs **for each release** (e.g., `br31p2`). It filters the de-identified data from Phase B to the release subject subset (par_visit participants minus all exclusions) and validates the results.

## Inputs

| Source | Files | Status |
|--------|-------|--------|
| Phase B (de-identified) | `onlyQc_deid.*`, GRM, PC-AiR, PC-Relate | De-IDed, all QC subjects |
| `data_handoff/` | CNV files (pscid-based) | External, needs de-ID |
| Imputation | `$HBCD_IMPUTATION_DIR/imputed/c{1,2,3}/` | Already de-IDed |
| Identifiers | `release_identifiers_20260628.csv` | Crosswalk |
| Exclusions | `HBCDexclusions.csv`, Excel QC sheet | Multi-source |
| par_visit | `par_visit_data_br21_1.tsv` | Participation filter |

## Steps

### Step 25: `25-filterGenotypeFiles.py` — Genotype De-ID + Filter

**Purpose**: Build release whitelist (`keep_list.txt`) and de-identified `temp.fam`

```bash
python 25-filterGenotypeFiles.py
```

**Logic**:
1. **Load data**: identifiers, batch info, par_visit, exclusions
2. **Map batch IIDs**: `{pscid}{C|M}` → `release_candid` via crosswalk
3. **Load `onlyQc.fam`**: Extract `release_candid` from IID (already de-IDed by Phase B)
4. **Merge**: `.fam` × identifiers × batch info
5. **Relationship filter**: Keep only rows where batch relationship matches original suffix
6. **Two-stage filter**:
   - Stage 1 (QC): par_visit subjects only
   - Stage 2 (Exclusions): Remove HBCDexclusions.csv + Excel exclusions
7. **Outputs**:
   - `temp.fam` — ALL subjects, de-IDed, preserves .bed row count
   - `keep_list.txt` — Release whitelist (FID + IID)
   - `GDA/batch.info` — Batch metadata for release subjects
   - `GDA/removed_individuals.txt` — Excluded IIDs

**Key Data Structures**:
```python
# valid_release_candids = (par_visit ∩ identifiers) - exclusions
# keep_list = valid subjects with non-missing visit + plate_number
```

### Step 26: `26-run_plink_filter.sh` — PLINK2 Filter

**Purpose**: Apply `--keep` to produce release PLINK files

```bash
./26-run_plink_filter.sh
```

**Command**:
```bash
plink2 --bfile onlyQc \
       --allow-extra-chr \
       --fam temp.fam \
       --keep keep_list.txt \
       --make-bed --out GDA/merged_chroms
```

**Critical**: `--fam temp.fam` supplies de-IDed FID/IID + original row count so .bed reads correctly.

**Output**: `GDA/merged_chroms.{bed,bim,fam}`

### Step 27: `27-filter_imputed_vcf.SLURM` — Filter Imputed VCFs

**Purpose**: Filter TOPMed imputed VCFs to release subjects

```bash
sbatch 27-filter_imputed_vcf.SLURM
# Wait: squeue -u $USER --name=filter_imputed
```

**SLURM Array**: 24 tasks (chr 1-22 + X)

**Per-chromosome logic**:
1. Find source VCF in `c1/`, `c2/`, or `c3/` subdirs
2. VCF samples already de-identified (`{release_candid}{C|M}`)
3. `bcftools query -l` → extract sample names
4. Intersect with `keep_list.txt` IIDs
5. `bcftools view --samples-file` → filtered VCF
6. Output: `imputed/chr{chr}_dose.vcf.gz` + `.tbi` + copy `.info.gz`
7. Task 0: Copy batch QC statistics files (batch*-prefixed)

**Outputs**: `imputed/chr{1-22}_dose.vcf.gz`, `imputed/chrX_dose.vcf.gz`, indices, INFO, QC stats

### Step 28: `28-cnv-deid.py` — CNV De-identification

**Purpose**: De-identify CNV files (only file with raw pscids)

```bash
python 28-cnv-deid.py
```

**Two inputs**:

1. **CNV slim clean**: `data_handoff/CNV_slim_clean.txt`
   - Sample IDs: `{array}_{channel}_{pscid}{C|M}` (e.g., `GSM0000000_Grn_12345C`)
   - Extract pscid + relationship → map to release_candid
   - Output: `staging/cnv/CNV_slim_deid.txt` (de-IDed, all QC subjects)

2. **CNV bookmarks**: `data_handoff/HBCD_CNV_bookmark_metrics_clean.csv`
   - Sample IDs already de-identified (`{release_candid}{C|M}`)
   - Copy as-is: `staging/cnv/CNV_bookmarks_deid.csv`

**Note**: Runs de-ID only (no filtering). Step 29 will filter to release.

### Step 29: `29-filter_release_outputs.py` — Filter All Derivatives

**Purpose**: Filter all de-identified handoff derivatives to release IIDs

```bash
python 29-filter_release_outputs.py
```

**Reads**: `keep_list.txt` → `release_iids` set
**Sources**: `data_handoff/` (de-IDed by Phase B) + `staging/cnv/` (de-IDed by step 28)
**Writes**: `genotype_microarray/genesis/` + `genotype_microarray/cnv/`

| Output | Source | Filter Column(s) | Method |
|--------|--------|------------------|--------|
| `genesis/pcrelate_relatedness.grm.{id,bin,N.bin}` | `hbcd_pcrelate_grm.grm.*` | IID | Binary GRM subset |
| `genesis/pcrelate_relatedness.grm.gz` | `hbcd_pcrelate_grm.grm.gz` | IID1, IID2 | Text GRM index mapping |
| `genesis/pcrelate_relatedness.tsv` | `hbcd_pcrelate_grm_pairwise.tsv` | ID1, ID2 | TSV filter |
| `genesis/pcair_weights.tsv` | `hbcd_pcair_32PCs_clean.tsv` | subject_id | TSV filter |
| `cnv/CNV_slim.txt` | `staging/cnv/CNV_slim_deid.txt` | sample_id | TSV filter |
| `cnv/CNV_bookmarks.csv` | `staging/cnv/CNV_bookmarks_deid.csv` | sample_id | CSV filter |

**Cross-form consistency check**: Verifies CNV and bookmarks have matching IID sets.

### Step 30: `30-validateExclusions.py` — Validate Exclusion Integrity

**Purpose**: Final quality gate — confirm no excluded subject in any output

```bash
python 30-validateExclusions.py
```

**Checks**:

1. **`GDA/merged_chroms.fam`** — FID (release_candid) + IID level
2. **`HBCDexclusions.csv`** — Per-column, mapped to release_candid
3. **All derivatives** — IID columns vs excluded-IID set
4. **CNV files** — `sample_id` columns
5. **Cross-check** — Warns of unexpected IIDs (neither release nor excluded)

**Exclusion sources**:
- `HBCDexclusions.csv` — PSCID-level, multi-column reasons
- Excel `Exclude_Summary` — release_candid-level
- `GDA/removed_individuals.txt` — IID-level (pipeline output)

**Log output**: Detailed overlap report per file + summary.

## Test Suite

```bash
python -m pytest tests/ -v
```

### `tests/test_release_data.py` (31 tests)

| Test | Purpose |
|------|---------|
| `test_no_fid_zero_in_output` | No FID=0 in release |
| `test_all_iids_match_deid_pattern` | IIDs match `^\d{10}[CM]$` |
| `test_all_fids_are_positive_10_digit_integers` | FIDs are 10-digit integers |
| `test_iid_first_10_digits_equal_fid` | IID prefix = FID |
| `test_relation_is_C_or_M` | Suffix is C or M |
| `test_fid_length_is_10` | Strict 10-digit FID |
| `test_iid_length_is_11` | IID length = 11 |
| `test_fid_is_numeric` | FID numeric (no raw ID leak) |
| `test_temp_fam_row_count_matches_onlyqc` | temp.fam rows = onlyQc.fam rows |
| `test_fam_and_batch_order_matches` | fam/batch IID sets + order match |
| `test_batch_info_iids_are_deidentified` | batch.info IIDs de-IDed |
| `test_all_output_iids_in_par_visit` | All IIDs pass par_visit + exclusions |
| `test_filter_correctness` | Re-derive expected set, confirm match |
| `test_variant_count_preserved` | .bim variant count unchanged |
| `test_pcrelate_grm_ids_are_release` | GRM .id IIDs are release |
| `test_pcrelate_grm_dimensions` | GRM .id rows = merged_chroms.fam |
| `test_pcair_32pcs_clean_ids_are_release` | PC-AiR subject_ids are release |
| `test_pcrelate_grm_pairwise_ids_are_release` | Pairwise ID1/ID2 are release |
| `test_pcrelate_grm_gz_ids_are_release` | GRM .gz indices map to release |
| `test_cnv_slim_clean_ids_are_release` | CNV sample_ids are release |
| `test_cnv_bookmark_ids_are_release` | Bookmark sample_ids are release |
| `test_imputed_vcf_sample_ids_deidentified` | VCF samples match de-ID pattern |
| `test_imputed_vcf_subject_count_matches_release` | VCF sample count = release total |
| `test_imputed_vcf_subject_set_matches_release` | VCF IID set = release IIDs |
| `test_imputed_output_completeness` | All expected VCF outputs exist |
| `test_release_tree_completeness` | Full release directory tree exists |
| `test_imputation_qc_files_content` | QC stats files non-empty, correct columns |
| `test_imputed_vcf_genotype_concordance` | PLINK/VCF genotype concordance (chr22) |

### `tests/test_exclusions.py` (4 tests)

| Test | Purpose |
|------|---------|
| `test_hbcdcsv_excluded_pscids_absent` | No HBCDexclusions.csv leaks |
| `test_each_hbcdcsv_exclusion_reason_individually` | Per-column exclusion check |
| `test_excel_excluded_release_candids_absent` | No Excel-excluded RC leaks |
| `test_removed_individuals_absent` | No removed_individuals.txt leaks |

## Wrapper Script

```bash
./run_release_pipeline.sh
```

**Orchestrates**:
1. Activates `gdcPipeline` conda env
2. Step 25 (de-ID + filter)
3. Step 26 (PLINK2 --keep)
4. Step 27 (SLURM array + polling)
5. Step 28 (CNV de-ID)
6. Step 29 (filter derivatives)
7. Step 30 (validate exclusions)

**Environment variables** (with defaults):
```bash
export HBCD_RELEASE=br_21p3
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
```

## Re-running After Exclusion Changes

If exclusion lists change, **only Phase C re-runs**:

```bash
python 25-filterGenotypeFiles.py      # Rebuild keep_list
./26-run_plink_filter.sh              # Re-filter PLINK
sbatch 27-filter_imputed_vcf.SLURM    # Re-filter VCFs (filter-only)
python 28-cnv-deid.py                 # Re-de-ID CNV (source always pscid)
python 29-filter_release_outputs.py   # Re-filter derivatives
python 30-validateExclusions.py       # Re-validate
python -m pytest tests/ -v            # Re-test
```

## Phase C Outputs

See [Output Files](outputs.md) for complete directory structure and file formats.

## Validation Checklist

Before considering release ready:

- [ ] All 35 tests pass (`pytest tests/ -v`)
- [ ] `30-validateExclusions.py` shows no overlaps
- [ ] `keep_list.txt` count matches expected release size
- [ ] All 23 imputed VCFs present with indices
- [ ] GRM dimensions match `merged_chroms.fam`
- [ ] CNV + bookmarks IID sets consistent
- [ ] No FID=0 in `merged_chroms.fam`
- [ ] All IIDs match `^\d{10}[CM]$` pattern
- [ ] `batch.info` 1:1 with `merged_chroms.fam`