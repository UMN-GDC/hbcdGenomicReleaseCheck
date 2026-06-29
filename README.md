# HBCD Genomics Release — De-identification & Validation

## Quick Start — De-ID then Filter

These steps produce the release dataset.  Run **Phase C** scripts in order:

```bash
# ── 1. De-identify ──────────────────────────────────────────────────
# Map raw pscid IDs → anonymous release_candid in every source file.
# This must happen BEFORE filtering so that a 3rd party can safely
# handle only de-identified data.

# Step 25 — Genotype de-ID + filter (one combined step):
python 25-filterGenotypeFiles.py

# Step 26 — PLINK2 --keep to produce release bed/bim/fam:
bash 26-run_plink_filter.sh

# Step 27 — Filter imputed VCFs:
sbatch 27-filter_imputed_vcf.SLURM

# Step 28 — CNV de-identification only (pscid → release_candid):
python 28-cnv-deid.py

# ── 2. Filter to release subjects ───────────────────────────────────
# All files are now de-identified.  Filter every derivative to only
# the IIDs listed in keep_list.txt.

# Step 29 — Filter all handoff derivatives (GRM, PCs, CNV):
python 29-filter_release_outputs.py

# ── 3. Validate ─────────────────────────────────────────────────────
# Confirm no excluded subject leaked into any output file.

# Step 30 — Exclusion validation:
python 30-validateExclusions.py

# Step 31 — CNV QC report (optional):
sbatch 31-cnv-qc-report.SLURM

# Step 32 — Run tests:
pytest tests/ -v
```

The same pattern applies to any new derivative: **de-ID first** (map pscid →
release_candid in `data_handoff/`), **then filter** (keep only release IIDs).

---

## Full pipeline

## Full pipeline

The pipeline has three phases.  All data sources live in
`/projects/standard/basu_hbcd/shared/data/` unless noted otherwise.

### Phase A — QC & derivative computation (scripts 01–16)

These scripts run on the **pre-release** data (raw pscid IDs).  They are
executed once to produce the `onlyQc` dataset and all downstream derivatives.

| Step(s) | What happens |
|---------|-------------|
| **00** | **Filter concept** — `00-filterConcept.qmd` Quarto document with 3-panel Venn diagrams showing the overlap between QC-removed, excluded, and onlyQc subjects, plus a pipeline flowchart. |
| **01–02** | **QC removal** — loads QC failures from Excel (`Exclude_Summary` sheet, `Study_ID` column) and removes them alongside control subjects.  Remaining subjects → `onlyQc.{bed,bim,fam}`. |
| **03–04** | **GRM computation** — genetic relatedness matrix via GCTA (`hbcd_gcta_grm.*`) and PLINK (`hbcd_plink_grm.rel.*`). |
| **05–06** | **PC-AiR** — ancestry-adjusted principal components on unrelated subjects using KING kinship.  Outputs: `hbcd_rsid_harmonized_pc_scores.txt`, PC-AiR RDS object, unrelated/related ID lists. |
| **07–09** | **PC-Relate** — pairwise relatedness estimation using GENESIS.  Outputs include IBD probabilities (`pcrelate_pairs`, `pcrelate_ibd`, `pcrelate_self`), kinship matrices (`kinmat_wide`, `kinmat_long`), and PC-AiR score export CSV/TSV. |
| **10–13** | **Imputation prep & submission** — prepares chromosomes 1–22 + X for the TOPMed imputation server, submits, downloads result chunks (c1/c2/c3/cX), and unzips. |
| **14** | **GP precomputation** — per-chromosome × per-ancestry group genetic-probability extraction from imputed VCFs (low-MAF variants only), split by ancestry with ancestry-appropriate race codes (WHT, BLK, AIAN, ASN, HPI, 2PLUS, OTH, UNK). |
| **15** | **GP aggregation** — combines per-chromosome GP results into final CSV files. |
| **16** | **Imputation QC report** — Quarto report with per-chunk SNP metrics, excluded-SNP breakdowns, typed-SNP chromosome counts, and embedded interactive quality control. |

All phase A outputs use the original pscid-based sample IDs and are **not yet
de-identified**.

### Phase B — De-identification (intermediate, run once)

After phase A, the pscid IDs in `onlyQc` and all derivative files are mapped to
anonymous `release_candid` integers using the identifiers crosswalk
(`release_identifiers_20260526.csv`).  This is a batch rename step that
produces de-identified versions of every file in place.

After this step:
- Every PLINK IID becomes `{release_candid}{C|M}` (e.g. `1234567890C`)
- Every derivative file (GRM, PC-AiR, PC-Relate) uses the same de-identified IDs
- The original pscid IDs are no longer present in any `data/` file

This is why step 29 can filter derivatives without re-de-identifying — they
already contain anonymous IDs.

### Phase C — De-identification & release filtering (scripts 25–32)

These scripts filter the de-identified data to the **release subject subset**
(par_visit participants minus all exclusion lists) and handle the one external
file that was never de-identified (CNV).

| # | Script | De-ID or Filter | Input | Output |
|   |--------|----------------|-------|--------|
| 25 | `filterGenotypeFiles.py` | **De-ID + Filter** | de-IDed `onlyQc.{bed,bim,fam}`, identifiers, exclusions | `keep_list.txt`, `temp.fam`, `batch.info`, `Removed_individuals.txt` |
| 26 | `run_plink_filter.sh` | Filter | `onlyQc`, `temp.fam`, `keep_list.txt` | `hbcd.{bed,bim,fam}` |
| 27 | `filter_imputed_vcf.SLURM` | Filter | imputed VCFs (c1/c2/c3/cX), `keep_list.txt` | filtered `imputed/chr*.dose.vcf.gz` |
| 28 | `cnv-deid.py` | **De-ID** | `data_handoff/CNV_slim_clean.txt`, identifiers | de-IDed `data_handoff/CNV_slim_clean_deid.txt` |
| 29 | `filter_release_outputs.py` | Filter | de-IDed handoff files, `keep_list.txt` | filtered GRM/PC-AiR/CNV files in release dir |
| 30 | `validateExclusions.py` | — | all release output files | console validation report |
| 31 | `cnv-qc-report.{qmd,SLURM}` | — | `data_handoff/CNV_slim.txt`, demographics, identifiers | HTML genomic profile report (rolling median per chr × race × type) |

## Release directory

All output goes to a release-specific directory — originals are never
modified.  Set the release tag via `HBCD_RELEASE` (default `br_21p2`):

```
export HBCD_RELEASE=br_21p2

# Output root:
#   /projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_21p2/
#   ├── data/            # filtered release files
#   │   ├── hbcd.bed/bim/fam
#   │   ├── batch.info
#   │   ├── hbcd_pcrelate_grm.{grm.id,grm.bin,grm.N.bin,gz}
#   │   ├── hbcd_pcair_32PCs_clean.tsv
#   │   ├── hbcd_pcrelate_grm_pairwise.tsv
#   │   ├── CNV_slim_clean.txt
#   │   ├── Removed_individuals.txt
#   │   └── imputed/     # per-chromosome dose VCFs
#   ├── keep_list.txt
#   └── temp.fam
```

To use a different release:
```bash
export HBCD_RELEASE=br_22p0
```

---

## Phase C, steps 25–26: Genotype filter & de-identification

These two steps together produce the release PLINK files.  Step 25 builds the
subject whitelist and the ID-remapping table; step 26 applies them via PLINK2.

### Step 25: `filterGenotypeFiles.py` — build keep list & temp.fam

Reads the de-identified source PLINK data (`onlyQc`), identifiers table, batch
info, par_visit, and exclusion lists.  Applies a **two-stage filter**:

1. **QC + control removal** — removes subjects flagged in the Excel QC sheet
   (`Exclude_Summary`, `Study_ID` column) and any control subjects.  The
   remaining set is called **onlyQc**.
2. **Exclusion list removal** — removes subjects found in any exclusion source:
   - `HBCDexclusions.csv` (pscid-level, multiple reason columns)
   - Excel `Exclude_Summary` sheet (release_candid-level)
   - Additional exclusion lists

The two overlapping filter stages produce three subject groups:
- **QC_removed** — subjects removed by stage 1 only
- **Excluded** — subjects removed by stage 2, possibly overlapping with QC_removed
- **onlyQc ∩ Excluded** — overlap between the two (see `00-filterConcept.qmd` for Venn diagram)

For the remaining (valid) subjects:

1. **Merge** — joins `.fam` rows with identifiers (by `pscid` to get
   `release_candid`) and batch info (by `release_candid`).
2. **Relationship filter** — when a `release_candid` has batch data for
   multiple relationships (C + M), keeps only the row matching the subject's
   original suffix.
3. **Output**:
   - `temp.fam` — **all** subjects (not just release) with de-identified
     FID/IID, preserving original row order (critical: must match .bed row
     count exactly to avoid PLINK size mismatch).
   - `keep_list.txt` — subjects that pass *all* filters and have non-missing
     batch metadata.  This is the final release whitelist.
   - `batch.info` — tab-delimited with columns IID, visit, plate_number.
   - `Removed_individuals.txt` — excluded IIDs for documentation.

### Step 26: `run_plink_filter.sh` — PLINK2 `--keep`

```bash
plink2 --bfile onlyQc \
       --allow-extra-chr \
       --fam temp.fam \
       --keep keep_list.txt \
       --make-bed --out hbcd
```

- `--bfile` points to the de-identified `onlyQc` data.  `--fam temp.fam`
  supplies the original row count + remapped IDs so the .bed file is read
  correctly and output gets de-identified FID/IID.
- `--keep` restricts output to the release whitelist.
- Post-processing ensures `batch.info` matches `hbcd.fam` 1:1.

> **Note:** step 26 reads `RELEASE_DIR` (not `HBCD_RELEASE`).  If using a
> non-default release:
> ```bash
> export HBCD_RELEASE=br_22p0
> export RELEASE_DIR=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}/data
> ./26-run_plink_filter.sh
> ```

---

## Phase C, step 27: Filter imputed VCFs

Runs as a SLURM array (24 tasks, one per chromosome).  For each chromosome
1–22 + X, finds the imputed VCF in one of the chunk directories (c1/c2/c3/cX),
uses `bcftools view --samples-file` with the release IID whitelist to subset
samples, and writes a filtered VCF + tabix index to `release/data/imputed/`.
Also copies the corresponding `.info.gz` file alongside.

```bash
sbatch 27-filter_imputed_vcf.SLURM
```

---

## Phase C, step 28: De-identify CNV calls (de-ID only)

CNV arrives with raw pscid-level IDs from `data_handoff/` (external source).
This step maps pscid → `release_candid` **without filtering**, so a 3rd party
can handle only de-identified data downstream.

Source file: `data_handoff/CNV_slim_clean.txt` — sample IDs are in the format
`{array}_{channel}_{pscid}{C|M}` (e.g. `GSM0000000_Grn_12345C`).

`28-cnv-deid.py`:

1. Reads `data_handoff/CNV_slim_clean.txt`
2. Extracts pscid + relationship suffix (C/M) from each `sample_id`
3. Maps pscid → `release_candid` via the identifiers crosswalk
4. Builds de-identified IIDs (`{release_candid}{C|M}`)
5. Writes de-identified `CNV_slim_clean_deid.txt` back to `data_handoff/`

```bash
conda run -n python python 28-cnv-deid.py
```

---

## Phase C, step 29: Filter all handoff derivatives to release subjects

All source files have been de-identified (by Phase B or step 28) and contain
**all** QC-passing subjects.  `29-filter_release_outputs.py` reads each
derivative from `data_handoff/`, keeps only rows whose sample IDs appear in
the release IID whitelist (`keep_list.txt`), and writes filtered copies to
the release directory.  **Originals are never modified.**

Filtering methods per file type:

| File type | Filter column(s) | Helper |
|-----------|---|--------|
| PC-Relate GRM binary (`hbcd_pcrelate_grm.grm.*`) | IID | `filter_binary_grm()` |
| PC-Relate GRM text (`hbcd_pcrelate_grm.gz`) | IID1, IID2 | `filter_grm_text_gz()` |
| PC-AiR 32 PCs (`hbcd_pcair_32PCs_clean.tsv`) | `subject_id` | `filter_text()` |
| PC-Relate pairwise (`hbcd_pcrelate_grm_pairwise.tsv`) | `ID1`, `ID2` | `filter_csv()` |
| CNV (`CNV_slim_clean.txt`) | `sample_id` | inline filter |

```bash
conda run -n python python 29-filter_release_outputs.py
```

---

## Phase C, step 31: CNV genomic profile report

A Quarto document rendered on SLURM that produces
an HTML report of CNV quality metrics.  For each CNV probe, the report plots
rolling-median (k = 51) Log R Ratio, B Allele Frequency, and CNV value across
each chromosome, faceted by chromosome with independent scales via
`ggh4x::facet_wrap2`.  Lines are colored by self-reported race and styled by
subject type (C = child, M = mother).

```bash
sbatch 31-cnv-qc-report.SLURM
# or render directly:
quarto render 31-cnv-qc-report.qmd
```

---

## Phase C, step 30: Validate exclusion integrity

The final quality gate (`30-validateExclusions.py`).  Reads every output file
and checks for contamination by excluded subjects:

- **hbcd.fam** — both FID (release_candid) and IID level
- **HBCDexclusions.csv** — per-column, mapped to release_candid
- **All derivative files** — IID columns checked against the excluded-IID set
  (`{exc_rc}C` / `{exc_rc}M`)
- **CNV_slim_clean.txt** — `sample_id` column
- **Cross-check** — warns if any IID is found that's neither in the release set
  nor the exclusion set (catches unexpected subjects)

```bash
conda run -n python python 30-validateExclusions.py
```

---

## Typical run sequence (Phase C only)

Assumes Phases A–B (QC, derivatives, de-ID) are complete.  Step 25 and 26 are
required first; steps 27–29 are independent of each other.

```bash
# 0. Environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate python

# 1. Build keep list + temp.fam (de-ID + filter logic)
python 25-filterGenotypeFiles.py

# 2. PLINK2 --keep → release bed/bim/fam
RELEASE_DIR=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}/data \
  ./26-run_plink_filter.sh

# 3. Filter imputed VCFs (SLURM array, 24 tasks — chromosomes 1–22 + X)
sbatch 27-filter_imputed_vcf.SLURM

# 4. De-identify CNV (the one file that still has raw pscids)
python 28-cnv-deid.py

# 5. Filter all de-identified handoff derivatives (GRM, PCs, CNV)
python 29-filter_release_outputs.py

# 6. Validate every output for excluded-subject contamination
python 30-validateExclusions.py

# 7. Run full test suite
python -m pytest tests/ -v
```

Steps 3–4 are de-identification; step 5 filters everything.  Steps 6–7 should run last.

---

## When the exclusion list changes

If `HBCDexclusions.csv`, the Excel QC sheet, or any other exclusion source is
updated, **Phase A and B do not need to be re-run** — the underlying
derivatives (GRM, PC-AiR, PC-Relate, imputation VCFs) still contain all
QC-passing subjects and remain valid.  Only the release subject subset changes,
so only Phase C scripts need to be re-executed:

```bash
# 1. Re-build keep list with updated exclusions
python 25-filterGenotypeFiles.py

# 2. Re-filter PLINK genotypes
RELEASE_DIR=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_${HBCD_RELEASE}/data \
  ./26-run_plink_filter.sh

# 3–5. Re-filter everything (independent of each other)
sbatch 27-filter_imputed_vcf.SLURM          # VCFs
python 28-cnv-deid.py                        # CNV de-ID (pscid → release_candid)
python 29-filter_release_outputs.py          # filter all derivatives + CNV

# 6–7. Validate
python 30-validateExclusions.py
python -m pytest tests/ -v
```

Step 28 must re-run because the CNV source file (`CNV_slim_clean.txt`) always
retains raw pscid IDs in `data_handoff/` — the de-identified copy is rebuilt
fresh from the original pscid-level data each time, then step 29 re-filters
everything together.

---

## Validation tests

### `tests/test_release_data.py`

De-identification integrity, row counts, filter correctness, and per-file
IID checks for release derivative outputs:

| Test | What it checks |
|------|----------------|
| `test_no_fid_zero_in_output` | No unmatched (FID=0) subjects in release |
| `test_all_iids_match_deid_pattern` | Every IID matches `^\d{10}[CM]$` |
| `test_all_fids_are_positive_10_digit_integers` | FIDs are 10-digit integers |
| `test_iid_first_10_digits_equal_fid` | IID prefix = FID |
| `test_relation_is_C_or_M` | Suffix is C or M |
| `test_fid_length_is_10` | Strict 10-digit FID check |
| `test_iid_length_is_11` | IID length = 11 |
| `test_fid_is_numeric` | FID is numeric (not raw ID leak) |
| `test_temp_fam_row_count_matches_onlyqc` | temp.fam rows = onlyQc.fam rows |
| `test_fam_and_batch_order_matches` | fam / batch IID sets match |
| `test_batch_info_iids_are_deidentified` | batch.info IIDs match pattern |
| `test_all_output_iids_in_par_visit` | All IIDs pass par_visit + exclusion filter |
| `test_filter_correctness` | Re-derives expected subject set, confirms match |
| `test_variant_count_preserved` | .bim variant count unchanged |
| `test_pcrelate_grm_ids_are_release` | All IIDs in pcrelate_grm.id are in release set |
| `test_pcrelate_grm_dimensions` | pcrelate_grm.id rows = hbcd.fam rows |
| `test_pcair_32pcs_clean_ids_are_release` | All subject_ids in 32PCs_clean are in release set |
| `test_pcrelate_grm_pairwise_ids_are_release` | All ID1/ID2 in pairwise are in release set |
| `test_cnv_slim_clean_ids_are_release` | All sample_ids in CNV are in release set |

### `tests/test_exclusions.py` — 4 tests

| Test | What it checks |
|------|----------------|
| `test_hbcdcsv_excluded_pscids_absent` | No HBCDexclusions.csv pscid leaks into hbcd.fam |
| `test_each_hbcdcsv_exclusion_reason_individually` | Each exclusion reason checked separately |
| `test_excel_excluded_release_candids_absent` | No Excel-excluded RC leaks |
| `test_removed_individuals_absent` | No Removed_individuals.txt IID leaks |

```bash
conda run -n python python -m pytest tests/ -v
```

---

## Dependencies

- Python ≥ 3.10 with `pandas`, `numpy`, `openpyxl`, `pytest`, `matplotlib`
- PLINK 2.00 (alpha) — `module load plink/2.00-alpha-091019`
- bcftools — `module load bcftools`
- R ≥ 4.3 with `tidyverse`, `ggh4x`, `zoo`, `scales` for `.qmd` reports
- Conda environments: `python` (Python tools), `gdcPipeline` (R tools)
