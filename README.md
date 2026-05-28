# HBCD Genomics Release — De-identification & Validation

De-identifies PLINK genotype data for public release: raw subject IDs are
mapped to anonymous `release_candid` integers, filtered to subjects present in
the parent-visit (par_visit) table (minus an exclusion list), and matched
with batch metadata.

## Pipeline

### 1. Filter & de-identify (`01-filterGenotypeFiles.py`)

Reads the source PLINK data, identifiers table, batch info, par_visit, and
exclusion list.  Performs:

1. **Merge** — joins raw `.fam` subjects with identifiers (by `pscid`) and
   batch info (by `release_candid`).
2. **Relationship filter** — when a `release_candid` has batch data for
   multiple relationships (e.g. C + M), keeps only the row matching the
   subject's original relationship.
3. **Par_visit filter** — computes `valid_release_candids` as the set of
   `release_candid` values whose `candid` links to par_visit, minus any
   excluded release_candids.
4. **Output**:
   - `temp.fam` — all subjects with de-identified FID / IID (preserving
     original row order, used by plink2 `--fam`).
   - `keep_list.txt` — subjects that pass *all* filters (par_visit,
     exclusion, relationship match, non-missing batch metadata).
   - `batch.info` — tab-delimited with columns IID, visit, plate_number.
   - `Removed_individuals.txt` — excluded subjects (documentation).

### 2. PLINK2 filter (`02-run_plink_filter.sh`)

```bash
plink2 --bfile $dataPrefix \
       --allow-extra-chr \
       --set-all-var-ids 'chr@_#_\\$r_\\$a_b38' \
       --fam temp.fam \
       --keep keep_list.txt \
       --make-bed --out hbcd
```

- `--bfile` points to the **raw** (source) data so `temp.fam` subject count
  matches `.bed` exactly.
- `--keep` restricts output to validated subjects only.
- An `awk` post-processing step re-writes `batch.info` to contain exactly the
  same subjects (and order) as `hbcd.fam`, preserving the header row.

## Output files

| File | Format | Description |
|------|--------|-------------|
| `hbcd.bed` | binary | Filtered genotype matrix |
| `hbcd.bim` | text | Variant metadata (identical row count to source .bim) |
| `hbcd.fam` | space-delimited | 6 columns: FID, IID, PAT, MAT, SEX, PHENO |
| `batch.info` | tab-delimited | 3 columns: IID, plate_number (float), visit (character) |
| `Removed_individuals.txt` | 1 column | Excluded IIDs for documentation |

**Intermediate files** (`temp.fam`, `keep_list.txt`): deleted after a
successful run.

## Validation tests (`tests/test_release_data.py`)

| Test | Description | Status |
|------|-------------|--------|
| `test_fam_and_batch_order_matches` | hbcd.fam and batch.info have same IID set, row count, and order | ✓ |
| `test_deid_fid_length_is_10_or_na` | FID = release_candid (10 digits) or NA for excluded | ✓ |
| `test_deid_iid_length_is_11` | IID = 10-digit FID + 1-char relationship | ✓ |
| `test_deid_fid_is_numeric` | FID column is numeric | ✓ |
| `test_iid_first_10_digits_equal_fid` | IID prefix matches FID exactly | ✓ |
| `test_relation_is_C_or_M` | Relationship is `C` (child) or `M` (mother) | ✓ |
| `test_all_output_iids_in_par_visit` | Every output release_candid maps to a par_visit entry; none are excluded | ✓ |
| `test_filter_correctness` | Re-derives the expected subject set from scratch and verifies exact match with hbcd.fam | ✓ |
| `test_variant_count_preserved` | Source .bim and output .bim have identical variant counts | ✓ |

## Running

### Validate (fast, no pipeline re-run)

```bash
# activate environment, then:
pytest --junitxml=tests.xml
```

### Dependencies

- Python ≥ 3.10 with `pandas`, `numpy`, `openpyxl`, `pytest`
- PLINK 2.00 (alpha)
- `awk` (for batch.info post-processing)
