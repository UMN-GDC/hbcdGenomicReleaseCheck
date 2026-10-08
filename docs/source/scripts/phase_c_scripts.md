# Phase C Scripts Detailed Reference

## 25-filterGenotypeFiles.py

**Purpose**: Genotype de-identification + filtering → `keep_list.txt`, `temp.fam`, `GDA/batch.info`, `GDA/removed_individuals.txt`

**Usage**:
```bash
python 25-filterGenotypeFiles.py
```

**Inputs** (via `_lib.py`):
- `$HBCD_DATA_DIR/onlyQc.fam` — De-identified PLINK .fam (Phase B output)
- `$HBCD_DATA_DIR/batch.info` — Batch metadata (IID={pscid}{C|M})
- `$HBCD_DATA_DIR/release_identifiers_20260628.csv` — PSCID → release_candid
- `$HBCD_DATA_DIR/HBCDexclusions.csv` — Additional PSCID exclusions
- `$HBCD_DATA_DIR/HBCD_genetics_QC1_missing_race_LORIS.xlsx` — Excel exclusions
- `$HBCD_DATA_DIR/par_visit_data_br21_1.tsv` — par_visit participation

**Outputs** (to `$RELEASE_BASE/` and `$RELEASE_DIR/GDA/`):
- `temp.fam` — ALL subjects, de-IDed, preserves .bed row count
- `keep_list.txt` — Release whitelist (FID + IID)
- `GDA/batch.info` — Batch metadata for release subjects
- `GDA/removed_individuals.txt` — Excluded IIDs

### Detailed Logic

```python
# 1. LOAD IDENTIFIERS & EXCLUSIONS
identifiers = load_identifiers()                    # pscid → release_candid
excluded_pscids = load_additional_excluded_pscids() # HBCDexclusions.csv
identifiers = identifiers[~identifiers.pscid.isin(excluded_pscids)]

# 2. LOAD & MAP BATCH INFO
batch = pd.read_csv("batch.info", sep=r"\s+")
# IID format: {array}_{channel}_{pscid}{C|M} → extract pscid + rel
batch["pscid"] = batch["IID"].str[:-1]
batch["relationship"] = batch["IID"].str[-1]
batch["release_candid"] = batch["pscid"].map(pscid_to_rc)

# 3. LOAD ONLYQC.FAM (already de-IDed by Phase B)
fam = pd.read_csv("onlyQc.fam", sep=r"\s+", ...)
# IID = {release_candid}{C|M}
fam["release_candid"] = pd.to_numeric(fam["IID"].str[:-1])
fam["_orig_rel"] = fam["IID"].str[-1]

# 4. MERGE: fam × identifiers × batch
combined = fam.merge(identifiers, on="release_candid", how="left")
combined = combined.merge(batch, on="release_candid", how="left")

# 5. RELATIONSHIP FILTER
# Keep only rows where batch relationship matches original suffix
rel_ok = combined["relationship"].isna() | (combined["relationship"] == combined["_orig_rel"])
combined = combined[rel_ok]

# 6. DE-IDENTIFIED IDs FOR ALL SUBJECTS
combined["new_FID"] = combined["release_candid"].fillna(0).astype(int)
combined["new_rel"] = combined["relationship"].fillna(combined["_orig_rel"])
combined["new_IID"] = combined["new_FID"].astype(str) + combined["new_rel"]
# Unmapped: FID=0, IID="0_unknown"

# 7. WRITE TEMP.FAM (preserves row order for PLINK .bed compatibility)
combined = combined.sort_values("_idx")  # Original row order
combined[["new_FID", "new_IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv("temp.fam", ...)

# 8. RELEASE FILTER: par_visit \ exclusions
par_candids = load_par_visit_candids()           # candid from par_visit
exc_rc = load_excluded_release_candids()         # Excel Exclude_Summary
valid_rc = (identifiers[identifiers.release_candid.isin(par_candids)].release_candid - exc_rc)

# 9. APPLY FILTER + METADATA CHECK
valid = combined[combined.release_candid.isin(valid_rc)]
valid = valid.dropna(subset=["visit", "plate_number"])
valid = valid.drop_duplicates(subset="new_IID")

# 10. OUTPUTS
valid[["new_FID", "new_IID"]].to_csv("keep_list.txt", sep=" ", ...)
valid[["new_IID", "visit", "plate_number"]].to_csv("GDA/batch.info", sep="\t", ...)
# removed_individuals.txt from Excel exclusions
```

### Key Data Structures

```python
# temp.fam: ALL QC subjects (preserves .bed row count)
#   FID=release_candid (0 if unmapped), IID={release_candid}{C|M}

# keep_list.txt: Release subjects ONLY
#   FID=release_candid, IID={release_candid}{C|M}

# GDA/batch.info: Release subjects with metadata
#   IID, visit, plate_number

# GDA/removed_individuals.txt: Excluded IIDs
#   IID={release_candid}{C|M}
```

---

## 26-run_plink_filter.sh

**Purpose**: PLINK2 `--keep` to produce release PLINK files

**Usage**:
```bash
export HBCD_RELEASE=br31p2
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

**Parameters**:
- `--bfile onlyQc` — De-identified `onlyQc` data (Phase B)
- `--allow-extra-chr` — Allow non-standard chr codes
- `--fam temp.fam` — Supplies de-IDed FID/IID + **original row count** (critical for .bed size match)
- `--keep keep_list.txt` — Restricts to release whitelist
- `--make-bed --out GDA/merged_chroms` — Output release PLINK files

**Post-processing** (in script):
```bash
# Ensure batch.info matches merged_chroms.fam 1:1
# Reorder batch.info to match .fam IID order
```

**Outputs**:
- `GDA/merged_chroms.{bed,bim,fam}` — Release genotype data
- `GDA/batch.info` — Reordered to match .fam

**Reads Environment**:
```bash
HBCD_RELEASE  # Used in output path construction
```

---

## 27-filter_imputed_vcf.SLURM

**Purpose**: Filter TOPMed imputed VCFs to release subjects (SLURM array)

**Usage**:
```bash
sbatch 27-filter_imputed_vcf.SLURM
# Monitor: squeue -u $USER --name=filter_imputed
```

**SLURM Configuration**:
```bash
#SBATCH --array=0-23          # 24 tasks: chr 1-22 + X (0=chr1, 22=chrX)
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --partition=msismall
#SBATCH --job-name=filter_imputed
```

**Per-Task Logic** (`task_id` 0-23):
```python
# Map task_id → chromosome
if task_id < 22:
    chr = task_id + 1      # 1-22
else:
    chr = "X"              # 23 = X

# 1. FIND SOURCE VCF
# Search order: c1/, c2/, c3/ subdirs under $HBCD_IMPUTATION_DIR/imputed/
src_vcf = find_vcf(chr, [c1, c2, c3])

# 2. EXTRACT SAMPLE NAMES
samples = bcftools query -l ${src_vcf}
# Samples already de-identified: {release_candid}{C|M}

# 3. INTERSECT WITH RELEASE IIDs
keep_list = load_keep_list()  # from keep_list.txt
keep_samples = intersection(samples, keep_list)

# 4. FILTER VCF
bcftools view --samples-file <(printf "%s\n" "${keep_samples[@]}") \
              ${src_vcf} -Oz -o ${OUT}/chr${chr}_dose.vcf.gz

# 5. INDEX
bcftools index ${OUT}/chr${chr}_dose.vcf.gz

# 6. COPY INFO FILE
cp ${src_vcf%.vcf.gz}.info.gz ${OUT}/chr${chr}.info.gz

# 7. TASK 0: COPY QC STATISTICS FILES
if task_id == 0:
    copy_batch_qc_files()  # batch*-snps-typed-only.txt, etc.
```

**Outputs** (to `$RELEASE_DIR/imputed/`):
- `chr{1-22}_dose.vcf.gz` + `.tbi`
- `chrX_dose.vcf.gz` + `.tbi`
- `chr{1-22}.info.gz`, `chrX.info.gz`
- `batch1/2/3-snps-typed-only.txt`
- `batch1/2/3-snps-excluded.txt`
- `batch1/2/3-chunks-excluded.txt`
- `batch1/2/3-quality-control.html`

---

## 28-cnv-deid.py

**Purpose**: CNV de-identification (pscid → release_candid)

**Usage**:
```bash
python 28-cnv-deid.py
```

**Inputs**:
- `$HBCD_DATA_DIR/data_handoff/CNV_slim_clean.txt` — CNV calls (pscid IIDs)
- `$HBCD_DATA_DIR/data_handoff/HBCD_CNV_bookmark_metrics_clean.csv` — Bookmarks (already de-IDed)
- `$HBCD_DATA_DIR/release_identifiers_20260628.csv` — Crosswalk

**Outputs** (to `$RELEASE_BASE/staging/cnv/`):
- `CNV_slim_deid.txt` — De-identified CNV calls
- `CNV_bookmarks_deid.csv` — Copied bookmarks

### CNV Slim Clean Processing

**Source IID format**: `{array}_{channel}_{pscid}{C|M}` (e.g., `GSM0000000_Grn_12345C`)

```python
# 1. LOAD CNV DATA
cnv = pd.read_csv("CNV_slim_clean.txt", sep="\t", dtype=str)

# 2. EXTRACT PSCID + RELATIONSHIP
# sample_id = GSM0000000_Grn_12345C
cnv["pscid"] = cnv["sample_id"].str.extract(r'_([A-Z0-9]+)[CM]$')[0]
cnv["relationship"] = cnv["sample_id"].str[-1]

# 3. MAP TO RELEASE_CANDID
identifiers = load_identifiers()
pscid_to_rc = identifiers.set_index("pscid")["release_candid"]
cnv["release_candid"] = cnv["pscid"].map(pscid_to_rc)

# 4. BUILD DE-IDENTIFIED IID
cnv["new_sample_id"] = cnv["release_candid"].astype(int).astype(str) + cnv["relationship"]

# 5. REPLACE & WRITE
cnv["sample_id"] = cnv["new_sample_id"]
cnv.drop(columns=["pscid", "relationship", "release_candid", "new_sample_id"]).to_csv(
    "staging/cnv/CNV_slim_deid.txt", sep="\t", index=False)
```

### CNV Bookmarks Processing

**Source**: Already de-identified (`{release_candid}{C|M}` format)

```python
# Simple copy - no de-ID needed
bm = pd.read_csv("HBCD_CNV_bookmark_metrics_clean.csv", dtype=str)
bm.to_csv("staging/cnv/CNV_bookmarks_deid.csv", index=False)
```

---

## 29-filter_release_outputs.py

**Purpose**: Filter all de-identified handoff derivatives to release IIDs

**Usage**:
```bash
python 29-filter_release_outputs.py [--release-dir /path]
```

**Inputs**:
- `$RELEASE_BASE/keep_list.txt` — Release IID whitelist
- `$HBCD_DATA_DIR/data_handoff/` — De-identified derivatives (Phase B)
- `$RELEASE_BASE/staging/cnv/` — De-identified CNV (step 28)

**Outputs** (to `$RELEASE_DIR/genesis/` and `$RELEASE_DIR/cnv/`):

| Output | Source | Filter Col(s) | Method |
|--------|--------|---------------|--------|
| `genesis/pcrelate_relatedness.grm.id` | `hbcd_pcrelate_grm.grm.id` | IID | `filter_binary_grm()` |
| `genesis/pcrelate_relatedness.grm.bin` | `hbcd_pcrelate_grm.grm.bin` | — | Binary subset |
| `genesis/pcrelate_relatedness.grm.N.bin` | `hbcd_pcrelate_grm.grm.N.bin` | — | Binary subset |
| `genesis/pcrelate_relatedness.grm.gz` | `hbcd_pcrelate_grm.grm.gz` | IID1, IID2 | `filter_grm_text_gz()` |
| `genesis/pcrelate_relatedness.tsv` | `hbcd_pcrelate_grm_pairwise.tsv` | ID1, ID2 | `filter_csv()` |
| `genesis/pcair_weights.tsv` | `hbcd_pcair_32PCs_clean.tsv` | subject_id | `filter_text()` |
| `cnv/CNV_slim.txt` | `staging/cnv/CNV_slim_deid.txt` | sample_id | Inline filter |
| `cnv/CNV_bookmarks.csv` | `staging/cnv/CNV_bookmarks_deid.csv` | sample_id | Inline filter |

### Filter Helpers

```python
def filter_binary_grm(in_prefix, out_prefix, id_ext=".grm.id", ...):
    """Subset binary GRM (.id, .bin, .N.bin) to release IIDs."""
    ids = pd.read_csv(in_id, sep="\t", header=None, names=["FID", "IID"])
    keep_mask = ids["IID"].isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    
    # Filter .id
    ids[keep_mask].to_csv(out_id, sep="\t", header=False, index=False)
    
    # Filter .bin (lower triangle, row-major)
    grm_flat = np.fromfile(in_bin, dtype=np.float32)
    grm_full = reconstruct_symmetric(grm_flat, n)
    grm_sub = grm_full[np.ix_(keep_idx, keep_idx)]
    grm_sub[np.tril_indices(n_keep)].tofile(out_bin)
    
    # Same for .N.bin

def filter_grm_text_gz(in_gz, in_id, out_gz):
    """Filter gzipped text GRM (1-based indices into .id)."""
    # Map indices to IIDs via .id file
    orig_ids = pd.read_csv(in_id, ...)
    keep_indices = {i+1 for i, row in orig_ids.iterrows() if row["IID"] in release_iids}
    
    # Filter rows where BOTH indices in keep_indices
    df = pd.read_csv(in_gz, sep=r"\s+", names=["IID1", "IID2", "N", "GRM"])
    mask = df["IID1"].isin(keep_indices) & df["IID2"].isin(keep_indices)
    df[mask].to_csv(out_gz, sep=" ", header=False, index=False, compression="gzip")

def filter_text(in_path, id_col, out_path):
    """Generic TSV filter by IID column."""
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    df = df[df[id_col].isin(release_iids)]
    df.to_csv(out_path, sep="\t", index=False)

def filter_csv(in_path, out_path, id_cols, sep=","):
    """Generic CSV filter by multiple IID columns."""
    df = pd.read_csv(in_path, sep=sep, dtype=str)
    mask = df[id_cols].apply(lambda c: c.isin(release_iids)).all(axis=1)
    df[mask].to_csv(out_path, sep=sep, index=False)
```

### Logging & Cross-form Check

- Logs to `$RELEASE_BASE/log_29_filter_release_outputs_*.txt`
- Prints before/after row counts per file
- Cross-form consistency: CNV vs Bookmarks IID sets

---

## 30-validateExclusions.py

**Purpose**: Final validation — confirm no excluded subject in any output

**Usage**:
```bash
python 30-validateExclusions.py [--release-dir /path]
```

**Checks**:

### 1. Load Exclusion Sets
```python
# HBCDexclusions.csv → pscids → release_candid
excluded_pscids = load_additional_excluded_pscids()
exc_rc = identifiers[identifiers.pscid.isin(excluded_pscids)].release_candid

# Excel Exclude_Summary → release_candid directly
exc_rc_excel = load_excluded_release_candids()

# Combined excluded release_candids
exc_rc = exc_rc ∪ exc_rc_excel

# Excluded IIDs (both C and M for each RC)
exc_iids = {f"{rc}C" for rc in exc_rc} ∪ {f"{rc}M" for rc in exc_rc}
```

### 2. GDA/merged_chroms.fam Check
```python
fam = pd.read_csv("GDA/merged_chroms.fam", ...)
fam_rc = set(fam[fam.FID != 0].FID.unique())
fam_iids = set(fam.IID.astype(str))

# FID-level
overlap_rc = exc_rc ∩ fam_rc

# IID-level
overlap_iid = exc_iids ∩ fam_iids
```

### 3. HBCDexclusions.csv Per-Column Check
```python
raw = pd.read_csv("HBCDexclusions.csv")
for col in raw.columns:
    col_pscids = set(raw[col].dropna().astype(str).str.strip())
    col_rc = identifiers[identifiers.pscid.isin(col_pscids)].release_candid
    overlap = col_rc ∩ fam_rc
```

### 4. Derivative Files Check (IID-level)
```python
def check_file(path, label, id_cols, reader="csv"):
    df = read_file(path, reader)
    all_ids = set()
    for col in id_cols:
        all_ids.update(df[col].dropna().astype(str))
    overlap = all_ids ∩ exc_iids
    # Also check for unexpected IIDs (not in release, not excluded)
    extra = all_ids - release_iids - exc_iids
```

**Files checked**:
- `genesis/pcrelate_relatedness.grm.id` (col 1 = IID)
- `genesis/pcair_weights.tsv` (subject_id)
- `genesis/pcrelate_relatedness.tsv` (ID1, ID2)
- `genesis/pcrelate_relatedness.grm.gz` (via .id index mapping)
- `cnv/CNV_slim.txt` (sample_id)
- `cnv/CNV_bookmarks.csv` (sample_id)

### 5. Cross-form Consistency
```python
# Compare IID sets across all derivative files
# Flag mismatches
```

### Output
- Log to `$RELEASE_BASE/log_30_validateExclusions_*.txt`
- Console summary with OK/OVERLAP per check
- Final summary statistics