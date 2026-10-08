# Troubleshooting Guide

## Common Issues by Pipeline Step

### Step 25: `25-filterGenotypeFiles.py`

#### `FileNotFoundError: onlyQc.fam`
**Cause**: Phase B de-identification not run, or wrong `$HBCD_DATA_DIR`
**Fix**: 
```bash
export HBCD_DATA_DIR=/correct/path
# Ensure Phase B completed: onlyQc.fam should have de-IDed IIDs
head $HBCD_DATA_DIR/onlyQc.fam
```

#### `KeyError: 'release_candid'` in identifiers merge
**Cause**: `release_identifiers_20260628.csv` missing or wrong format
**Fix**:
```bash
# Check file exists and has correct columns
head $HBCD_DATA_DIR/release_identifiers_20260628.csv
# Should have: pscid, candid, release_candid, ...
```

#### `ValueError: could not convert string to float: 'release_candid'`
**Cause**: Header row not skipped in identifiers file
**Fix**: `_lib.py` handles this with `.query("release_candid != 'release_candid'")` — ensure using latest `_lib.py`

#### `keep_list.txt` empty or too small
**Cause**: 
- par_visit file wrong version
- Exclusion lists too aggressive
- Batch info mapping failed
**Fix**:
```bash
# Check par_visit
wc -l $HBCD_DATA_DIR/par_visit_data_br21_1.tsv
# Check exclusions
wc -l $HBCD_DATA_DIR/HBCDexclusions.csv
# Check batch mapping
python -c "
import pandas as pd
batch = pd.read_csv('$HBCD_DATA_DIR/batch.info', sep=r'\s+')
print(batch.head())
print(batch['IID'].str.match(r'^[A-Za-z]{5}\d{4}[CM]\$').sum(), 'valid IIDs')
"
```

#### `temp.fam` row count ≠ `onlyQc.fam`
**Cause**: Duplicate handling or merge issues
**Fix**: Check step 25 logs for "Duplicate IIDs removed" and "Relationship mismatch removed" messages

---

### Step 26: `26-run_plink_filter.sh`

#### `PLINK2 --fam temp.fam` error: "Sample size mismatch"
**Cause**: `temp.fam` row count ≠ `onlyQc.bed` sample count
**Fix**: Re-run step 25 — `temp.fam` MUST preserve original row order and count

#### `--keep keep_list.txt` produces empty output
**Cause**: `keep_list.txt` IIDs don't match `temp.fam` IIDs
**Fix**:
```bash
# Check IID format consistency
awk '{print $2}' temp.fam | head
awk '{print $2}' keep_list.txt | head
# Both should be {release_candid}{C|M} format
```

#### `GDA/merged_chroms.fam` has FID=0 subjects
**Cause**: Step 25 didn't filter unmapped subjects before PLINK
**Fix**: Step 25 should filter `valid = combined[combined._valid]` before writing `keep_list.txt` — verify logic

---

### Step 27: `27-filter_imputed_vcf.SLURM`

#### SLURM jobs stuck in queue
**Cause**: Partition/account issues, resource limits
**Fix**:
```bash
# Check queue status
squeue -u $USER --name=filter_imputed
# Submit with explicit account/partition
sbatch --account=hbcd_genomics --partition=msismall 27-filter_imputed_vcf.SLURM
```

#### `bcftools view --samples-file` error: "Sample not found"
**Cause**: VCF sample IDs don't match `keep_list.txt` format
**Fix**:
```bash
# Check VCF sample format
bcftools query -l imputed/c1/chr1.dose.vcf.gz | head
# Check keep_list format
head keep_list.txt
# Both should be {release_candid}{C|M}
```

#### Missing chromosome outputs
**Cause**: Source VCF not found in c1/c2/c3
**Fix**:
```bash
# Check imputation directory structure
ls $HBCD_IMPUTATION_DIR/imputed/
ls $HBCD_IMPUTATION_DIR/imputed/c1/
# Should have chr1-22.dose.vcf.gz in one of c1/c2/c3
```

#### Job fails with "No space left on device"
**Cause**: `/tmp` or work directory full
**Fix**: Set `TMPDIR` to larger filesystem:
```bash
# In SLURM script or before sbatch
export TMPDIR=/large/scratch/$USER
```

---

### Step 28: `28-cnv-deid.py`

#### `FileNotFoundError: CNV_slim_clean.txt`
**Cause**: `data_handoff/` not populated or wrong path
**Fix**:
```bash
ls $HBCD_DATA_DIR/data_handoff/CNV_slim_clean.txt
# If missing, check Phase A CNV processing completed
```

#### `KeyError: 'pscid'` in mapping
**Cause**: CNV `sample_id` format changed
**Fix**:
```bash
# Check CNV sample_id format
head $HBCD_DATA_DIR/data_handoff/CNV_slim_clean.txt
# Expected: {array}_{channel}_{pscid}{C|M} e.g., GSM0000000_Grn_12345C
# If format differs, update extraction regex in 28-cnv-deid.py
```

#### Unmapped pscids in CNV
**Cause**: pscids in CNV not in identifiers crosswalk
**Fix**:
```bash
# Check crosswalk coverage
python -c "
import pandas as pd
cnv = pd.read_csv('$HBCD_DATA_DIR/data_handoff/CNV_slim_clean.txt', sep='\t')
ids = pd.read_csv('$HBCD_DATA_DIR/release_identifiers_20260628.csv')
cnv_pscids = set(cnv['sample_id'].str.extract(r'_([A-Z0-9]+)[CM]\$')[0].dropna())
id_pscids = set(ids['pscid'].astype(str))
print('CNV pscids not in crosswalk:', cnv_pscids - id_pscids)
"
```

---

### Step 29: `29-filter_release_outputs.py`

#### `FileNotFoundError: hbcd_pcrelate_grm.grm.id`
**Cause**: Phase B de-identification not run, handoff files missing
**Fix**: Ensure Phase B completed — all `data_handoff/` files should exist

#### GRM binary filter fails: "shape mismatch"
**Cause**: `.grm.id` row count ≠ `.grm.bin` matrix dimension
**Fix**:
```bash
# Verify GRM files consistent
wc -l data_handoff/hbcd_pcrelate_grm.grm.id
# Should match matrix dimension from .bin size
```

#### Cross-form CNV consistency warning
**Cause**: CNV and bookmarks have different subject sets
**Fix**: This is a warning — check if expected. CNV may have more subjects than bookmarks or vice versa.

---

### Step 30: `30-validateExclusions.py`

#### `[OVERLAP]` in validation output
**Cause**: Excluded subject found in release output
**Fix**:
```bash
# 1. Check which file has overlap
# 2. Re-run step 25 with current exclusions
# 3. Re-run steps 26-29
# 4. Re-validate
```

#### `WARN: N IID(s) not in release set nor excluded`
**Cause**: Subject in derivative not in keep_list.txt and not in exclusions
**Fix**: 
- Check if subject should be in release (par_visit + not excluded)
- If yes: re-run step 29 after fixing keep_list
- If no: investigate source of extra subject

---

### Test Failures

#### `test_temp_fam_row_count_matches_onlyqc` FAILED
**Cause**: Step 25 `temp.fam` row count ≠ `onlyQc.fam`
**Fix**: Re-run step 25, check for duplicate handling logic

#### `test_all_iids_match_deid_pattern` FAILED
**Cause**: Non-de-identified IIDs in output
**Fix**: 
```bash
# Find bad IIDs
python -c "
import pandas as pd
df = pd.read_csv('GDA/merged_chroms.fam', sep=r'\s+', header=None, names=['FID','IID',...])
bad = df[~df.IID.astype(str).str.match(r'^\d{10}[CM]\$')]
print(bad.IID.tolist())
"
```

#### `test_imputed_vcf_subject_set_matches_release` FAILED
**Cause**: VCF sample set ≠ release IIDs
**Fix**: Re-run step 27, check SLURM job logs for failed chromosomes

#### `test_imputed_vcf_genotype_concordance` FAILED
**Cause**: Sample/genotype mismatch between PLINK and VCF
**Fix**: 
- Indicates sample shuffling during imputation pipeline
- Check step 27 filtering logic
- Verify VCF sample order matches PLINK

---

## Environment Issues

### `conda: command not found` / `module: command not found`
**Fix**: Source initialization scripts
```bash
source ~/miniconda3/etc/profile.d/conda.sh
# For modules (MSI):
source /etc/profile.d/modules.sh
```

### `plink2: command not found`
**Fix**: Load PLINK module
```bash
module load plink/2.00-alpha-091019
# Or check PATH
which plink2
```

### `bcftools: command not found`
**Fix**: Load bcftools module
```bash
module load bcftools
```

### `Rscript: command not found`
**Fix**: Load R module
```bash
module load R/4.4.2-openblas-rocky8
```

### Permission denied on SLURM submit
**Fix**: Check SLURM account
```bash
# List valid accounts
sacctmgr show user $USER format=account%30
# Submit with account
sbatch --account=your_account script.SLURM
```

---

## Data Path Issues

### Wrong `$HBCD_DATA_DIR`
**Symptom**: FileNotFoundError for core data files
**Fix**:
```bash
# Verify path
export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
ls $HBCD_DATA_DIR/HBCD.fam
ls $HBCD_DATA_DIR/onlyQc.fam
```

### Wrong `$HBCD_IMPUTATION_DIR`
**Symptom**: Step 27 can't find imputed VCFs
**Fix**:
```bash
export HBCD_IMPUTATION_DIR=/projects/standard/basu_hbcd/shared/hbcdSandboxData
ls $HBCD_IMPUTATION_DIR/imputed/c1/chr1.dose.vcf.gz
```

### Release tag mismatch
**Symptom**: Outputs written to wrong directory
**Fix**:
```bash
export HBCD_RELEASE=br31p2
# Verify
echo $HBCD_RELEASE
python -c "from _lib import get_release_dir; print(get_release_dir())"
```

---

## Debugging Commands

### Inspect Release Directory
```bash
RELEASE_DIR=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br31p2/genotype_microarray
tree $RELEASE_DIR
```

### Check Log Files
```bash
RELEASE_BASE=/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br31p2
ls -la $RELEASE_BASE/log_*
cat $RELEASE_BASE/log_30_validateExclusions_*.txt
```

### Verify De-ID Pattern
```bash
# All IIDs should match ^\d{10}[CM]$
for f in GDA/merged_chroms.fam GDA/batch.info genesis/pcair_weights.tsv; do
    echo "=== $f ==="
    awk '{print $1}' $f | sort -u | head -20
done
```

### Verify Exclusion Compliance
```bash
# Quick exclusion check
python 30-validateExclusions.py 2>&1 | grep -E "OVERLAP|WARN"
```

### Compare Subject Counts
```bash
echo "keep_list: $(wc -l < keep_list.txt)"
echo "merged_chroms.fam: $(wc -l < GDA/merged_chroms.fam)"
echo "GRM .id: $(wc -l < genesis/pcrelate_relatedness.grm.id)"
echo "PC-AiR: $(wc -l < genesis/pcair_weights.tsv)"
for f in imputed/chr*_dose.vcf.gz; do
    echo "$f: $(bcftools query -l $f | wc -l)"
done
```

---

## Getting Help

### Log Files to Collect
When reporting issues, include:
1. `$RELEASE_BASE/log_25_filterGenotypeFiles_*.txt`
2. `$RELEASE_BASE/log_26_run_plink_filter_*.txt`
3. `$RELEASE_BASE/log_29_filter_release_outputs_*.txt`
4. `$RELEASE_BASE/log_30_validateExclusions_*.txt`
5. SLURM job outputs: `slurm-<JOB_ID>.out`
6. Test output: `python -m pytest tests/ -v 2>&1 | tail -50`

### Key Contacts
- **HBCD Genomics Team**: hbcd-genomics@umn.edu
- **HDCC Liaison**: Jen Zink (for release workflow)
- **Genomics WG Chair**: For scientific questions

### Useful Resources
- [Pipeline README](https://github.com/hbcd-genomics/hbcd-genomic-release)
- [HBCD Processing Documentation](https://hbcd-cbrain-processing.readthedocs.io/)
- [NMIND Standards](https://www.nmind.org/standards-checklist/)