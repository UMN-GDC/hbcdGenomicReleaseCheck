# Phase B: De-identification

## Overview

Phase B runs **once** after Phase A completes. It maps all `pscid` IDs in Phase A outputs to anonymous `release_candid` integers using the identifiers crosswalk (`release_identifiers_20260628.csv`).

This is a batch rename operation that produces de-identified versions of every file in place.

## Inputs

| File | Location | Description |
|------|----------|-------------|
| All Phase A outputs | `$HBCD_DATA_DIR` | Genotype, GRM, PC-AiR, PC-Relate |
| `release_identifiers_20260628.csv` | `$HBCD_DATA_DIR` | PSCID ↔ release_candid crosswalk |

## Crosswalk Format

`release_identifiers_20260628.csv`:

| Column | Type | Description |
|--------|------|-------------|
| `pscid` | string | Original PSCID (e.g., `12345`) |
| `candid` | int | Internal candid |
| `release_candid` | int | Anonymous release ID (10-digit) |
| ... | | Additional columns |

**Example**:
```csv
pscid,candid,release_candid,...
12345,1001,1234567890,...
23456,1002,2345678901,...
```

## De-identification Rules

### PLINK IIDs

**Before**: `{array}_{channel}_{pscid}{C|M}` (e.g., `GSM0000000_Grn_12345C`)
**After**: `{release_candid}{C|M}` (e.g., `1234567890C`)

- Extract `pscid` + relationship suffix (C/M) from original IID
- Map `pscid` → `release_candid` via crosswalk
- Construct new IID: `{release_candid}{relationship}`

### FID

**Before**: Same as IID (pscid-based)
**After**: 10-digit `release_candid` (e.g., `1234567890`)

### Derivative Files

All derivative files (GRM, PC-AiR, PC-Relate) use the same de-identified IIDs:
- `{release_candid}{C|M}` format consistently
- GRM `.id` files: FID = release_candid, IID = release_candid + C/M
- PC-AiR `subject_id`: release_candid + C/M
- PC-Relate `ID1`/`ID2`: release_candid + C/M

## Implementation

### Current Status

**No standalone script in repository** — Phase B is currently a manual/one-off process. The de-identification logic is embedded in:

1. **Step 25** (`25-filterGenotypeFiles.py`): Contains the core mapping logic for genotypes
2. **Step 28** (`28-cnv-deid.py`): Handles CNV de-identification separately
3. **Step 29** (`29-filter_release_outputs.py`): Assumes derivatives already de-identified

### Recommended Implementation

Create a dedicated script `deidentify_phase_a.py`:

```python
#!/usr/bin/env python3
"""
Phase B: De-identify all Phase A outputs using identifiers crosswalk.
Run once after Phase A completes, before any Phase C runs.
"""

import pandas as pd
from pathlib import Path
from _lib import DATA_DIR, load_identifiers

def deidentify_plink_fam(fam_path, identifiers, output_path):
    """De-identify PLINK .fam file."""
    fam = pd.read_csv(fam_path, sep=r'\s+', header=None,
                      names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    
    # Extract pscid + relationship from IID
    # Format: {array}_{channel}_{pscid}{C|M}
    fam['pscid'] = fam['IID'].str.extract(r'_([A-Z0-9]+)[CM]$')[0]
    fam['relationship'] = fam['IID'].str[-1]
    
    # Map to release_candid
    pscid_to_rc = identifiers.set_index('pscid')['release_candid']
    fam['release_candid'] = fam['pscid'].map(pscid_to_rc)
    
    # Build new IDs
    fam['new_FID'] = fam['release_candid'].astype(int).astype(str)
    fam['new_IID'] = fam['new_FID'] + fam['relationship']
    
    # Handle unmapped
    unmapped = fam['release_candid'].isna()
    if unmapped.any():
        fam.loc[unmapped, 'new_FID'] = '0'
        fam.loc[unmapped, 'new_IID'] = '0_' + fam.loc[unmapped, 'pscid'].fillna('unknown')
    
    # Write
    fam[['new_FID', 'new_IID', 'PAT', 'MAT', 'SEX', 'PHENO']].to_csv(
        output_path, sep=' ', index=False, header=False)
    
    return fam

def deidentify_grm_id(grm_id_path, identifiers, output_path):
    """De-identify GRM .id file."""
    grm_id = pd.read_csv(grm_id_path, sep=r'\s+', header=None,
                         names=['FID', 'IID'], dtype=str)
    # Similar mapping logic...
    pass

def deidentify_text_file(input_path, id_col, identifiers, output_path, sep='\t'):
    """Generic de-identification for TSV/CSV files with IID column."""
    df = pd.read_csv(input_path, sep=sep, dtype=str)
    # Extract pscid + rel from id_col, map, rebuild
    pass

def main():
    identifiers = load_identifiers()
    pscid_to_rc = identifiers.set_index('pscid')['release_candid']
    
    # Files to de-identify
    files = [
        ('onlyQc.fam', 'GDA/onlyQc_deid.fam'),
        ('data_handoff/hbcd_pcrelate_grm.grm.id', 'data_handoff/hbcd_pcrelate_grm_deid.grm.id'),
        ('data_handoff/hbcd_pcair_32PCs_clean.tsv', 'data_handoff/hbcd_pcair_32PCs_clean_deid.tsv'),
        ('data_handoff/hbcd_pcrelate_grm_pairwise.tsv', 'data_handoff/hbcd_pcrelate_grm_pairwise_deid.tsv'),
    ]
    
    for in_rel, out_rel in files:
        in_path = DATA_DIR / in_rel
        out_path = DATA_DIR / out_rel
        if in_path.exists():
            print(f"De-identifying {in_path} → {out_path}")
            # Call appropriate function
        else:
            print(f"SKIP: {in_path} not found")

if __name__ == '__main__':
    main()
```

## Expected Outputs After Phase B

```
$HBCD_DATA_DIR/
├── onlyQc_deid.{bed,bim,fam}                    # De-identified genotypes
├── hbcd_gcta_grm_deid.grm.{id,bin,N.bin}       # De-identified GCTA GRM
├── hbcd_plink_grm_rel_deid.rel.{id,bin}        # De-identified PLINK GRM
├── hbcd_rsid_harmonized_pc_scores_deid.txt     # De-identified PC scores
├── data_handoff/
│   ├── hbcd_pcrelate_grm_deid.grm.{id,bin,N.bin}
│   ├── hbcd_pcrelate_grm_deid.grm.gz
│   ├── hbcd_pcrelate_grm_pairwise_deid.tsv
│   ├── hbcd_pcair_32PCs_clean_deid.tsv
│   ├── CNV_slim_clean.txt                      # Still pscid (external)
│   └── HBCD_CNV_bookmark_metrics_clean.csv     # Already de-IDed
```

## Verification

```bash
# Check de-ID pattern in .fam
awk '{print $2}' onlyQc_deid.fam | head -20
# Should show: 1234567890C, 1234567890M, 2345678901C, ...

# Check no pscid remnants
grep -E 'GSM|pscid' onlyQc_deid.fam && echo "FAIL: pscid found" || echo "OK"

# Check GRM .id
head data_handoff/hbcd_pcrelate_grm_deid.grm.id
# FID (10-digit)  IID (10-digit + C/M)

# Check PC-AiR
head data_handoff/hbcd_pcair_32PCs_clean_deid.tsv
# subject_id should be 1234567890C format
```

## Notes

- Phase B should run **exactly once** after Phase A completes
- After Phase B, all `data/` files contain only anonymous IDs
- Original `pscid` IDs are **never** in any `data/` file after Phase B
- CNV files in `data_handoff/` remain pscid-based (external source) — handled in Phase C step 28
- Step 29 assumes Phase B is complete (filters de-identified derivatives)