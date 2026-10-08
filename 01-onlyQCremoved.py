#!/usr/bin/env python3
"""
Step 01: QC Removal — Map HBCD.fam IIDs to release_candid, Produce QC_removed.txt

Maps raw HBCD genotype IIDs to anonymous release_candid identifiers and identifies
QC failures + control subjects for removal. This is the first step of Phase A.

Pipeline Phase: A (QC & Derivative Computation) — Run FIRST
Prerequisites: Raw HBCD transfer data, identifiers crosswalk, QC Excel

Inputs:
    - $HBCD_DATA_DIR/../HST_HBCD_Transfer_July2025/HBCD/HBCD.fam: Raw PLINK .fam
      IID format: {prefix}_{pscid}{C|M} (pscid extracted as 9-char substring)
    - $HBCD_DATA_DIR/release_identifiers_20260628.csv: PSCID -> release_candid crosswalk
    - $HBCD_DATA_DIR/HBCD_genetics_QC1_missing_race_LORIS.xlsx: QC exclusions (Exclude_Summary)

Outputs (to $HBCD_DATA_DIR/):
    - HBCD_remapped.fam: Remapped .fam with release_candid-based FID/IID
    - QC_removed.txt: PSCIDs to remove (QC failures from Excel)
    - Removed_controls.txt: Control subject IIDs (non-standard IID format)

Processing:
    1. Load HBCD.fam, extract pscid (9 chars) + relationship (C/M) from IID
    2. Merge with identifiers crosswalk (pscid -> release_candid)
    3. Assign FID = release_candid (or negative for unmapped), IID = FID + rel
    4. Disambiguate duplicate IIDs with _1, _2 suffixes
    5. Identify controls: IIDs not matching ^\d{10}[CM]$ pattern
    6. Load QC Excel, extract Study_ID -> pscid + rel, map to release_candid
    7. Write QC_removed.txt (QC failures) and Removed_controls.txt

Usage:
    python 01-onlyQCremoved.py

Environment Variables:
    HBCD_DATA_DIR: Base data directory (default: /projects/standard/basu_hbcd/shared/data)

Example:
    export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
    python 01-onlyQCremoved.py
"""

import re
from pathlib import Path

import pandas as pd

from _lib import DATA_DIR, load_identifiers

HBCD_DIR = DATA_DIR / ".." / "HST_HBCD_Transfer_July2025" / "HBCD"

fam = pd.read_csv(
    HBCD_DIR / "HBCD.fam",
    sep=r"\s+",
    header=None,
)
fam.columns = ["bad_fid", "bad_iid", "PAT", "MAT", "Sex", "Pheno"]
fam["iid"] = fam["bad_iid"].str[21:]
fam["pscid"] = fam["iid"].str[:9]
fam["rel"] = fam["iid"].str[9]

identifiers = load_identifiers()

fam = fam.merge(
    identifiers[["pscid", "release_candid"]].drop_duplicates(subset="pscid"),
    on="pscid",
    how="left",
)

# -- assign FID/IID --
has_rc = fam["release_candid"].notna()

fam["FID"] = fam["release_candid"].fillna(0).astype(int)
fam.loc[~has_rc, "FID"] = range(-1, -((~has_rc).sum()) - 1, -1)

fam["IID_final"] = fam["FID"].astype(str) + fam["rel"]

# -- disambiguate any duplicate IIDs --
dups = fam["IID_final"].duplicated(keep=False)
if dups.any():
    counter = {}
    for idx in fam[dups].index:
        iid = fam.at[idx, "IID_final"]
        counter[iid] = counter.get(iid, 0) + 1
        fam.at[idx, "IID_final"] = f"{iid}_{counter[iid]}"

# -- identify controls (IID not matching 10 digits + C/M) --
pat = re.compile(r"^\d{10}[CM]$")
controls = fam[~fam["IID_final"].str.match(pat)]
controls[["FID", "IID_final"]].to_csv(
    str(DATA_DIR / "Removed_controls.txt"),
    sep=" ",
    index=False,
    header=False,
)

# -- identify QC failures from xlsx --
qc_xlsx = pd.read_excel(
    Path("../data/HBCD_genetics_QC1_missing_race_LORIS.xlsx"),
    sheet_name="Exclude_Summary",
)
qc_raw = qc_xlsx["Study_ID"].dropna().astype(str).str.strip()
qc_df = qc_raw.str.extract(r"(?P<pscid>.+)(?P<rel>[CM])$")
qc_df = qc_df.merge(
    identifiers[["pscid", "release_candid"]].drop_duplicates(subset="pscid"),
    on="pscid", how="inner",
)
qc_iids = qc_df["release_candid"].astype(int).astype(str) + qc_df["rel"]
qc = fam.loc[fam["IID_final"].isin(qc_iids), ["FID", "IID_final"]]
qc.to_csv(str(DATA_DIR / "QC_removed.txt"), sep=" ", index=False, header=False)

print(f"  QC failures (xlsx)     : {len(qc)}")
print(f"  Written: {DATA_DIR}/QC_removed.txt")

fam[["FID", "IID_final", "PAT", "MAT", "Sex", "Pheno"]].to_csv(
    str(DATA_DIR / "HBCD_remapped.fam"),
    sep=" ",
    index=False,
    header=False,
    na_rep="NA",
)

n_ok = fam["release_candid"].notna().sum()
n_miss = fam["release_candid"].isna().sum()
n_dup = dups.sum()
n_ctrl = len(controls)
print(f"  HBCD.fam rows (total)   : {len(fam)}")
print(f"  Mapped to release_candid: {n_ok}")
print(f"  No identifier match     : {n_miss}")
print(f"  Controls (non-10d+C/M)  : {n_ctrl}")
print(f"  Duplicate IIDs fixed    : {n_dup}")
print(f"  Written: {DATA_DIR}/HBCD_remapped.fam")
print(f"  Written: {DATA_DIR}/Removed_controls.txt")
