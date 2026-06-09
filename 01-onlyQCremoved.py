#!/usr/bin/env python3

# step 1: remap HBCD.fam IIDs from raw format to release_candid-based IIDs
# source /projects/standard/basu_hbcd/shared/.venv/bin/activate

import re

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
