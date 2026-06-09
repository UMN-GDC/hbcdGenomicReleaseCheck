#!/usr/bin/env python3

# step 1: remap HBCD.fam IIDs from raw format to release_candid-based IIDs
# source /projects/standard/basu_hbcd/shared/.venv/bin/activate

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

fam["release_candidMC"] = fam["release_candid"].astype(int).astype(str) + fam["rel"]

fam[["release_candid", "release_candidMC", "PAT", "MAT", "Sex", "Pheno"]].to_csv(
    str(DATA_DIR / "HBCD_remapped.fam"),
    sep=" ",
    index=False,
    na_rep="NA",
)

n_ok = fam["release_candid"].notna().sum()
n_miss = fam["release_candid"].isna().sum()
print(f"  HBCD.fam rows (total)   : {len(fam)}")
print(f"  Mapped to release_candid: {n_ok}")
print(f"  No identifier match     : {n_miss}")
print(f"  Written: {DATA_DIR}/HBCD_remapped.fam")
