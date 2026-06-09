#!/usr/bin/env python3

# step 5: validate no excluded subjects appear in the final hbcd.fam

import sys
import pandas as pd
from _lib import DATA_DIR, get_release_dir, load_additional_excluded_pscids

# Accept --release-dir / -r override, default from env/HBCD_RELEASE
if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

excluded = load_additional_excluded_pscids()
print(f"Excluded pscids  : {len(excluded)}")

identifiers = pd.read_csv(DATA_DIR / "release_identifiers_20260526.csv")
identifiers = identifiers[identifiers["release_candid"] != "release_candid"]
identifiers["release_candid"] = pd.to_numeric(identifiers["release_candid"])
identifiers = identifiers.dropna(subset=["release_candid"])
identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()

exc_rc = set(
    identifiers.loc[identifiers["pscid"].isin(excluded), "release_candid"]
    .dropna()
    .astype(int)
    .unique()
)
print(f"Mapped to RC     : {len(exc_rc)}")
print()

hbcd = pd.read_csv(
    RELEASE_DIR / "hbcd.fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
fam_rc = set(hbcd.loc[hbcd["FID"] != 0, "FID"].dropna().unique())

overlap = exc_rc & fam_rc
status = "OK" if not overlap else "OVERLAP"
print(f"  [{status:>7}] hbcd.fam")
print(f"           unique RC in fam  : {len(fam_rc):>6}")
print(f"           excluded RC in fam: {len(overlap):>6}")
if overlap:
    print(f"           overlapping RC    : {sorted(overlap)[:20]}")
print()

raw = pd.read_csv(DATA_DIR / "HBCDexclusions.csv")
for col in raw.columns:
    pscids = set(raw[col].dropna().astype(str).str.strip())
    pscids = {p for p in pscids if p and p != "nan"}
    col_rc = set(
        identifiers.loc[identifiers["pscid"].isin(pscids), "release_candid"]
        .dropna()
        .astype(int)
        .unique()
    )
    in_fam = col_rc & fam_rc
    s = "OK" if not in_fam else "OVERLAP"
    print(f"  [{s:>7}] {col}")
    print(f"           list size    : {len(pscids):>6}")
    print(f"           in hbcd.fam  : {len(in_fam):>6}")
    if in_fam:
        print(f"           overlapping RC: {sorted(in_fam)[:15]}")
    print()
