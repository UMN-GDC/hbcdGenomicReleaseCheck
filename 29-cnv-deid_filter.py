#!/usr/bin/env python3

# step 29: de-identify + filter CNV_slim_clean.txt to only include release subjects
# Maps pscid → release_candid via identifiers crosswalk, then filters to
# release IIDs from keep_list.txt.
#
# CNV sample_id format:  ..._<pscid><C|M>
#   e.g.  "GSM0000000_Grn_12345C" → pscid=12345, type=C
#
# Usage: python 29-cnv-deid_filter.py [--release-dir /path/to/release/data]

import sys
import pandas as pd
from pathlib import Path

from _lib import DATA_DIR, get_release_dir, load_identifiers

CNV_SOURCE_DIR = DATA_DIR.parent / "data_handoff"
CNV_INPUT = CNV_SOURCE_DIR / "CNV_slim_clean.txt"

if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

release_base = RELEASE_DIR.parent

# ── load release IIDs ───────────────────────────────────────────────────
keep = pd.read_csv(
    release_base / "keep_list.txt",
    sep=r"\s+",
    header=None,
    names=["FID", "IID"],
)
release_iids = set(keep["IID"].astype(str))
print(f"Release IIDs: {len(release_iids)}")

# ── load identifiers crosswalk ─────────────────────────────────────────
identifiers = load_identifiers()
pscid_to_rc = dict(
    zip(identifiers["pscid"].astype(str), identifiers["release_candid"].astype(int))
)
print(f"Crosswalk entries: {len(pscid_to_rc)}")

# ── read CNV file ──────────────────────────────────────────────────────
if not CNV_INPUT.exists():
    print(f"ERROR: CNV input not found: {CNV_INPUT}")
    sys.exit(1)

cnv = pd.read_csv(CNV_INPUT, sep="\t", dtype=str)
before = len(cnv)
print(f"\nCNV rows : {before}")
print(f"CNV cols : {list(cnv.columns)}")

# ── extract pscid + type from sample_id ────────────────────────────────
# sample_id:  <anything>_<anything>_<pscid><C|M>
cnv["type"] = cnv["sample_id"].str[-1]
cnv["_pscid_orig"] = cnv["sample_id"].str.split("_").str[-1].str[:-1]

# ── map pscid → release_candid ─────────────────────────────────────────
cnv["release_candid"] = cnv["_pscid_orig"].map(pscid_to_rc)
unmapped = cnv["release_candid"].isna().sum()
if unmapped:
    print(f"  Warning: {unmapped} sample(s) with unmapped pscid — dropping")
    cnv = cnv.dropna(subset=["release_candid"])

cnv["release_candid"] = cnv["release_candid"].astype(int)

# ── build de-identified IID ────────────────────────────────────────────
cnv["IID"] = cnv["release_candid"].astype(str) + cnv["type"]

# ── filter to release IIDs ──────────────────────────────────────────────
cnv = cnv[cnv["IID"].isin(release_iids)]
print(f"  After filter: {len(cnv)} rows")

# ── replace sample_id with de-identified IID ───────────────────────────
cnv["sample_id"] = cnv["IID"]
cnv = cnv.drop(columns=["type", "_pscid_orig", "release_candid", "IID"])

# ── write ──────────────────────────────────────────────────────────────
RELEASE_DIR.mkdir(parents=True, exist_ok=True)
out_path = RELEASE_DIR / "CNV_slim_clean.txt"
cnv.to_csv(out_path, sep="\t", index=False)
print(f"\nWrote: {out_path}")
print(f"  {before} → {len(cnv)} rows ({(before - len(cnv)) / before * 100:.1f}% removed)")
