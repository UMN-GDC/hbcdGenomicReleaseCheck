#!/usr/bin/env python3

# step 28: de-identify CNV_slim_clean.txt — map pscid → release_candid
# Writes de-identified CNV to release_base for downstream filtering.
# CNV sample_id format:  ..._<pscid><C|M>
#   e.g.  "GSM0000000_Grn_12345C" → pscid=12345, type=C
#
# Usage: python 28-cnv-deid.py

import sys
from pathlib import Path
import pandas as pd
from _lib import DATA_DIR, get_release_base, get_release_dir, load_identifiers

CNV_SOURCE_DIR = DATA_DIR.parent / "data_handoff"
CNV_INPUT = CNV_SOURCE_DIR / "CNV_slim_clean.txt"
CNV_DEID_DIR = get_release_dir() / "cnv"
CNV_DEID_DIR.mkdir(parents=True, exist_ok=True)
CNV_DEID_OUT = CNV_DEID_DIR / "CNV_slim_deid.txt"

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
    print(f"  Dropping {unmapped} sample(s) with unmapped pscid")
    cnv = cnv.dropna(subset=["release_candid"])

cnv["release_candid"] = cnv["release_candid"].astype(int)

# ── build de-identified IID ────────────────────────────────────────────
cnv["sample_id"] = cnv["release_candid"].astype(str) + cnv["type"]
cnv = cnv.drop(columns=["type", "_pscid_orig", "release_candid"])
print(f"  De-identified rows: {len(cnv)}")

# ── write de-identified CNV (no filtering) ────────────────────────────
cnv.to_csv(CNV_DEID_OUT, sep="\t", index=False)
print(f"\nWrote de-identified: {CNV_DEID_OUT}")
print(f"  {before} → {len(cnv)} rows ({((before - len(cnv)) / before * 100):.1f}% removed)")

# ===========================================================================
# CNV bookmark metrics (already de-identified — copy as-is)
# ===========================================================================
print("\n--- CNV bookmark metrics (already de-identified) ---")

BM_INPUT = CNV_SOURCE_DIR / "HBCD_CNV_bookmark_metrics_clean.csv"
BM_DEID_OUT = CNV_DEID_DIR / "HBCD_CNV_bookmark_metrics_clean_deid.csv"

if not BM_INPUT.exists():
    print(f"  SKIP (not found): {BM_INPUT}")
else:
    import shutil
    CNV_DEID_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(BM_INPUT, BM_DEID_OUT)
    print(f"  Copied: {BM_INPUT} → {BM_DEID_OUT}")
