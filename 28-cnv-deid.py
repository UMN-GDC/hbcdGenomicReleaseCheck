#!/usr/bin/env python3
"""
Step 28: CNV De-identification (pscid → release_candid)

De-identifies CNV call files by mapping PSCIDs to anonymous release_candid integers.
This is the ONLY de-identification step in Phase C — all other derivatives are
already de-identified by Phase B.

Pipeline Phase: C (Release Filtering & Validation)
Input Source: data_handoff/ (external CNV pipeline outputs, still pscid-based)

Inputs:
    - $HBCD_DATA_DIR/data_handoff/CNV_slim_clean.txt: CNV calls with pscid IIDs
      Format: sample_id = {array}_{channel}_{pscid}{C|M} (e.g., GSM0000000_Grn_12345C)
    - $HBCD_DATA_DIR/data_handoff/HBCD_CNV_bookmark_metrics_clean.csv: CNV bookmarks
      Already de-identified (sample_id = {release_candid}{C|M})
    - $HBCD_DATA_DIR/release_identifiers_20260628.csv: PSCID -> release_candid crosswalk

Outputs (to $RELEASE_BASE/staging/cnv/):
    - CNV_slim_deid.txt: De-identified CNV calls (all QC subjects, NO filtering)
    - CNV_bookmarks_deid.csv: Copied bookmarks (already de-identified)

Processing:
    1. CNV Slim Clean:
       - Extract pscid + relationship (C/M) from sample_id
       - Map pscid -> release_candid via crosswalk
       - Build de-identified IID: {release_candid}{C|M}
       - Drop unmapped pscids
       - Write all rows (no filtering - step 29 handles release filtering)
    2. CNV Bookmarks:
       - Already de-identified, copy as-is

Note: This step does NOT filter to release subjects. Step 29 applies the
keep_list.txt filter to both files.

Usage:
    python 28-cnv-deid.py

Environment Variables:
    HBCD_RELEASE: Release tag (default: br_21p3)
    HBCD_DATA_DIR: Base data directory (default: /projects/standard/basu_hbcd/shared/data)

Example:
    export HBCD_RELEASE=br31p2
    python 28-cnv-deid.py
"""

import sys
from pathlib import Path
import pandas as pd
from _lib import DATA_DIR, get_release_staging_dir, load_identifiers

CNV_SOURCE_DIR = DATA_DIR.parent / "data_handoff"
CNV_INPUT = CNV_SOURCE_DIR / "CNV_slim_clean.txt"
CNV_DEID_DIR = get_release_staging_dir() / "cnv"
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
cnv_unique = cnv["sample_id"].unique()
cnv_n_subjects = len(cnv_unique)
print(f"  De-identified rows: {len(cnv)}, unique subjects: {cnv_n_subjects}")

# ── write de-identified CNV (no filtering) ────────────────────────────
cnv.to_csv(CNV_DEID_OUT, sep="\t", index=False)
print(f"\nWrote de-identified: {CNV_DEID_OUT}")
print(f"  {before} → {len(cnv)} rows ({((before - len(cnv)) / before * 100):.1f}% removed)")

# ===========================================================================
# CNV bookmark metrics (already de-identified — copy as-is)
# ===========================================================================
print("\n--- CNV bookmark metrics (already de-identified) ---")

BM_INPUT = CNV_SOURCE_DIR / "HBCD_CNV_bookmark_metrics_clean.csv"
BM_DEID_OUT = CNV_DEID_DIR / "CNV_bookmarks_deid.csv"

if not BM_INPUT.exists():
    print(f"  SKIP (not found): {BM_INPUT}")
else:
    import shutil
    bm_src = pd.read_csv(BM_INPUT, sep=",", dtype=str)
    bm_unique = bm_src["sample_id"].unique()
    bm_n_subjects = len(bm_unique)
    CNV_DEID_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(BM_INPUT, BM_DEID_OUT)
    print(f"  Copied: {BM_INPUT} → {BM_DEID_OUT}")
    print(f"  Rows: {len(bm_src)}, unique subjects: {bm_n_subjects}")
    print(f"  Subject overlap with CNV slim: {len(set(cnv_unique) & set(bm_unique))} / "
          f"{cnv_n_subjects} CNV slim subjects")
