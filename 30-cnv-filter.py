#!/usr/bin/env python3

# step 30: filter de-identified CNV to only include release subjects
# Reads de-identified CNV from DATA_DIR (produced by 29-cnv-deid.py),
# filters to release IIDs from keep_list.txt, and writes to release dir.
#
# Usage: python 30-cnv-filter.py [--release-dir /path/to/release/data]

import sys
import pandas as pd
from pathlib import Path
from _lib import DATA_DIR, get_release_dir

CNV_DEID_PATH = DATA_DIR / "CNV_slim_clean_deid.txt"

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

# ── read de-identified CNV ──────────────────────────────────────────────
if not CNV_DEID_PATH.exists():
    print(f"ERROR: de-identified CNV not found — run 29-cnv-deid.py first: {CNV_DEID_PATH}")
    sys.exit(1)

cnv = pd.read_csv(CNV_DEID_PATH, sep="\t", dtype=str)
before = len(cnv)
print(f"\nDe-identified CNV rows: {before}")

# ── filter to release IIDs ──────────────────────────────────────────────
cnv = cnv[cnv["sample_id"].isin(release_iids)]
print(f"  After filter: {len(cnv)} rows")

# ── write ──────────────────────────────────────────────────────────────
RELEASE_DIR.mkdir(parents=True, exist_ok=True)
out_path = RELEASE_DIR / "CNV_slim_clean.txt"
cnv.to_csv(out_path, sep="\t", index=False)
print(f"\nWrote: {out_path}")
print(f"  {before} → {len(cnv)} rows ({(before - len(cnv)) / before * 100:.1f}% removed)")
