#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

DATA_DIR = Path("/projects/standard/basu_hbcd/shared/data")

identifiers = pd.read_csv(DATA_DIR / "release_identifiers_20260526.csv")
identifiers = identifiers[identifiers["release_candid"] != "release_candid"]
identifiers["release_candid"] = pd.to_numeric(identifiers["release_candid"])
identifiers = identifiers.dropna(subset=["release_candid"])
identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()
release_pscids = set(identifiers["pscid"].unique())
print(f"Release pscids  : {len(release_pscids)}")

# ── load exclusions (CSV: columns = exclusion reasons, values = pscids) ──
raw = pd.read_csv(DATA_DIR / "HBCDexclusions.csv")
excl_lists = {
    col: set(raw[col].dropna().astype(str).str.strip())
    for col in raw.columns
}

# ── check overlap ────────────────────────────────────────────────────────
print(f"\n{len(excl_lists)} exclusion lists\n")

for header, pscids in excl_lists.items():
    overlap = release_pscids & pscids
    n = len(overlap)
    status = "OK" if n == 0 else "OVERLAP"
    print(f"  [{status:>7}] {header}")
    print(f"           list size: {len(pscids):>6}  |  intersection: {n:>6}")
    if n > 0:
        print(f"           overlapping: {sorted(overlap)[:10]}")
    print()

n_bad = sum(1 for v in excl_lists.values() if len(release_pscids & v) > 0)
if n_bad:
    print(f"WARNING: {n_bad} list(s) overlap with release subjects")
else:
    print("All clean — zero overlap")
