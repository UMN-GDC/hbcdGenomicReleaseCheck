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

# ── load exclusions ──────────────────────────────────────────────────────
raw = pd.read_excel(
    DATA_DIR / "HBCD_exclusions20250526.xlsx",
    sheet_name=0,
    skiprows=13,
    header=None,
)
col = raw.iloc[:, 0].dropna().astype(str).str.strip()

excl_lists = {}
header = None
items = []

for val in col:
    if not val or val == "nan":
        continue
    if val.isdigit() and len(val) <= 12:
        if header is not None:
            items.append(val)
    else:
        if header is not None and items:
            excl_lists[header] = set(items)
        header = val
        items = []

if header is not None and items:
    excl_lists[header] = set(items)

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
