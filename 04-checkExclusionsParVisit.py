#!/usr/bin/env python3

import pandas as pd
from _lib import DATA_DIR, load_par_visit_candids, load_identifiers

identifiers = load_identifiers()
identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()

# map pscid → candid
pscid_to_candid = (
    identifiers[["pscid", "candid"]]
    .dropna(subset=["candid"])
    .drop_duplicates(subset="pscid")
    .set_index("pscid")["candid"]
    .astype(int)
    .to_dict()
)

par_candids = load_par_visit_candids()
print(f"par_visit candids  : {len(par_candids)}")
print()

raw = pd.read_csv(DATA_DIR / "HBCDexclusions.csv")
print(f"Exclusion lists    : {len(raw.columns)}\n")

for col in raw.columns:
    pscids = set(raw[col].dropna().astype(str).str.strip())
    pscids = {p for p in pscids if p and p != "nan"}

    # map to candids
    mapped = {p: pscid_to_candid.get(p) for p in pscids}
    mapped_candids = {c for c in mapped.values() if c is not None}
    missing = {p for p, c in mapped.items() if c is None}
    overlap = mapped_candids & par_candids
    pct = len(overlap) / len(mapped_candids) * 100 if mapped_candids else 0

    status = "OK" if not overlap else "OVERLAP"
    print(f"  [{status:>7}] {col}")
    print(f"           list pscids          : {len(pscids):>6}")
    print(f"           mapped to candid     : {len(mapped_candids):>6}")
    if missing:
        print(f"           unmapped pscids      : {len(missing):>6}")
    print(f"           in par_visit         : {len(overlap):>6}  ({pct:.0f}%)")
    if overlap:
        print(f"           overlapping candids  : {sorted(overlap)[:15]}")
    print()
