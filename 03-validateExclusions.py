#!/usr/bin/env python3

from _lib import load_identifiers, load_additional_excluded_pscids

identifiers = load_identifiers()
identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()
release_pscids = set(identifiers["pscid"].unique())
print(f"Release pscids  : {len(release_pscids)}")

excluded = load_additional_excluded_pscids()
print(f"Excluded pscids : {len(excluded)}")

overlap = release_pscids & excluded
print(f"Overlap         : {len(overlap)}")
if overlap:
    print(f"Overlapping     : {sorted(overlap)[:20]}")
else:
    print("All clean — zero overlap")
