#!/usr/bin/env python3

"""Load HBCD_exclusions.xlsx, parse stacked lists of pscids, and check
overlap against the release subjects (release_candid-level).
"""

from pathlib import Path

from _lib import (
    DATA_DIR,
    load_identifiers,
    parse_additional_exclusion_lists,
    ADDITIONAL_EXCLUSIONS_FILE,
)


def main():
    identifiers = load_identifiers()
    identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()
    release_pscids = set(identifiers["pscid"].unique())
    print(f"Release pscids (from identifiers)  : {len(release_pscids)}")

    path = DATA_DIR / ADDITIONAL_EXCLUSIONS_FILE
    print(f"Loading exclusion lists from       : {path}")
    print()

    excl_lists = parse_additional_exclusion_lists(path)

    print(f"Found {len(excl_lists)} exclusion lists\n")
    all_intersections = {}

    for header, pscids in excl_lists.items():
        overlap = release_pscids & pscids
        n = len(overlap)
        all_intersections[header] = n
        status = "OK" if n == 0 else "OVERLAP"
        print(f"  [{status:>7}] {header}")
        print(f"           list size: {len(pscids):>6}  |  "
              f"intersection: {n:>6}")
        if n > 0:
            print(f"           overlapping: {sorted(overlap)[:10]}")
        print()

    n_bad = sum(1 for v in all_intersections.values() if v > 0)
    if n_bad:
        print(f"WARNING: {n_bad} list(s) have non-zero intersection with "
              f"release subjects")
    else:
        print("All exclusion lists are clean — zero overlap with release "
              "subjects")

    return all_intersections


if __name__ == "__main__":
    main()
