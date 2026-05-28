#!/usr/bin/env python3
"""Filter a GCTA-format binary GRM to subjects in par_visit minus the
Excel exclusion list.

Reads ``.grm.bin`` + ``.grm.id``, identifies subjects whose release_candid
appears in (par_visit \\ excluded), and writes a smaller GRM.

Usage
-----
python futureScripts/03-filter_grm.py \\
    --grm-prefix /path/to/full_grm \\
    --out-prefix /path/to/filtered_grm

Based on MASH's ReadGRMBin()
(https://github.com/anomalyco/MASH/blob/main/src/Estimate/data_input/load_data.py).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _lib import (
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
)


def read_grm_ids(prefix):
    path = f"{prefix}.grm.id"
    ids = pd.read_csv(path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    print(f"  GRM subjects: {len(ids)}")
    return ids


def read_grm_bin(prefix, ids):
    path = f"{prefix}.grm.bin"
    n = len(ids)
    dt = np.dtype("f4")
    raw = np.fromfile(path, dtype=dt)
    expected = n * (n + 1) // 2
    if len(raw) != expected:
        raise ValueError(
            f"{path}: expected {expected} elements for N={n}, got {len(raw)}"
        )
    GRM = np.zeros((n, n), dtype=dt)
    tri = np.tril_indices(n)
    GRM[tri] = raw
    GRM = GRM + GRM.T - np.diag(np.diag(GRM))
    return GRM


def write_grm(ids, GRM, prefix):
    n = len(ids)
    ids.to_csv(f"{prefix}.grm.id", sep="\t", header=False, index=False)
    tri = np.tril_indices(n)
    GRM[tri].tofile(f"{prefix}.grm.bin")
    print(f"  Written: {prefix}.grm.id  ({n} subjects)")
    print(f"  Written: {prefix}.grm.bin  ({GRM[tri].nbytes / 1e6:.1f} MB)")


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GCTA binary GRM to par_visit minus exclusion"
    )
    parser.add_argument("--grm-prefix", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument(
        "--grm-id-col", default="FID",
        choices=("FID", "IID"),
        help="Column in .grm.id matching release_candid (default: '%(default)s')"
    )
    args = parser.parse_args()

    # ---- 1. read the full GRM ----------------------------------------------
    print("Reading GRM …")
    grm_ids = read_grm_ids(args.grm_prefix)
    GRM = read_grm_bin(args.grm_prefix, grm_ids)

    # ---- 2. single inclusive filter ----------------------------------------
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    par_rc = set(
        identifiers[identifiers["candid"].isin(par_candids)]["release_candid"]
        .dropna().astype(int).unique()
    )
    valid_rc = par_rc - exc_rc
    print(f"\n  Valid release_candids: {len(valid_rc)} "
          f"(par_visit={len(par_rc)}, excluded={len(exc_rc)})")

    # ---- 3. subset ---------------------------------------------------------
    grm_ids["_rc"] = pd.to_numeric(grm_ids[args.grm_id_col], errors="coerce").astype("Int64")
    keep = grm_ids["_rc"].isin(valid_rc)
    n_keep = keep.sum()
    print(f"  GRM subjects kept: {n_keep} / {len(grm_ids)}")

    if n_keep == 0:
        print("ERROR: no subjects overlap — check ID format.", file=sys.stderr)
        sys.exit(1)

    grm_ids_sub = grm_ids.loc[keep, ["FID", "IID"]].reset_index(drop=True)
    idx = np.where(keep.values)[0]
    GRM_sub = GRM[np.ix_(idx, idx)]

    print(f"\nWriting filtered GRM ({n_keep} × {n_keep}) …")
    write_grm(grm_ids_sub, GRM_sub, args.out_prefix)
    print("Done.")


if __name__ == "__main__":
    main()
