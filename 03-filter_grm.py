#!/usr/bin/env python3
"""Filter a GCTA-format binary GRM to only subjects in the par_visit table.

Reads the full GRM (``.grm.bin`` + ``.grm.id``) into memory, identifies the
subset of subjects that appear in par_visit (via the identifiers mapping),
and writes a new, smaller GRM.

Usage
-----
python 03-filter_grm.py \\
    --grm-prefix /path/to/full_grm \\
    --par-visit ../HBCD_genomics_release_br_20p2/par_visit_data_br20.2.tsv \\
    --identifiers ../HBCD_genomics_release_br_20p2/release_identifiers_20251211.csv \\
    --out-prefix /path/to/filtered_grm

The filtering logic follows the same approach as MASH's ReadGRMBin():
https://github.com/anomalyco/MASH/blob/main/src/Estimate/data_input/load_data.py
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_grm_ids(prefix):
    """Read the ``.grm.id`` file (FID, IID), return as a DataFrame."""
    path = f"{prefix}.grm.id"
    ids = pd.read_csv(path, sep=r"\s+", header=None, names=["FID", "IID"], dtype=str)
    print(f"  GRM subjects: {len(ids)}")
    return ids


def read_grm_bin(prefix, ids):
    """Read the lower-triangle binary GRM and expand to a full symmetric matrix.

    Parameters
    ----------
    prefix : str
        Path prefix (without ``.grm.bin``).
    ids : DataFrame
        Subject IDs (used to determine N).

    Returns
    -------
    np.ndarray
        N × N symmetric matrix (dtype float32).
    """
    path = f"{prefix}.grm.bin"
    n = len(ids)
    dt = np.dtype("f4")          # GCTA stores relatedness as 4-byte float

    raw = np.fromfile(path, dtype=dt)
    expected = n * (n + 1) // 2
    if len(raw) != expected:
        raise ValueError(
            f"{path}: expected {expected} lower-triangle elements for N={n}, "
            f"got {len(raw)}"
        )

    GRM = np.zeros((n, n), dtype=dt)
    tri = np.tril_indices(n)
    GRM[tri] = raw
    GRM = GRM + GRM.T - np.diag(np.diag(GRM))
    return GRM


def write_grm(ids, GRM, prefix):
    """Write filtered GRM files (``.grm.id``, ``.grm.bin``).

    Parameters
    ----------
    ids : DataFrame
        Columns ``FID``, ``IID``.
    GRM : np.ndarray
        Symmetric N × N matrix.
    prefix : str
        Output path prefix.
    """
    n = len(ids)
    ids.to_csv(f"{prefix}.grm.id", sep="\t", header=False, index=False)
    tri = np.tril_indices(n)
    GRM[tri].tofile(f"{prefix}.grm.bin")
    print(f"  Written: {prefix}.grm.id  ({n} subjects)")
    print(f"  Written: {prefix}.grm.bin  ({GRM[tri].nbytes / 1e6:.1f} MB)")


def collect_valid_ids(par_visit_path, identifiers_path):
    """Return a set of ``release_candid`` (as int) from par_visit.

    This is the ID format expected in the GRM ``.grm.id`` file's FID column
    after the de-identification pipeline.
    """
    par_visit = pd.read_csv(par_visit_path, sep="\t")
    par_candid = par_visit["participant_id"].dropna().unique()

    ids = pd.read_csv(identifiers_path)
    ids["candid"] = pd.to_numeric(ids["candid"], errors="coerce")
    ids["release_candid"] = pd.to_numeric(ids["release_candid"], errors="coerce")

    par = ids[ids["candid"].isin(par_candid)]
    valid = par["release_candid"].dropna().astype(int).unique()
    return set(valid)


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GCTA binary GRM to par_visit subjects"
    )
    parser.add_argument(
        "--grm-prefix", required=True,
        help="Prefix of the input GRM files (*.grm.bin, *.grm.id)"
    )
    parser.add_argument(
        "--par-visit", required=True,
        help="Path to par_visit_data.tsv"
    )
    parser.add_argument(
        "--identifiers",
        default="../HBCD_genomics_release_br_20p2/release_identifiers_20251211.csv",
        help="Release identifiers CSV (default: %(default)s)"
    )
    parser.add_argument(
        "--out-prefix", required=True,
        help="Prefix for the filtered output GRM files"
    )
    parser.add_argument(
        "--grm-id-col", default="FID",
        choices=("FID", "IID"),
        help="Which column in .grm.id to match against release_candid "
             "(default: '%(default)s')"
    )

    args = parser.parse_args()

    # ---- 1. read GRM -------------------------------------------------------
    print("Reading GRM …")
    grm_ids = read_grm_ids(args.grm_prefix)
    GRM = read_grm_bin(args.grm_prefix, grm_ids)

    # ---- 2. collect valid release_candid values from par_visit -------------
    print("\nCollecting valid subjects from par_visit …")
    valid = collect_valid_ids(args.par_visit, args.identifiers)
    print(f"  {len(valid)} unique subjects in par_visit")

    # ---- 3. find overlap ---------------------------------------------------
    grm_ids["release_candid"] = pd.to_numeric(
        grm_ids[args.grm_id_col], errors="coerce"
    ).astype("Int64")

    keep = grm_ids["release_candid"].isin(valid)
    n_keep = keep.sum()
    print(f"\n  GRM subjects in par_visit: {n_keep} / {len(grm_ids)}")

    if n_keep == 0:
        print("ERROR: no subjects overlap — check ID formats.", file=sys.stderr)
        sys.exit(1)

    # ---- 4. subset ---------------------------------------------------------
    grm_ids_sub = grm_ids.loc[keep, ["FID", "IID"]].reset_index(drop=True)
    idx = np.where(keep.values)[0]
    GRM_sub = GRM[np.ix_(idx, idx)]

    # ---- 5. write ----------------------------------------------------------
    print(f"\nWriting filtered GRM ({n_keep} × {n_keep}) …")
    write_grm(grm_ids_sub, GRM_sub, args.out_prefix)
    print("\nDone.")


if __name__ == "__main__":
    main()
