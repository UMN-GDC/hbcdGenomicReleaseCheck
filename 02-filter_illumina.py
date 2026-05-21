#!/usr/bin/env python3
"""Filter Illumina Log R Ratio and BAF text files to subjects in par_visit.

Expects tabular text files (typically TSV) exported from GenomeStudio or a
similar pipeline.  Each file must contain a sample/subject identifier column
(configurable via --id-column).  Rows whose ID is NOT found in par_visit
(via the identifiers mapping) are dropped.

Usage example
-------------
python 02-filter_illumina.py \\
    lrr_data.tsv baf_data.tsv \\
    --par-visit ../HBCD_genomics_release_br_20p2/par_visit_data_br20.2.tsv \\
    --identifiers ../HBCD_genomics_release_br_20p2/release_identifiers_20251211.csv \\
    --id-column "Sample ID" \\
    --out-dir filtered/
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def load_valid_ids(par_visit_path, identifiers_path=None):
    """Return a set of valid subject IDs from the par_visit table.

    When *identifiers_path* is given, all known ID types (pscid,
    release_candid, …) for each par_visit participant are collected so
    that the filtering works regardless of which ID format the Illumina
    file uses.
    """
    par_visit = pd.read_csv(par_visit_path, sep="\t")
    par_candid = par_visit["participant_id"].dropna().unique()

    if identifiers_path is None:
        return set(str(c) for c in par_candid)

    ids = pd.read_csv(identifiers_path)
    ids["candid"] = pd.to_numeric(ids["candid"], errors="coerce")
    par = ids[ids["candid"].isin(par_candid)]

    valid = set()
    for col in ids.columns:
        if col in ("candid",):
            continue
        vals = par[col].dropna()
        if pd.api.types.is_numeric_dtype(vals):
            vals = vals.astype(int).astype(str)
        else:
            vals = vals.astype(str)
        valid |= set(vals)

    return valid


def filter_file(in_path, out_path, valid_ids, id_column):
    """Read a TSV, keep rows whose *id_column* value is in *valid_ids*,
    write to *out_path*."""
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    if id_column not in df.columns:
        available = [c for c in df.columns]
        raise KeyError(
            f"Column '{id_column}' not found in {in_path.name}. "
            f"Available columns: {available}"
        )

    n_before = len(df)
    df = df[df[id_column].isin(valid_ids)]
    n_after = len(df)
    print(f"  {in_path.name}: {n_before} rows → {n_after} rows "
          f"({n_before - n_after} removed)")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Filter Illumina LRR/BAF text files by par_visit"
    )
    parser.add_argument(
        "input_files", nargs="+",
        help="Illumina text file(s) to filter (space-separated or glob pattern)"
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
        "--id-column", default="Sample ID",
        help="Name of the sample/subject ID column in the Illumina files "
             "(default: '%(default)s')"
    )
    parser.add_argument(
        "--out-dir", default=".",
        help="Output directory for filtered files (default: current dir)"
    )

    args = parser.parse_args()

    # Resolve input files (support glob patterns)
    input_files = []
    for pattern in args.input_files:
        matched = sorted(Path().glob(pattern))
        if not matched:
            print(f"Warning: no files matched '{pattern}'", file=sys.stderr)
        input_files.extend(matched)

    if not input_files:
        print("No input files to process.", file=sys.stderr)
        sys.exit(1)

    print("Loading valid IDs from par_visit …")
    valid_ids = load_valid_ids(args.par_visit, args.identifiers)
    print(f"  {len(valid_ids)} unique identifiers found")

    out_dir = Path(args.out_dir)
    print(f"\nFiltering {len(input_files)} file(s):")
    for fpath in input_files:
        out_path = out_dir / fpath.name
        filter_file(fpath, out_path, valid_ids, args.id_column)

    print("\nDone.")


if __name__ == "__main__":
    main()
