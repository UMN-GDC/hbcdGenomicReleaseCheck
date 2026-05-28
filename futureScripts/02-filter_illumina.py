#!/usr/bin/env python3
"""Filter Illumina Log R Ratio and BAF text files to subjects in par_visit
minus the Excel exclusion list.

Applies a single inclusive filter:
    valid = par_visit_candids → identifiers → {all ID formats}
            \\ {excluded_release_candids → identifiers → {all ID formats}}

Usage
-----
python futureScripts/02-filter_illumina.py \\
    lrr_data.tsv baf_data.tsv \\
    --out-dir filtered/
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _lib import (
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
    DATA_DIR,
)


def _collect_all_id_formats(release_candids, identifiers):
    """Return every known ID format for the given set of release_candids.

    Collects values from *every* column in identifiers except ``candid``,
    so the returned set works regardless of which ID format the Illumina
    file uses (pscid, release_candid as string, etc.).
    """
    sub = identifiers[identifiers["release_candid"].isin(release_candids)]
    ids = set()
    for col in identifiers.columns:
        if col == "candid":
            continue
        vals = sub[col].dropna()
        if pd.api.types.is_numeric_dtype(vals):
            vals = vals.astype(int).astype(str)
        else:
            vals = vals.astype(str)
        ids |= set(vals)
    return ids


def filter_file(in_path, out_path, valid_ids, id_column):
    """Read a TSV, keep rows whose *id_column* is in *valid_ids*; write."""
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    if id_column not in df.columns:
        raise KeyError(
            f"Column '{id_column}' not found in {in_path.name}. "
            f"Available: {df.columns.tolist()}"
        )
    n_before = len(df)
    df = df[df[id_column].isin(valid_ids)]
    print(f"  {in_path.name}: {n_before} → {len(df)} rows")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Filter Illumina LRR/BAF files by par_visit minus exclusion"
    )
    parser.add_argument(
        "input_files", nargs="+",
        help="Illumina text file(s) (space-separated or glob pattern)"
    )
    parser.add_argument(
        "--id-column", default="Sample ID",
        help="Sample ID column name in Illumina files (default: '%(default)s')"
    )
    parser.add_argument(
        "--out-dir", default=".",
        help="Output directory (default: current dir)"
    )

    args = parser.parse_args()

    # Resolve glob patterns
    files = []
    for p in args.input_files:
        matched = sorted(Path().glob(p))
        if not matched:
            print(f"Warning: no files matched '{p}'", file=sys.stderr)
        files.extend(matched)
    if not files:
        print("No input files.", file=sys.stderr)
        sys.exit(1)

    # Single inclusive filter set
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    par_rc = set(
        identifiers[identifiers["candid"].isin(par_candids)]["release_candid"]
        .dropna().astype(int).unique()
    )
    valid_rc = par_rc - exc_rc
    valid_ids = _collect_all_id_formats(valid_rc, identifiers)
    print(f"Valid IDs (par_visit \\ excluded): {len(valid_ids)}")

    out_dir = Path(args.out_dir)
    for fpath in files:
        filter_file(fpath, out_dir / fpath.name, valid_ids, args.id_column)
    print("Done.")


if __name__ == "__main__":
    main()
