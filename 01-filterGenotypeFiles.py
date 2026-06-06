#!/usr/bin/env python3

# source /projects/standard/basu_hbcd/shared/.venv/bin/activate

"""Filter batch and genomics data.

Rewrites the entire PLINK .fam with de-identified IDs (ALL original
subjects), then writes a keep_list for plink2 --keep so only subjects
in par_visit *and* not on the exclusion list survive in the final .bed.
"""
import pandas as pd
import numpy as np
from pathlib import Path

from _lib import (
    DATA_DIR,
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
    load_excluded_with_relationship,
    load_additional_excluded_pscids,
)


def main():
    data_prefix = DATA_DIR / "HBCD"
    release_dir = Path(
        "/projects/standard/basu_hbcd/shared/HBCD_genomics_release_br_21p2/data/"
    )
    release_dir.mkdir(parents=True, exist_ok=True)
    release_base = release_dir.parent

    # -- identifiers --
    identifiers = load_identifiers()

    # -- batch info --
    batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
    batch["relationship"] = batch["IID"].str[-1]
    batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
    batch = batch.drop(columns=["IID"])

    # -- original PLINK .fam (ALL subjects) --
    fam = pd.read_csv(
        str(data_prefix) + ".fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    fam["_idx"] = range(len(fam))
    fam["pscid"] = fam["IID"].str[-10:-1]
    fam["_orig_rel"] = fam["IID"].str[-1]
    fam["PHENO"] = "NONE"

    # -- merge to get de-identified IDs --
    combined = fam.merge(identifiers, how="left", on="pscid")
    combined = combined.merge(
        batch, how="left", on=["release_candid"]
    )
    # Drop cross-product rows — when a release_candid has batch data for
    # multiple relationships (e.g. C + M), keep only the row matching the
    # subject's original relationship.  Subjects without any batch match
    # survive with NaN and are dropped later by dropna.
    rel_ok = combined["relationship"].isna() | (
        combined["relationship"] == combined["_orig_rel"]
    )
    combined = combined[rel_ok]
    # Safety: one row per original .fam subject
    combined = combined.drop_duplicates(subset="_idx")

    # -- remove subjects listed in HBCDexclusions.csv --
    excluded_pscids = load_additional_excluded_pscids()
    n_before = len(combined)
    combined = combined[~combined["pscid"].isin(excluded_pscids)]
    n_removed = n_before - len(combined)
    if n_removed:
        print(f"  Removed {n_removed} subject(s) via HBCDexclusions.csv")

    # -- de-identified FID / IID for ALL subjects --
    combined["new_FID"] = combined["release_candid"].fillna(0).astype(int)
    combined["new_rel"] = combined["relationship"].fillna(combined["_orig_rel"])
    combined["new_IID"] = combined["new_FID"].astype(str) + combined["new_rel"]

    # restore original .fam row order
    combined = combined.sort_values("_idx")

    # -- write temp.fam with ALL subjects (de-identified, same order) --
    combined[["new_FID", "new_IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv(
        release_base / "temp.fam",
        sep=" ",
        index=False,
        header=False,
        na_rep="NA",
    )

    # -- single inclusive filter: valid = par_visit \\ excluded --
    par_candids = load_par_visit_candids()
    exc_release_candids = load_excluded_release_candids()

    valid_release_candids = (
        set(
            int(v)
            for v in identifiers[identifiers["release_candid"].isin(par_candids)][
                "release_candid"
            ]
            .dropna()
            .unique()
        )
        - exc_release_candids
    )

    combined["_valid"] = combined["release_candid"].isin(valid_release_candids)
    # require valid subjects to also have batch metadata so the output
    # .fam and batch.info contain exactly the same subjects
    valid = combined[combined["_valid"]].dropna(subset=["visit", "plate_number"])
    # keep the original .fam row order so batch.info matches hbcd.fam exactly
    valid = valid.sort_values("_idx")
    # safety net: one row per IID (catches any remaining merge artefacts)
    valid = valid.drop_duplicates(subset="new_IID")

    # -- write keep_list.txt for plink2 --keep --
    valid[["new_FID", "new_IID"]].to_csv(
        release_base / "keep_list.txt",
        sep=" ",
        index=False,
        header=False,
    )

    # -- write batch.info (same subjects as keep_list / hbcd.fam) --
    valid[["new_IID", "visit", "plate_number"]].rename(
        columns={"new_IID": "IID"}
    ).to_csv(
        release_dir / "batch.info",
        sep="\t",
        index=False,
    )

    # -- write Removed_individuals.txt (excluded subjects, documentation) --
    exc = load_excluded_with_relationship()
    exc["IID"] = (
        exc["release_candid"].astype(int).astype(str)
        + exc["relationship"].fillna("")
    )
    exc[["IID"]].drop_duplicates().to_csv(
        release_base / "Removed_individuals.txt",
        sep=" ",
        index=False,
        header=False,
    )


if __name__ == "__main__":
    main()
