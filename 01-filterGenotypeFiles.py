#!/usr/bin/env python3

# source /projects/standard/basu_hbcd/shared/.venv/bin/activate

"""Filter batch and genomics data.

Reads de-identified PLINK .fam, batch info, identifiers, and an exclusion
list, keeps only subjects that appear in par_visit *and* are not on the
exclusion list, then writes filtered output files.

The par_visit / exclusion filter is applied as a single inclusive set:
    valid_release_candids = par_visit_candids → identifiers → release_candid
                           \\ excluded_release_candids
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
)


def main():
    # -- configuration --
    data_prefix = DATA_DIR / "HBCD"
    release = "br_21p2"
    release_dir = Path(
        f"/projects/standard/basu_hbcd/shared/HBCD_genomics_release_{release}/data/"
    )
    release_dir.mkdir(parents=True, exist_ok=True)
    release_base = release_dir.parent

    # -- Genotype array batch info --
    batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
    batch["relationship"] = batch["IID"].str[-1]
    batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
    batch = batch.drop(columns=["IID"])

    # -- identifiers --
    identifiers = load_identifiers()

    # -- raw PLINK .fam --
    fam = pd.read_csv(
        str(data_prefix) + ".fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    fam["pscid"] = fam["IID"].str[-10:-1]
    fam["relationship"] = fam["IID"].str[-1]
    fam["PHENO"] = "NONE"
    fam = fam.drop(columns=["FID", "IID"])

    # -- combine data --
    combined = fam.merge(identifiers, how="left", on="pscid")
    combined = combined.merge(batch, how="left", on=["release_candid", "relationship"])

    # -- single inclusive filter: valid_release_candid = par_visit \\ excluded --
    par_candids = load_par_visit_candids()
    exc_release_candids = load_excluded_release_candids()

    valid_release_candids = (
        set(
            identifiers[identifiers["candid"].isin(par_candids)][
                "release_candid"
            ]
            .astype(int)
            .unique()
        )
        - exc_release_candids
    )
    combined = combined[combined["release_candid"].isin(valid_release_candids)]

    # -- de-identified FID / IID --
    combined["FID"] = combined["release_candid"]
    combined["IID"] = combined["FID"].astype(str) + combined["relationship"].fillna("")
    combined["FID"] = combined["FID"].fillna(0).astype(int)

    # -- write temp.fam (space-delimited, no header) --
    combined[["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv(
        release_base / "temp.fam",
        sep=" ",
        index=False,
        header=False,
        na_rep="NA",
    )

    # -- write batch.info (tab-delimited) --
    combined[["IID", "visit", "plate_number"]].to_csv(
        release_dir / "batch.info",
        sep="\t",
        index=False,
    )

    # -- write exclusion file (Removed_individuals.txt) --
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
