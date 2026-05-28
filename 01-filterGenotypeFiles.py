#!/usr/bin/env python3

# source /projects/standard/basu_hbcd/shared/.venv/bin/activate


"""Filter batch and genomics data.

Input:
    - pre deIdentified Plink .fam 
    - pre deIdentified demographics
    - deIdenfified batch info 
    - release deIdentification mapping
    - par_visit table (not in use yet)
Output:
    - Subjects passing QC filter


Reads deIdentified PLINK .fam, batch info, identifiers, and an exclusion list, then (though it is not implemented yet)
right-joins against the par_visit table so that only subjects present in
par_visit survive.  Non-matching par_visit subjects get placeholder IDs
(NA1, NA2, …) and are written into the removal file.
"""
import pandas as pd
import numpy as np
from pathlib import Path


def main():
    # -- configuration --
    projectDir = Path("/projects/standard/basu_hbcd/shared")
    dataDir = projectDir / "data"
    data_prefix = dataDir / "HBCD"
    release = "br_21p2"
    parVisit = "par_visit_data_br21_1.tsv"
    release_dir = Path(f"/projects/standard/basu_hbcd/shared/HBCD_genomics_release_{release}/data/")
    release_dir.mkdir(parents=True, exist_ok=True)
    release_base = release_dir.parent  # ../HBCD_genomics_release_br_20p2/

    # -- Genotype array batch info --
    batch = pd.read_csv(dataDir / "batch.info", sep=r"\s+")
    batch["relationship"] = batch["IID"].str[-1]
    batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
    batch = batch.drop(columns=["IID"])

    # -- exclusion list --
    exc = pd.read_excel(
        dataDir / "HBCD_genetics_QC1_missing_race_LORIS.xlsx",
        sheet_name="Exclude_Summary",
    )
    exc[["something", "release_candid", "pscidR"]] = exc["Sampled ID"].str.split(
        "_", n=2, expand=True
    )
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    exc["relationship"] = exc["Study_ID"].str[-1]
    exc["pscid"] = exc["Study_ID"].str[:-1]
    exc = (
        exc.drop(columns=["Study_ID", "something", "pscidR"])
        .dropna(subset=["release_candid"])
    )

    # -- de-identifiers --
    identifiers = (
            pd.read_csv(dataDir / "release_identifiers_20260526.csv")
            .query("release_candid != 'release_candid'") # filter some placeholders
            .assign(release_candid = lambda x: pd.to_numeric(x["release_candid"]))
            .assign(candid = lambda x: pd.to_numeric(x["candid"]))
            .dropna(subset=["release_candid"])
            .drop_duplicates(subset=['pscid', 'candid', 'release_candid'])
    )

    # -- raw PLINK .fam --
    fam = pd.read_csv(
        str(data_prefix) + ".fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    fam["pscid"] = fam["IID"].str[-10:-1]   # str_sub(IID, -10, -2)
    fam["relationship"] = fam["IID"].str[-1]  # str_sub(IID, -1, -1)
    fam["PHENO"] = "NONE"
    fam = fam.drop(columns=["FID", "IID"])

    # -- parent-visit table (used as a filter) --
    par_visit = (
            pd.read_csv(
            dataDir / "par_visit_data_br21_1.tsv", sep="\t"
        )
            .query("par_visit_data_visit_missed == 'No'")
            .rename(columns={"participant_id": "candid"})
            .assign(candid=lambda x: x["candid"].apply(
                lambda val: int(val[4:])))
    )[["candid"]].drop_duplicates().dropna()

    # -- filter --
    combined = fam.merge(identifiers, how="left", on = "pscid")
    combined = combined.merge(batch, how="left", on = ["release_candid", "relationship"])
    combined = combined[combined['release_candid'].isin(par_visit["candid"])]

    # de-identified FID / IID
    combined["FID"] = combined["release_candid"]
    combined["IID"] = (
        combined["FID"].astype(str) + combined["relationship"].fillna("")
    )
    combined["FID"] = combined["FID"].fillna(0).astype(int)

    # Exclude
    combined = combined[~combined["FID"].isin(exc.release_candid)].dropna(subset = ["FID", "IID"])

    # Subjects in par_visit but NOT in the genomics data get placeholder
    # IIDs (NA1, NA2, …) so plink2 --remove can drop them later.
    combined.loc[na_mask, "IID"] = [
        f"NA{i+1}" for i in range(na_mask.sum())
    ]

    sum(~combined["FID"].isin(par_visit["candid"]))

    # -- write temp.fam (space-delimited, no header) --
    combined[["FID", "IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv(
        release_base / "temp.fam",
        sep=" ",
        index=False,
        header=False,
        na_rep="NA",
    )

    # -- write batch.info (tab-delimited, excludes NA placeholders) --
    valid = combined[~na_mask]
    valid[["IID", "visit", "plate_number"]].to_csv(
        release_dir / "batch.info",
        sep="\t",
        index=False,
    )

    # -- write exclusion file (Removed_individuals.txt) --
    exc_out = combined.merge(
        exc, on=["release_candid", "relationship"], how="outer"
    )
    exc_out["FID"] = exc_out["release_candid"]
    exc_out["IID"] = (
        exc_out["FID"]
        .fillna(0)
        .astype(int)
        .astype(str)
        + exc_out["relationship"].fillna("")
    )
    exc_out = exc_out.loc[
        exc_out["Exclusion Reason"].notna() | (exc_out["FID"].fillna(0) == 0)
    ]
    exc_out = exc_out["IID"].drop_duplicates().to_frame()
    exc_out.to_csv(
        release_base / "Removed_individuals.txt",
        sep=" ",
        index=False,
        header=False,
    )


if __name__ == "__main__":
    main()
