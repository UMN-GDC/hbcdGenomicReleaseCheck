#!/usr/bin/env python3
"""Filter batch and genomics data.

Python translation of the R data-processing sections in 01-setupNfilter.qmd.

Reads raw PLINK .fam, batch info, identifiers, and an exclusion list, then
right-joins against the par_visit table so that only subjects present in
par_visit survive.  Non-matching par_visit subjects get placeholder IDs
(NA1, NA2, …) and are written into the removal file.
"""
import pandas as pd
import numpy as np
from pathlib import Path


def main():
    # -- configuration --
    data_prefix = Path("/scratch.global/hbcd/full/full.QC8")
    release = "br_20p2"
    release_dir = Path(f"../HBCD_genomics_release_{release}/data/")
    release_dir.mkdir(parents=True, exist_ok=True)
    release_base = release_dir.parent  # ../HBCD_genomics_release_br_20p2/

    # -- batch info --
    batch = pd.read_csv("/projects/standard/basu_hbcd/shared/batch.info", sep=r"\s+")
    batch["relationship"] = batch["IID"].str[-1]
    batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
    batch = batch.drop(columns=["IID"])

    # -- exclusion list --
    exc = pd.read_excel(
        release_base / "HBCD_genetics_QC1_missing_race_LORIS.xlsx",
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

    # -- identifiers --
    identifiers = pd.read_csv(release_base / "release_identifiers_20251211.csv")
    identifiers["release_candid"] = pd.to_numeric(identifiers["release_candid"])
    identifiers["candid"] = pd.to_numeric(identifiers["candid"])
    identifiers = identifiers.dropna(subset=["release_candid"])

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
    par_visit = pd.read_csv(
        release_base / "par_visit_data_br20.2.tsv", sep="\t"
    )
    par_visit = (
        par_visit[["participant_id"]]
        .drop_duplicates()
        .dropna()
        .rename(columns={"participant_id": "candid"})
    )

    # -- combine / filter --
    combined = fam.merge(identifiers, how="left")
    combined = combined.merge(batch, how="left")
    combined = combined.merge(par_visit, how="right", on="candid")

    # de-identified FID / IID
    combined["FID"] = combined["release_candid"]
    combined["IID"] = (
        combined["FID"].astype(str) + combined["relationship"].fillna("")
    )
    combined["FID"] = combined["FID"].fillna(0).astype(int)

    # Subjects in par_visit but NOT in the genomics data get placeholder
    # IIDs (NA1, NA2, …) so plink2 --remove can drop them later.
    na_mask = combined["release_candid"].isna()
    combined.loc[na_mask, "IID"] = [
        f"NA{i+1}" for i in range(na_mask.sum())
    ]

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
