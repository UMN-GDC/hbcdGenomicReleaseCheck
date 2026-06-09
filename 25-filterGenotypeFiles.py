#!/usr/bin/env python3

# step 3: filter batch and genomics data for the HBCD release
# source /projects/standard/basu_hbcd/shared/.venv/bin/activate

import pandas as pd
from pathlib import Path

from _lib import (
    DATA_DIR,
    get_release_dir,
    get_release_base,
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
    load_excluded_with_relationship,
    load_additional_excluded_pscids,
)

data_prefix = Path("onlyQc")
release_dir = get_release_dir()
release_dir.mkdir(parents=True, exist_ok=True)
release_base = get_release_base()

# -- identifiers --
identifiers = load_identifiers()
excluded_pscids = load_additional_excluded_pscids()
n_exc = identifiers["pscid"].isin(excluded_pscids).sum()
if n_exc:
    print(f"  Excluding {n_exc} subject(s) from identifiers via HBCDexclusions.csv")
    identifiers = identifiers[~identifiers["pscid"].isin(excluded_pscids)]

# -- batch info --
batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
batch["relationship"] = batch["IID"].str[-1]
batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
batch = batch.drop(columns=["IID"])

# -- PLINK .fam (from onlyQc) --
fam = pd.read_csv(
    str(data_prefix) + ".fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
fam["_idx"] = range(len(fam))
fam["pscid"] = fam["IID"].astype(str)
fam["_orig_rel"] = fam["IID"].str[-1]
fam["PHENO"] = "NONE"

# -- merge to get de-identified IDs --
combined = fam.merge(identifiers, how="left", on="pscid")
combined = combined.merge(batch, how="left", on=["release_candid"])
rel_ok = combined["relationship"].isna() | (
    combined["relationship"] == combined["_orig_rel"]
)
combined = combined[rel_ok]
combined = combined.drop_duplicates(subset="_idx")

# -- de-identified FID / IID for ALL subjects --
combined["new_FID"] = combined["release_candid"].fillna(0).astype(int)
combined["new_rel"] = combined["relationship"].fillna(combined["_orig_rel"])
has_rc = combined["release_candid"].notna()
combined["new_IID"] = combined["new_FID"].astype(str) + combined["new_rel"]
combined.loc[~has_rc, "new_IID"] = (
    "0_" + combined.loc[~has_rc, "pscid"] + "_" + combined.loc[~has_rc, "_orig_rel"]
)

combined = combined.sort_values("_idx")

# -- write temp.fam --
combined[["new_FID", "new_IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv(
    release_base / "temp.fam",
    sep=" ",
    index=False,
    header=False,
    na_rep="NA",
)

# -- single inclusive filter: valid = par_visit \ excluded --
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
valid = combined[combined["_valid"]].dropna(subset=["visit", "plate_number"])
valid = valid.sort_values("_idx")
valid = valid.drop_duplicates(subset="new_IID")

# -- write keep_list.txt --
valid[["new_FID", "new_IID"]].to_csv(
    release_base / "keep_list.txt",
    sep=" ",
    index=False,
    header=False,
)

# -- write batch.info --
valid[["new_IID", "visit", "plate_number"]].rename(
    columns={"new_IID": "IID"}
).to_csv(
    release_dir / "batch.info",
    sep="\t",
    index=False,
)

# -- write Removed_individuals.txt --
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
