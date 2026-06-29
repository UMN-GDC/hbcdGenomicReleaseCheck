#!/usr/bin/env python3

# step 3: filter batch and genomics data for the HBCD release
# source /projects/standard/basu_hbcd/shared/.venv/bin/activate

import pandas as pd

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

data_prefix = DATA_DIR / "onlyQc"
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
print(f"  Identifiers total              : {len(identifiers):>6}")

# -- batch info --
batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
batch["relationship"] = batch["IID"].str[-1]
batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
batch = batch.drop(columns=["IID"])
print(f"  Batch info total               : {len(batch):>6}")

# -- PLINK .fam (from onlyQc) --
fam = pd.read_csv(
    str(data_prefix) + ".fam",
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
fam["_idx"] = range(len(fam))
fam["PHENO"] = "NONE"

# IID = {release_candid}{C|M}; FID = {release_candid}
fam["release_candid"] = pd.to_numeric(fam["IID"].astype(str).str[:-1], errors="coerce")
fam["_orig_rel"] = fam["IID"].str[-1]
print(f"  onlyQc.fam total               : {len(fam):>6}")

# -- merge to get de-identified IDs --
combined = fam.merge(identifiers, how="left", on="release_candid")
n_no_id = combined["pscid"].isna().sum()
print(f"  Missing identifiers merge      : {n_no_id:>6}")

combined = combined.merge(batch, how="left", on="release_candid")
n_no_batch = combined["relationship"].isna().sum()
print(f"  Missing batch.info merge       : {n_no_batch:>6}")

rel_ok = combined["relationship"].isna() | (
    combined["relationship"] == combined["_orig_rel"]
)
n_rel_bad = (~rel_ok).sum()
combined = combined[rel_ok]
if n_rel_bad:
    print(f"  Relationship mismatch removed  : {n_rel_bad:>6}")

combined = combined.drop_duplicates(subset="_idx")

# -- de-identified FID / IID for ALL subjects --
combined["new_FID"] = pd.to_numeric(combined["release_candid"], errors="coerce").fillna(0).astype(int)
combined["new_rel"] = combined["relationship"].fillna(combined["_orig_rel"])
has_rc = combined["release_candid"].notna()
combined["new_IID"] = combined["new_FID"].astype(str) + combined["new_rel"]
combined.loc[~has_rc, "new_IID"] = (
    "0_" + combined.loc[~has_rc, "release_candid"].fillna("unknown").astype(str)
)

combined = combined.sort_values("_idx")

# -- write temp.fam (all subjects, preserves row count) --
combined[["new_FID", "new_IID", "PAT", "MAT", "SEX", "PHENO"]].to_csv(
    release_base / "temp.fam",
    sep=" ",
    index=False,
    header=False,
    na_rep="NA",
)
print(f"  Wrote temp.fam                 : {len(combined):>6}")

# -- single inclusive filter: valid = par_visit \ excluded --
par_candids = load_par_visit_candids()
exc_release_candids = load_excluded_release_candids()
print(f"  Par_visit candids              : {len(par_candids):>6}")
print(f"  Excel-excluded release_candids : {len(exc_release_candids):>6}")

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
print(f"  Valid release_candids          : {len(valid_release_candids):>6}")

n_not_par = (~combined["release_candid"].isin(valid_release_candids)).sum()
combined["_valid"] = combined["release_candid"].isin(valid_release_candids)
print(f"  Not in valid_release_candids   : {n_not_par:>6}")

valid = combined[combined["_valid"]]
n_no_meta = valid["visit"].isna().sum() + valid["plate_number"].isna().sum()
valid = valid.dropna(subset=["visit", "plate_number"])
if n_no_meta:
    print(f"  Missing visit/plate_number     : {n_no_meta:>6}")

valid = valid.sort_values("_idx")
n_dup = len(valid) - len(valid.drop_duplicates(subset="new_IID"))
valid = valid.drop_duplicates(subset="new_IID")
if n_dup:
    print(f"  Duplicate IIDs removed         : {n_dup:>6}")

print(f"  Release subjects (keep_list)   : {len(valid):>6}")

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
print(f"  Removed_individuals.txt        : {len(exc):>6}")
