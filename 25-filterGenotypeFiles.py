#!/usr/bin/env python3
"""
Step 25: Genotype De-identification & Filtering for HBCD Release

This script performs the core de-identification and filtering logic for the HBCD
genomic release pipeline. It builds the release subject whitelist (keep_list.txt)
and creates a de-identified temp.fam file that preserves the original PLINK .bed
row count for compatibility with PLINK2 --fam.

Pipeline Phase: C (Release Filtering & Validation)
Dependencies: Phase A (QC & derivatives) + Phase B (de-identification) complete

Inputs (via _lib.py):
    - $HBCD_DATA_DIR/onlyQc.fam: De-identified PLINK .fam from Phase B
    - $HBCD_DATA_DIR/batch.info: Batch metadata (IID={pscid}{C|M})
    - $HBCD_DATA_DIR/release_identifiers_20260628.csv: PSCID -> release_candid crosswalk
    - $HBCD_DATA_DIR/HBCDexclusions.csv: Additional PSCID exclusions (multi-column)
    - $HBCD_DATA_DIR/HBCD_genetics_QC1_missing_race_LORIS.xlsx: Excel exclusions (Exclude_Summary)
    - $HBCD_DATA_DIR/par_visit_data_br21_1.tsv: par_visit participation data

Outputs (to $RELEASE_BASE/ and $RELEASE_DIR/GDA/):
    - temp.fam: ALL QC subjects with de-identified FID/IID (preserves .bed row count)
    - keep_list.txt: Release whitelist (FID + IID for valid subjects)
    - GDA/batch.info: Batch metadata for release subjects (IID, visit, plate_number)
    - GDA/removed_individuals.txt: Excluded IIDs for documentation

Filter Logic:
    1. Load identifiers, apply HBCDexclusions.csv PSCID filter
    2. Map batch.info IIDs (pscid -> release_candid)
    3. Load onlyQc.fam (already de-identified by Phase B)
    4. Merge: fam x identifiers x batch_info
    5. Relationship filter: keep batch row matching original C/M suffix
    6. Two-stage inclusive filter:
       - Stage 1: par_visit subjects only (completed visits)
       - Stage 2: Remove Excel Exclude_Summary release_candids
    7. Require non-missing visit + plate_number metadata
    8. Deduplicate by IID (keep first)

Usage:
    python 25-filterGenotypeFiles.py

Environment Variables:
    HBCD_RELEASE: Release tag (default: br_21p3)
    HBCD_DATA_DIR: Base data directory (default: /projects/standard/basu_hbcd/shared/data)

Example:
    export HBCD_RELEASE=br31p2
    export HBCD_DATA_DIR=/projects/standard/basu_hbcd/shared/data
    python 25-filterGenotypeFiles.py
"""

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
gda_dir = release_dir / "GDA"
gda_dir.mkdir(parents=True, exist_ok=True)
release_base = get_release_base()

# -- identifiers --
identifiers = load_identifiers()
excluded_pscids = load_additional_excluded_pscids()
n_exc = identifiers["pscid"].isin(excluded_pscids).sum()
if n_exc:
    print(f"  Excluding {n_exc} subject(s) from identifiers via HBCDexclusions.csv")
    identifiers = identifiers[~identifiers["pscid"].isin(excluded_pscids)]
print(f"  Identifiers total              : {len(identifiers):>6}")

# -- batch info (IID = {pscid}{C|M}; map pscid → release_candid) --
batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
n_raw = len(batch)
valid_iid = batch["IID"].astype(str).str.match(r"^[A-Za-z]{5}\d{4}[CM]$")
n_dropped = n_raw - valid_iid.sum()
if n_dropped:
    print(f"  Dropped {n_dropped} batch row(s) with non-standard IID")
    batch = batch[valid_iid].copy()
batch["relationship"] = batch["IID"].str[-1]
batch["pscid"] = batch["IID"].str[:-1]
pscid_to_rc = (
    identifiers[["pscid", "release_candid"]]
    .drop_duplicates()
    .assign(pscid=lambda x: x["pscid"].astype(str))
    .set_index("pscid")["release_candid"]
)
batch["release_candid"] = batch["pscid"].map(pscid_to_rc)
n_unmapped = batch["release_candid"].isna().sum()
if n_unmapped:
    print(f"  Warning: {n_unmapped} batch pscid(s) unmapped to release_candid")
    batch = batch.dropna(subset=["release_candid"])
batch = batch.drop(columns=["IID", "pscid"])
print(f"  Batch info usable              : {len(batch):>6}")

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
    gda_dir / "batch.info",
    sep="\t",
    index=False,
)

# -- write GDA/removed_individuals.txt --
exc = load_excluded_with_relationship()
exc["IID"] = (
    exc["release_candid"].astype(int).astype(str)
    + exc["relationship"].fillna("")
)
exc[["IID"]].drop_duplicates().to_csv(
    gda_dir / "removed_individuals.txt",
    sep=" ",
    index=False,
    header=False,
)
print(f"  GDA/removed_individuals.txt    : {len(exc):>6}")
