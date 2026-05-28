"""Validate the de-identified HBCD genomics release data.

Python / pytest translation of tests/testthat/test-ids.R.
"""
import pandas as pd
import numpy as np
from pathlib import Path

from _lib import load_par_visit_candids, load_identifiers, load_excluded_release_candids

RELEASE = "br_21p2"
RELEASE_DIR = (
    Path(__file__).resolve().parent.parent
    / f"HBCD_genomics_release_{RELEASE}"
    / "data"
)


def _load_data():
    fam = pd.read_csv(
        RELEASE_DIR / "hbcd.fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    batch = pd.read_csv(RELEASE_DIR / "batch.info", sep="\t")
    excluded = pd.read_csv(
        RELEASE_DIR.parent / "Removed_individuals.txt",
        delim_whitespace=True,
        header=None,
        names=["IID"],
    )
    return fam, batch, excluded


def _prepare_iid_test(fam, batch, excluded):
    iid_test = (
        fam.merge(batch, on="IID", how="outer")
        .merge(excluded, on="IID", how="left", suffixes=("", "_exc"))
        .assign(
            IID2=lambda df: pd.to_numeric(df["IID"].str[:10], errors="coerce"),
            Relation=lambda df: df["IID"].str[10],
        )
    )
    iid_test = iid_test[~iid_test["IID"].isin(excluded["IID"])].copy()
    return iid_test


def _id_lengths(iid_test):
    fid_len = (
        iid_test["FID"]
        .astype(str)
        .str.len()
        .where(iid_test["FID"].notna(), other=np.nan)
    )
    iid_len = iid_test["IID"].astype(str).str.len()
    lengths = (
        pd.DataFrame({"FID_len": fid_len, "IID_len": iid_len})
        .groupby(["FID_len", "IID_len"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    return lengths


# ── formatting tests ──────────────────────────────────────────────────────

def test_fam_and_batch_order_matches():
    fam, batch, _ = _load_data()
    assert (fam["IID"] == batch["IID"]).mean() == 1.0


def test_deid_fid_length_is_10_or_na():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    lengths = _id_lengths(iid_test)
    ok = lengths["FID_len"].isna() | (lengths["FID_len"] == 10.0)
    assert ok.mean() == 1.0


def test_deid_iid_length_is_11():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    lengths = _id_lengths(iid_test)
    assert (lengths["IID_len"] == 11.0).mean() == 1.0


def test_deid_fid_is_numeric():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    assert pd.api.types.is_numeric_dtype(iid_test["FID"])


def test_iid_first_10_digits_equal_fid():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    assert (iid_test["IID2"] == iid_test["FID"]).mean() == 1.0


def test_relation_is_C_or_M():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    assert ((iid_test["Relation"] == "C") | (iid_test["Relation"] == "M")).mean() == 1.0


# ── par_visit + exclusion constraint tests ────────────────────────────────

def _load_temp_fam():
    return pd.read_csv(
        RELEASE_DIR.parent / "temp.fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )


def test_all_output_iids_in_par_visit():
    """Every output IID maps (via identifiers) to a par_visit candid, and is
    NOT in the exclusion list."""
    fam, _, _ = _load_data()
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    fids = set(fam["FID"].unique()) - {0}

    # every FID must appear in identifiers
    ids_sub = identifiers[identifiers["release_candid"].isin(fids)]
    missing_fids = fids - set(ids_sub["release_candid"].astype(int))
    assert len(missing_fids) == 0, (
        f"{len(missing_fids)} FID(s) missing from identifiers"
    )

    # every FID's candid must be in par_visit
    candids_out = set(ids_sub["candid"].dropna().astype(int))
    orphans = candids_out - par_candids
    assert len(orphans) == 0, (
        f"{len(orphans)} output subject(s) NOT in par_visit"
    )

    # no excluded subject appears in output
    excluded_in_output = fids & exc_rc
    assert len(excluded_in_output) == 0, (
        f"{len(excluded_in_output)} excluded subject(s) found in output"
    )


def test_pipeline_no_subjects_lost_or_gained():
    """set(temp.fam IIDs) == set(hbcd.fam IIDs) ∪ set(Removed_individuals IIDs)."""
    fam, _, excluded = _load_data()
    temp_fam = _load_temp_fam()

    output = set(fam["IID"])
    removed = set(excluded["IID"])
    temp = set(temp_fam["IID"])
    union = output | removed

    extra = union - temp
    assert len(extra) == 0, f"IIDs in output/removed NOT in temp.fam: {sorted(extra)[:10]}"
    missing = temp - union
    assert len(missing) == 0, f"IIDs in temp.fam missing from both: {sorted(missing)[:10]}"
