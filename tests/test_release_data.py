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
    assert set(fam["IID"]) == set(batch["IID"]), "fam / batch IID sets differ"


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
    """Every output IID (fam *and* batch) maps to a par_visit candid and is
    NOT in the exclusion list."""
    fam, batch, _ = _load_data()
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    # hbcd.fam: use FID column (release_candid)
    fids = set(fam["FID"].unique()) - {0}
    ids_sub = identifiers[identifiers["release_candid"].isin(fids)]
    missing = fids - set(ids_sub["release_candid"].astype(int))
    assert len(missing) == 0, (
        f"{len(missing)} FID(s) in hbcd.fam missing from identifiers: "
        f"{sorted(missing)[:10]}"
    )
    candids_out = set(ids_sub["candid"].dropna().astype(int))
    orphans = candids_out - par_candids
    assert len(orphans) == 0, (
        f"{len(orphans)} subject(s) in hbcd.fam NOT in par_visit: "
        f"{sorted(orphans)[:10]}"
    )
    excluded_in = fids & exc_rc
    assert len(excluded_in) == 0, (
        f"{len(excluded_in)} excluded subject(s) found in hbcd.fam: "
        f"{sorted(excluded_in)[:10]}"
    )

    # batch.info: no FID column; extract release_candid from IID
    batch_rc = set(pd.to_numeric(batch["IID"].str[:-1], errors="coerce").dropna().astype(int))
    ids_sub = identifiers[identifiers["release_candid"].isin(batch_rc)]
    missing = batch_rc - set(ids_sub["release_candid"].astype(int))
    assert len(missing) == 0, (
        f"{len(missing)} subject(s) in batch.info missing from identifiers: "
        f"{sorted(missing)[:10]}"
    )
    candids_out = set(ids_sub["candid"].dropna().astype(int))
    orphans = candids_out - par_candids
    assert len(orphans) == 0, (
        f"{len(orphans)} subject(s) in batch.info NOT in par_visit: "
        f"{sorted(orphans)[:10]}"
    )
    excluded_in = batch_rc & exc_rc
    assert len(excluded_in) == 0, (
        f"{len(excluded_in)} excluded subject(s) found in batch.info: "
        f"{sorted(excluded_in)[:10]}"
    )


def test_filter_correctness():
    """Re-derive the expected subject set from par_visit \\ excluded and
    verify it exactly matches hbcd.fam — the if-and-only-if constraint."""
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    # subjects in par_visit, mapped to release_candid, minus excluded
    expected_rc = (
        set(
            int(v)
            for v in identifiers[identifiers["candid"].isin(par_candids)][
                "release_candid"
            ]
            .dropna()
            .unique()
        )
        - exc_rc
    )

    # restrict to release_candids that actually exist in the data
    temp = _load_temp_fam()
    rc_in_data = set(temp["FID"].unique()) - {0}
    expected_rc &= rc_in_data

    # actual FIDs in the output
    fam, _, _ = _load_data()
    actual_rc = set(fam["FID"].unique()) - {0}

    missing = expected_rc - actual_rc
    extra = actual_rc - expected_rc
    assert len(missing) == 0, (
        f"{len(missing)} subject(s) expected in output but missing: "
        f"{sorted(missing)[:10]}"
    )
    assert len(extra) == 0, (
        f"{len(extra)} subject(s) in output but not expected: "
        f"{sorted(extra)[:10]}"
    )
