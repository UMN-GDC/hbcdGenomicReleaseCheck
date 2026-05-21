"""Validate the de-identified HBCD genomics release data.

Python / pytest translation of tests/testthat/test-ids.R.
"""
import pandas as pd
import numpy as np
from pathlib import Path

RELEASE = "br_20p2"
# Test file lives under tests/ ; the release lives at ../HBCD_genomics_release_<RELEASE>/data/
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
    # Drop rows whose IID appears in the exclusion list
    iid_test = iid_test[~iid_test["IID"].isin(excluded["IID"])].copy()
    return iid_test


def _id_lengths(iid_test):
    # replicate R nchar(): numeric FID -> string length; NaN stays NaN
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


def test_fam_and_batch_order_matches():
    fam, batch, _ = _load_data()
    assert (fam["IID"] == batch["IID"]).mean() == 1.0


def test_deid_fid_length_is_10_or_na():
    """FID should be 10 characters when present, or NA for excluded subjects."""
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


# ---------------------------------------------------------------------------
# par_visit constraint tests
# ---------------------------------------------------------------------------

def _load_par_visit_and_ids():
    """Load par_visit and identifiers for cross-checking."""
    par_visit = pd.read_csv(
        RELEASE_DIR.parent / "par_visit_data_br20.2.tsv", sep="\t"
    )
    par_candid = set(
        par_visit["participant_id"].dropna().astype(int).unique()
    )

    identifiers = pd.read_csv(
        RELEASE_DIR.parent / "release_identifiers_20251211.csv"
    )
    identifiers["release_candid"] = pd.to_numeric(
        identifiers["release_candid"], errors="coerce"
    ).astype("Int64")
    identifiers["candid"] = pd.to_numeric(
        identifiers["candid"], errors="coerce"
    ).astype("Int64")
    identifiers = identifiers.dropna(subset=["release_candid", "candid"])

    # Build mapping: release_candid → candid
    candid_map = (
        identifiers[["release_candid", "candid"]]
        .drop_duplicates(subset="release_candid")
        .set_index("release_candid")["candid"]
        .to_dict()
    )
    return par_candid, candid_map


def _load_temp_fam():
    """Load the intermediate temp.fam produced by 01-filter_and_prepare_data."""
    return pd.read_csv(
        RELEASE_DIR.parent / "temp.fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )


def test_all_output_iids_in_par_visit():
    """Every IID in the final output maps to a candid that exists in par_visit.

    This is the *only-if* direction of the iff constraint.
    """
    fam, _, _ = _load_data()
    par_candid, candid_map = _load_par_visit_and_ids()

    # FID = release_candid (0 for placeholder subjects that were removed)
    fids_in_output = set(fam["FID"].unique()) - {0}

    missing_from_ids = fids_in_output - set(candid_map.keys())
    assert len(missing_from_ids) == 0, (
        f"{len(missing_from_ids)} FID(s) in hbcd.fam have no entry in "
        f"identifiers: {sorted(missing_from_ids)[:10]}"
    )

    candids_in_output = {candid_map[fid] for fid in fids_in_output}
    orphans = candids_in_output - par_candid
    assert len(orphans) == 0, (
        f"{len(orphans)} subject(s) in hbcd.fam have a candid that is NOT "
        f"in par_visit: {sorted(orphans)[:10]}"
    )


def test_pipeline_no_subjects_lost_or_gained():
    """The set of IIDs in temp.fam equals hbcd.fam ∪ Removed_individuals.

    This is the *if* direction of the iff constraint — every subject that
    entered through the right_join is accounted for in either final output
    or removal file.
    """
    fam, _, excluded = _load_data()
    temp_fam = _load_temp_fam()

    output_iids = set(fam["IID"])
    removed_iids = set(excluded["IID"])
    temp_iids = set(temp_fam["IID"])

    union = output_iids | removed_iids

    extra = union - temp_iids
    assert len(extra) == 0, (
        f"{len(extra)} IID(s) appear in output or removed but NOT in "
        f"temp.fam: {sorted(extra)[:10]}"
    )

    missing = temp_iids - union
    assert len(missing) == 0, (
        f"{len(missing)} IID(s) from temp.fam are missing from both "
        f"hbcd.fam and Removed_individuals: {sorted(missing)[:10]}"
    )
