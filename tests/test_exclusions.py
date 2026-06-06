"""Ensure all pscids in HBCDexclusions.csv are absent from release subjects."""

from _lib import load_identifiers, load_additional_excluded_pscids


def test_no_additional_excluded_pscids_in_release():
    identifiers = load_identifiers()
    identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()
    release_pscids = set(identifiers["pscid"].unique())

    excluded = load_additional_excluded_pscids()
    assert len(excluded) > 0, "No pscids loaded from HBCDexclusions.csv"

    overlap = release_pscids & excluded
    assert len(overlap) == 0, (
        f"{len(overlap)} excluded pscid(s) found in release subjects: "
        f"{sorted(overlap)[:20]}"
    )
