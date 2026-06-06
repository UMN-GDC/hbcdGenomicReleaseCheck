"""Shared helpers for the HBCD genomics release pipeline."""
from pathlib import Path

import pandas as pd

DATA_DIR = Path("/projects/standard/basu_hbcd/shared/data")
RELEASE = "br_21p2"
PAR_VISIT_FILE = "par_visit_data_br21_1.tsv"
IDENTIFIERS_FILE = "release_identifiers_20260526.csv"
EXCLUSION_FILE = "HBCD_genetics_QC1_missing_race_LORIS.xlsx"


def load_par_visit_candids(path=None):
    """Return a set of integral participant IDs (*candid*) from par_visit.

    Filters to completed visits (``par_visit_data_visit_missed == 'No'``)
    and strips the ``sub-`` prefix from ``participant_id`` per
    ``int(val[4:])``.
    """
    if path is None:
        path = DATA_DIR / PAR_VISIT_FILE
    pv = pd.read_csv(path, sep="\t").query("par_visit_data_visit_missed == 'No'")
    def _to_int(v):
        s = str(v)
        return int(s[4:] if s.startswith("sub-") else s)
    ids = pv["participant_id"].dropna().apply(_to_int).unique()
    return set(int(x) for x in ids)


def load_identifiers(path=None):
    """Load the release-identifiers CSV, parse numeric columns, return DataFrame."""
    if path is None:
        path = DATA_DIR / IDENTIFIERS_FILE
    df = (
        pd.read_csv(path)
        .query("release_candid != 'release_candid'")  # drop placeholder rows
        .assign(release_candid=lambda x: pd.to_numeric(x["release_candid"]))
        .assign(candid=lambda x: pd.to_numeric(x["candid"]))
        .dropna(subset=["release_candid"])
        .drop_duplicates(subset=["pscid", "candid", "release_candid"])
    )
    return df


def load_excluded_release_candids(path=None):
    """Return a set of ``release_candid`` (int) from the Excel exclusion list."""
    if path is None:
        path = DATA_DIR / EXCLUSION_FILE
    exc = pd.read_excel(path, sheet_name="Exclude_Summary")
    exc[["_", "release_candid", "_pscid"]] = exc["Sampled ID"].str.split(
        "_", n=2, expand=True
    )
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    exc = exc.dropna(subset=["release_candid"])
    return set(int(x) for x in exc["release_candid"].unique())


def load_excluded_with_relationship(path=None):
    """Return a DataFrame with ``release_candid`` and ``relationship`` from
    the Excel exclusion list, for constructing IIDs of removed individuals."""
    if path is None:
        path = DATA_DIR / EXCLUSION_FILE
    exc = pd.read_excel(path, sheet_name="Exclude_Summary")
    exc[["_", "release_candid", "_pscid"]] = exc["Sampled ID"].str.split(
        "_", n=2, expand=True
    )
    exc["release_candid"] = pd.to_numeric(exc["release_candid"])
    exc["relationship"] = exc["Study_ID"].str[-1]
    return exc.dropna(subset=["release_candid"])[["release_candid", "relationship"]]


ADDITIONAL_EXCLUSIONS_FILE = "HBCD_exlcustions.xlsx"


def parse_additional_exclusion_lists(path=None):
    """Read the first sheet of HBCD_exlcustions.xlsx, skip first 13 rows,
    then parse stacked lists of pscids stacked in column A.

    The column contains alternating patterns::

        [empty lines]
        Header text describing the list
        pscid-1
        pscid-2
        ...
        [empty lines]
        Next header
        pscid-a
        ...

    Returns a dict of ``{header: set_of_pscid_strings}``.
    """
    if path is None:
        path = DATA_DIR / ADDITIONAL_EXCLUSIONS_FILE
    raw = pd.read_excel(path, sheet_name=0, skiprows=13, header=None)
    col = raw.iloc[:, 0].dropna().astype(str).str.strip()

    lists = {}
    current_header = None
    current_items = []

    for val in col:
        if not val or val == "nan":
            continue
        is_pscid = val.isdigit() and len(val) <= 12
        if is_pscid:
            if current_header is not None:
                current_items.append(val)
        else:
            if current_header is not None and current_items:
                lists[current_header] = set(current_items)
            current_header = val
            current_items = []

    if current_header is not None and current_items:
        lists[current_header] = set(current_items)

    return lists
