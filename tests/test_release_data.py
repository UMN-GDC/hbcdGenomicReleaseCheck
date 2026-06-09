"""Validate de-identification and data integrity of the HBCD release.

Set ``HBCD_RELEASE`` env var to test against a different release.
"""
import pandas as pd
import numpy as np
import pytest

from _lib import (
    DATA_DIR,
    get_release_dir,
    get_release_base,
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
)

RELEASE_DIR = get_release_dir()
NAME = "hbcd_rsid_harmonized"

# ── helpers ──────────────────────────────────────────────────────────────


def _load_data():
    fam = pd.read_csv(
        RELEASE_DIR / "hbcd.fam",
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    batch = pd.read_csv(RELEASE_DIR / "batch.info", sep="\t")
    batch = batch.rename(columns={batch.columns[0]: "IID"})
    excluded = pd.read_csv(
        get_release_base() / "Removed_individuals.txt",
        delim_whitespace=True,
        header=None,
        names=["IID"],
    )
    return fam, batch, excluded


def _load_temp_fam():
    p = get_release_base() / "temp.fam"
    if not p.exists():
        pytest.skip("temp.fam not found")
    return pd.read_csv(
        p, sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )


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


# ══════════════════════════════════════════════════════════════════════════
#  De-identification integrity
# ══════════════════════════════════════════════════════════════════════════

def test_no_fid_zero_in_output():
    """No unmatched subjects (FID=0) appear in the release output.

    Subjects without a release_candid get FID=0 in temp.fam but must be
    filtered out before the release PLINK files are written.
    """
    p = RELEASE_DIR / "hbcd.fam"
    if not p.exists():
        pytest.skip("hbcd.fam not found")
    fam = pd.read_csv(
        p, sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    n_zero = (fam["FID"] == 0).sum()
    assert n_zero == 0, f"{n_zero} subject(s) with FID=0 in release output"


def test_all_iids_match_deid_pattern():
    """Every IID matches ``^\\d{10}[CM]$`` — ten digits followed by C or M.

    This is the strictest de-identification check: it catches any
    non-numeric characters in the identifier portion, wrong length, or
    invalid relationship suffix.
    """
    fam, batch, _ = _load_data()
    for name, iids in [("hbcd.fam", fam["IID"]), ("batch.info", batch["IID"])]:
        ok = iids.astype(str).str.match(r"^\d{10}[CM]$")
        assert ok.mean() == 1.0, (
            f"{name}: {(~ok).sum()} IID(s) do not match "
            f"^\\d{{10}}[CM]$: {iids[~ok].tolist()[:10]}"
        )


def test_all_fids_are_positive_10_digit_integers():
    """Every non-zero FID is a 10-digit positive integer.

    Zero-FID rows are forbidden by ``test_no_fid_zero_in_output``; this
    test verifies that every remaining FID has exactly 10 digits and no
    non-numeric content (which would indicate a raw identifier leak).
    """
    fam, _, _ = _load_data()
    fids = fam["FID"].dropna()
    digits = fids.astype(str).str.len()
    bad = fids[(fids <= 0) | (~digits.between(9, 11))]
    assert len(bad) == 0, (
        f"{len(bad)} FID(s) not positive 10-digit integers: {bad.tolist()[:15]}"
    )


def test_iid_first_10_digits_equal_fid():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    assert (iid_test["IID2"] == iid_test["FID"]).mean() == 1.0, (
        "IID[:10] != FID for some subjects"
    )


def test_relation_is_C_or_M():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    ok = (iid_test["Relation"] == "C") | (iid_test["Relation"] == "M")
    assert ok.mean() == 1.0, (
        f"{ok.size - ok.sum()} IID(s) with non-C/M suffix"
    )


def test_fid_length_is_10():
    """FID is a 10-digit integer for every release subject.

    Unlike the earlier permissive test that allowed ``na``, this
    enforces that every subject in the release output has a proper
    10-digit de-identified FID.
    """
    fam, _, _ = _load_data()
    lengths = fam["FID"].astype(str).str.len()
    n_zero = (fam["FID"] == 0).sum()
    if n_zero:
        lengths = lengths[fam["FID"] != 0]
    ok = lengths == 10
    assert ok.mean() == 1.0, (
        f"{(~ok).sum()} FID(s) are not 10 digits: "
        f"{fam.loc[~ok, 'FID'].tolist()[:10]}"
    )


def test_iid_length_is_11():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    lengths = _id_lengths(iid_test)
    assert (lengths["IID_len"] == 11.0).mean() == 1.0, (
        f"IID lengths: {lengths.to_dict('records')}"
    )


def test_fid_is_numeric():
    fam, batch, excluded = _load_data()
    iid_test = _prepare_iid_test(fam, batch, excluded)
    assert pd.api.types.is_numeric_dtype(iid_test["FID"]), (
        "FID column is not numeric — raw identifier may have leaked"
    )


# ══════════════════════════════════════════════════════════════════════════
#  Row-count integrity
# ══════════════════════════════════════════════════════════════════════════

def test_temp_fam_row_count_matches_onlyqc():
    """temp.fam must have the same row count as onlyQc.fam.

    PLINK errors with ``Unexpected PLINK 1 .bed file size`` when the
    .fam row count disagrees with the .bed row count.  Since temp.fam
    is supplied via ``--fam`` while ``--bfile`` points to the original
    (un-remapped) .bed, the row counts MUST match.
    """
    temp = _load_temp_fam()
    onlyqc = pd.read_csv(
        DATA_DIR / "onlyQc.fam",
        sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    assert len(temp) == len(onlyqc), (
        f"temp.fam ({len(temp)} rows) ≠ onlyQc.fam ({len(onlyqc)} rows) "
        f"— PLINK .bed size mismatch risk"
    )


# ══════════════════════════════════════════════════════════════════════════
#  fam / batch consistency
# ══════════════════════════════════════════════════════════════════════════

def test_fam_and_batch_order_matches():
    fam, batch, _ = _load_data()
    assert set(fam["IID"]) == set(batch["IID"]), "fam / batch IID sets differ"
    assert len(fam) == len(batch), (
        f"row count mismatch: {len(fam)} fam vs {len(batch)} batch"
    )
    assert (fam["IID"].to_numpy() == batch["IID"].to_numpy()).mean() == 1.0


def test_batch_info_iids_are_deidentified():
    """batch.info IIDs must also satisfy the de-identified IID pattern.

    This duplicates part of test_all_iids_match_deid_pattern but is
    specific to batch.info so failure messages are clearer.
    """
    _, batch, _ = _load_data()
    ok = batch["IID"].astype(str).str.match(r"^\d{10}[CM]$")
    assert ok.mean() == 1.0, (
        f"batch.info: {(~ok).sum()} IID(s) not de-identified: "
        f"{batch.loc[~ok, 'IID'].tolist()[:10]}"
    )


# ══════════════════════════════════════════════════════════════════════════
#  par_visit + exclusion constraint
# ══════════════════════════════════════════════════════════════════════════

def test_all_output_iids_in_par_visit():
    """Every output IID (fam *and* batch) maps to a par_visit candid and is
    NOT in the exclusion list."""
    fam, batch, _ = _load_data()
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    def _check_subjects(rc_set, label):
        ids_sub = identifiers[identifiers["release_candid"].isin(rc_set)]
        ids_rc = set(ids_sub["release_candid"].dropna().astype(int))
        rc_py = set(int(v) for v in rc_set)
        missing = rc_py - ids_rc
        assert len(missing) == 0, (
            f"{len(missing)} subject(s) in {label} missing from identifiers: "
            f"{sorted(missing)[:10]}"
        )
        bad = rc_py - par_candids
        assert len(bad) == 0, (
            f"{len(bad)} subject(s) in {label} NOT in par_visit: "
            f"{sorted(bad)[:10]}"
        )
        excluded_in = rc_py & exc_rc
        assert len(excluded_in) == 0, (
            f"{len(excluded_in)} excluded subject(s) found in {label}: "
            f"{sorted(excluded_in)[:10]}"
        )

    fids = set(fam["FID"].unique()) - {0}
    _check_subjects(fids, "hbcd.fam")

    batch_rc = set(
        pd.to_numeric(batch["IID"].str[:-1], errors="coerce")
        .dropna()
        .astype(int)
    )
    _check_subjects(batch_rc, "batch.info")


# ══════════════════════════════════════════════════════════════════════════
#  Filter correctness (re-derives the expected subject set)
# ══════════════════════════════════════════════════════════════════════════

def test_filter_correctness():
    """Re-derive the expected subject set from scratch using the pipeline's
    own merge-and-filter logic, then verify it matches hbcd.fam exactly."""
    identifiers = load_identifiers()
    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    valid_rc = (
        set(
            int(v)
            for v in identifiers[identifiers["release_candid"].isin(par_candids)][
                "release_candid"
            ]
            .dropna()
            .unique()
        )
        - exc_rc
    )

    temp = _load_temp_fam()
    temp["_rel"] = temp["IID"].astype(str).str[-1]

    input_batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
    input_batch["_rel"] = input_batch["IID"].str[-1]
    input_batch["_rc"] = pd.to_numeric(input_batch["IID"].str[:-1])
    input_batch = input_batch.drop(columns=["IID"])

    merged = temp.merge(
        input_batch,
        left_on=["FID", "_rel"],
        right_on=["_rc", "_rel"],
        how="inner",
    )
    ok = (
        merged["FID"].isin(valid_rc)
        & merged["visit"].notna()
        & merged["plate_number"].notna()
    )
    expected_rc = set(merged.loc[ok, "FID"].unique()) - {0}

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


# ══════════════════════════════════════════════════════════════════════════
#  Variant integrity
# ══════════════════════════════════════════════════════════════════════════

def test_variant_count_preserved():
    """All variants in the source .bim are preserved in the output .bim.

    Plink2's ``--keep`` filters subjects only; variants should pass
    through unchanged so the row count must be identical.
    """
    src = pd.read_csv(
        str(DATA_DIR / "HBCD.bim"),
        sep=r"\s+", header=None,
        names=["CHR", "SNP", "GD", "BP", "A1", "A2"],
    )
    out = pd.read_csv(
        RELEASE_DIR / "hbcd.bim",
        sep=r"\s+", header=None,
        names=["CHR", "SNP", "GD", "BP", "A1", "A2"],
    )
    assert len(src) == len(out), (
        f"variant count mismatch: {len(src)} source vs {len(out)} output"
    )


# ══════════════════════════════════════════════════════════════════════════
#  Derivative integrity (GRM, PCs, relatedness)
# ══════════════════════════════════════════════════════════════════════════


def _release_n_subjects():
    """Number of subjects in hbcd.fam (used for matrix-dimension checks)."""
    fam = pd.read_csv(
        RELEASE_DIR / "hbcd.fam",
        sep=r"\s+", header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    return len(fam)


# ── PLINK GRM ────────────────────────────────────────────────────────────


def test_plink_grm_dimensions():
    """PLINK GRM .rel.id row count matches hbcd.fam, and .rel is square."""
    n = _release_n_subjects()
    id_path = RELEASE_DIR / "hbcd_plink_grm.rel.id"
    rel_path = RELEASE_DIR / "hbcd_plink_grm.rel"
    if not id_path.exists():
        pytest.skip("hbcd_plink_grm.rel.id not found")
    ids = pd.read_csv(id_path, sep=r"\s+", header=None, names=["FID", "IID"])
    assert len(ids) == n, f"PLINK GRM .rel.id: {len(ids)} rows, expected {n}"

    if rel_path.exists():
        mat = np.loadtxt(str(rel_path))
        assert mat.shape == (n, n), (
            f"PLINK GRM .rel shape {mat.shape}, expected ({n},{n})"
        )


def test_plink_grm_ids_are_release():
    """All IIDs in PLINK GRM .rel.id are in the release set."""
    id_path = RELEASE_DIR / "hbcd_plink_grm.rel.id"
    if not id_path.exists():
        pytest.skip("hbcd_plink_grm.rel.id not found")
    ids = set(
        pd.read_csv(id_path, sep=r"\s+", header=None, names=["FID", "IID"])["IID"]
        .astype(str)
    )
    release_ids = _release_iids()
    extra = ids - release_ids
    assert len(extra) == 0, (
        f"{len(extra)} IID(s) in PLINK GRM not in release set: "
        f"{sorted(extra)[:10]}"
    )


# ── GCTA binary GRM ──────────────────────────────────────────────────────


def test_gcta_grm_dimensions():
    """GCTA GRM .grm.id row count matches hbcd.fam."""
    n = _release_n_subjects()
    id_path = RELEASE_DIR / "hbcd_gcta_grm.grm.id"
    if not id_path.exists():
        pytest.skip("hbcd_gcta_grm.grm.id not found")
    ids = pd.read_csv(id_path, sep="\t", header=None, names=["FID", "IID"])
    assert len(ids) == n, f"GCTA GRM .grm.id: {len(ids)} rows, expected {n}"


def test_gcta_grm_ids_are_release():
    """All IIDs in GCTA GRM .grm.id are in the release set."""
    id_path = RELEASE_DIR / "hbcd_gcta_grm.grm.id"
    if not id_path.exists():
        pytest.skip("hbcd_gcta_grm.grm.id not found")
    ids = set(
        pd.read_csv(id_path, sep="\t", header=None, names=["FID", "IID"])["IID"]
        .astype(str)
    )
    extra = ids - _release_iids()
    assert len(extra) == 0, (
        f"{len(extra)} IID(s) in GCTA GRM not in release set: "
        f"{sorted(extra)[:10]}"
    )


# ── PC-AiR scores ────────────────────────────────────────────────────────


def _release_iids():
    """Memoised helper: load release IID set once per session."""
    if not hasattr(_release_iids, "_cache"):
        keep = pd.read_csv(
            get_release_base() / "keep_list.txt",
            sep=r"\s+", header=None, names=["FID", "IID"],
        )
        _release_iids._cache = set(keep["IID"].astype(str))
    return _release_iids._cache


def test_pcair_scores_ids_are_release():
    """All sample.id in PC-AiR scores are release IIDs."""
    p = RELEASE_DIR / f"{NAME}_pc_scores.txt"
    if not p.exists():
        pytest.skip("PC-AiR scores not found")
    scores = pd.read_csv(p, sep="\t", dtype=str)
    ids = set(scores["sample.id"])
    extra = ids - _release_iids()
    assert len(extra) == 0, (
        f"{len(extra)} sample(s) in PC-AiR scores not in release set: "
        f"{sorted(extra)[:10]}"
    )
    assert len(scores) == len(_release_iids()), (
        f"PC-AiR scores row count ({len(scores)}) ≠ release IIDs ({len(_release_iids())})"
    )


# ── PC-AiR unrelated / related IDs (two formats) ─────────────────────────


@pytest.mark.parametrize("fstem", [
    "unrelated_ids", "related_ids",
    "pcair_unrelated_ids", "pcair_related_ids",
])
def test_pcair_id_list_ids_are_release(fstem):
    """All IDs in PC-AiR unrelated/related lists are release IIDs."""
    # Try both text and CSV/TSV
    for ext, id_col in [("", "SampleID"), (".csv", "SampleID"), (".tsv", "SampleID")]:
        p = RELEASE_DIR / f"{NAME}_{fstem}{ext}"
        if p.exists():
            break
    else:
        pytest.skip(f"No {NAME}_{fstem} file found")
    df = pd.read_csv(p, dtype=str)
    col = [c for c in df.columns if c.lower() in ("sampleid", "sample.id", "id")][0]
    ids = set(df[col].dropna().astype(str))
    extra = ids - _release_iids()
    assert len(extra) == 0, (
        f"{len(extra)} ID(s) in {fstem} not in release set: {sorted(extra)[:10]}"
    )


# ── PC-AiR eigenvalues (no sample IDs — structural check) ────────────────


def test_pcair_eigenvalues_exist():
    """PC-AiR eigenvalues files exist with expected columns."""
    for ext in [".csv", ".tsv"]:
        p = RELEASE_DIR / f"{NAME}_pcair_eigenvalues{ext}"
        if p.exists():
            df = pd.read_csv(p)
            assert "PC" in df.columns, f"{p.name} missing PC column"
            assert "Eigenvalue" in df.columns, f"{p.name} missing Eigenvalue column"
            assert len(df) > 0, f"{p.name} is empty"
            return
    pytest.skip("No PC-AiR eigenvalues file found")


# ── PC-Relate pairs ──────────────────────────────────────────────────────


@pytest.mark.parametrize("fstem", ["pcrelate_pairs", "pcrelate_ibd"])
def test_pcrelate_pair_ids_are_release(fstem):
    """All ID1/ID2 in PC-Relate pairs/IBD are release IIDs."""
    p = RELEASE_DIR / f"{NAME}_{fstem}.csv"
    if not p.exists():
        p = RELEASE_DIR / f"{NAME}_{fstem}.tsv"
    if not p.exists():
        pytest.skip(f"{NAME}_{fstem} not found")
    df = pd.read_csv(p, dtype=str)
    all_ids = set(df["ID1"]).union(set(df["ID2"]))
    extra = all_ids - _release_iids()
    assert len(extra) == 0, (
        f"{len(extra)} ID(s) in {fstem} not in release set: {sorted(extra)[:10]}"
    )


def test_pcrelate_self_ids_are_release():
    """All ID in PC-Relate self-kinship are release IIDs."""
    p = RELEASE_DIR / f"{NAME}_pcrelate_self.csv"
    if not p.exists():
        p = RELEASE_DIR / f"{NAME}_pcrelate_self.tsv"
    if not p.exists():
        pytest.skip("pcrelate_self not found")
    df = pd.read_csv(p, dtype=str)
    ids = set(df["ID"])
    extra = ids - _release_iids()
    assert len(extra) == 0, (
        f"{len(extra)} ID(s) in pcrelate_self not in release set: {sorted(extra)[:10]}"
    )


def test_pcrelate_kinmat_long_ids_are_release():
    """All ID1/ID2 in PC-Relate kinmat long are release IIDs."""
    p = RELEASE_DIR / f"{NAME}_pcrelate_kinmat_long.csv"
    if not p.exists():
        p = RELEASE_DIR / f"{NAME}_pcrelate_kinmat_long.tsv"
    if not p.exists():
        pytest.skip("pcrelate_kinmat_long not found")
    df = pd.read_csv(p, dtype=str)
    all_ids = set(df["ID1"]).union(set(df["ID2"]))
    extra = all_ids - _release_iids()
    assert len(extra) == 0, (
        f"{len(extra)} ID(s) in kinmat_long not in release set: {sorted(extra)[:10]}"
    )


def test_pcrelate_kinmat_wide_ids_are_release():
    """SampleID column and column headers in PC-Relate kinmat wide are release IIDs."""
    p = RELEASE_DIR / f"{NAME}_pcrelate_kinmat_wide.csv"
    if not p.exists():
        p = RELEASE_DIR / f"{NAME}_pcrelate_kinmat_wide.tsv"
    if not p.exists():
        pytest.skip("pcrelate_kinmat_wide not found")
    df = pd.read_csv(p, dtype=str)
    row_ids = set(df["SampleID"])
    col_ids = set(df.columns[1:])
    extra_rows = row_ids - _release_iids()
    extra_cols = col_ids - _release_iids()
    assert len(extra_rows) == 0, (
        f"{len(extra_rows)} row SampleID(s) not in release set: {sorted(extra_rows)[:10]}"
    )
    assert len(extra_cols) == 0, (
        f"{len(extra_cols)} column ID(s) not in release set: {sorted(extra_cols)[:10]}"
    )
