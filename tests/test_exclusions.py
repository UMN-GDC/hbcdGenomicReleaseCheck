"""Ensure all exclusion-list pscids are absent from the release subjects."""

from _lib import (
    DATA_DIR,
    load_identifiers,
    parse_additional_exclusion_lists,
    ADDITIONAL_EXCLUSIONS_FILE,
)


def test_all_exclusion_lists_have_zero_overlap():
    identifiers = load_identifiers()
    identifiers["pscid"] = identifiers["pscid"].astype(str).str.strip()
    release_pscids = set(identifiers["pscid"].unique())

    path = DATA_DIR / ADDITIONAL_EXCLUSIONS_FILE
    excl_lists = parse_additional_exclusion_lists(path)

    assert len(excl_lists) > 0, "No exclusion lists were parsed"

    failures = []
    for header, pscids in excl_lists.items():
        overlap = release_pscids & pscids
        if overlap:
            failures.append(f"{header}: {sorted(overlap)[:10]}")

    assert not failures, (
        f"{len(failures)} exclusion list(s) overlap with release subjects:\n"
        + "\n".join(failures)
    )
