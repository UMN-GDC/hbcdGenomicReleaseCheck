#!/usr/bin/env python3

"""Explore overlap between data sources and filters using Venn diagrams.

Reuses the same merge + filter logic as 01-filterGenotypeFiles.py but
visualises how many release_candids / pscids survive each stage.
"""

import pandas as pd
import numpy as np
from pathlib import Path

from _lib import (
    DATA_DIR,
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
    load_excluded_with_relationship,
)

try:
    from matplotlib_venn import venn3, venn3_circles
    import matplotlib.pyplot as plt

    HAS_VENN = True
except ImportError:
    HAS_VENN = False
    print("matplotlib-venn not available — skipping Venn diagram.")


def main():
    # ── data sources ──────────────────────────────────────────────────────
    identifiers = load_identifiers()
    batch = pd.read_csv(DATA_DIR / "batch.info", sep=r"\s+")
    batch["release_candid"] = pd.to_numeric(batch["IID"].str[:-1])
    fam = pd.read_csv(
        str(DATA_DIR / "HBCD.fam"),
        sep=r"\s+",
        header=None,
        names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
    )
    fam["pscid"] = fam["IID"].str[-10:-1]

    par_candids = load_par_visit_candids()
    exc_rc = load_excluded_release_candids()

    # ── build key sets (by release_candid) ─────────────────────────────────
    id_rc = set(identifiers["release_candid"].unique())
    batch_rc = set(batch["release_candid"].dropna().unique())
    fam_pscid = set(fam["pscid"].unique())
    id_pscid = set(identifiers["pscid"].unique())

    # release_candids corresponding to fam pscids
    fam_rc = set(
        identifiers.loc[identifiers["pscid"].isin(fam_pscid), "release_candid"].unique()
    )

    # ── print stage-by-stage counts ────────────────────────────────────────
    print("═" * 60)
    print("Filtering cascade (release_candid level)")
    print("═" * 60)
    print(f"  identifiers (unique release_candids)   : {len(id_rc):>6}")
    print(f"  batch.info (unique release_candids)    : {len(batch_rc):>6}")
    print(f"  par_visit (completed visits)           : {len(par_candids):>6}")
    print(f"  excluded                               : {len(exc_rc):>6}")
    print()

    # identifiers ∩ batch
    id_and_batch = id_rc & batch_rc
    id_only = id_rc - batch_rc
    batch_only = batch_rc - id_rc
    print(f"  identifiers + batch overlap            : {len(id_and_batch):>6}")
    print(f"  identifiers only                       : {len(id_only):>6}")
    print(f"  batch only                             : {len(batch_only):>6}")
    print()

    # in par_visit
    in_par = id_and_batch & par_candids
    not_par = id_and_batch - par_candids
    print(f"  overlapping + in par_visit             : {len(in_par):>6}")
    print(f"  overlapping but NOT in par_visit        : {len(not_par):>6}")
    print()

    # minus excluded
    valid_rc = in_par - exc_rc
    excluded_in_valid = in_par & exc_rc
    print(f"  valid (par_visit - excluded)            : {len(valid_rc):>6}")
    print(f"  excluded from overlapping + par_visit   : {len(excluded_in_valid):>6}")
    print()

    # pscid overlap (fam)
    print("═" * 60)
    print("pscid-level overlap (fam × identifiers)")
    print("═" * 60)
    print(f"  fam (unique pscids)                    : {len(fam_pscid):>6}")
    print(f"  identifiers (unique pscids)            : {len(id_pscid):>6}")
    print(f"  fam ∩ identifiers                      : {len(fam_pscid & id_pscid):>6}")
    print(f"  fam only                               : {len(fam_pscid - id_pscid):>6}")
    print(f"  identifiers only                       : {len(id_pscid - fam_pscid):>6}")
    print()

    # ── Venn diagram (release_candid) ──────────────────────────────────────
    if HAS_VENN:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

        # Panel 1: identifiers × batch × par_visit
        v = venn3(
            subsets=(len(id_rc - batch_rc - par_candids),
                     len(batch_rc - id_rc - par_candids),
                     len((id_rc & batch_rc) - par_candids),
                     len(par_candids - id_rc - batch_rc),
                     len((id_rc & par_candids) - batch_rc),
                     len((batch_rc & par_candids) - id_rc),
                     len(id_rc & batch_rc & par_candids)),
            set_labels=("identifiers", "batch.info", "par_visit"),
            ax=ax1,
        )
        ax1.set_title("release_candid overlap", fontsize=12)
        for text in v.set_labels:
            if text:
                text.set_fontsize(10)
        if v.subset_labels:
            for text in v.subset_labels:
                if text:
                    text.set_fontsize(9)

        # Panel 2: add excluded annotation
        v2 = venn3(
            subsets=(len(id_rc - batch_rc - par_candids),
                     len(batch_rc - id_rc - par_candids),
                     len((id_rc & batch_rc) - par_candids),
                     len(par_candids - id_rc - batch_rc),
                     len((id_rc & par_candids) - batch_rc),
                     len((batch_rc & par_candids) - id_rc),
                     len(id_rc & batch_rc & par_candids)),
            set_labels=("identifiers", "batch.info", "par_visit"),
            ax=ax2,
        )
        ax2.set_title("with excluded highlighted", fontsize=12)
        for text in v2.set_labels:
            if text:
                text.set_fontsize(10)
        if v2.subset_labels:
            for text in v2.subset_labels:
                if text:
                    text.set_fontsize(9)

        # highlight the excluded count in the triple-overlap region
        triple_rc = id_rc & batch_rc & par_candids
        exc_in_triple = triple_rc & exc_rc
        valid_final = triple_rc - exc_rc
        # find the patch for the triple-overlap region (index 6, the 7th subset)
        # venn3 subset order: 100,010,110,001,101,011,111  (binary-coded)
        # region 111 is subsets[6]
        if len(exc_rc) > 0:
            # annotate in the bottom-left corner of panel 2
            ax2.text(
                -0.6, -0.7,
                f"Excluded: {len(exc_rc)}\n"
                f"Valid (final): {len(valid_final)}",
                fontsize=10,
                bbox=dict(facecolor="lightcoral", alpha=0.4),
            )
        # Also show valid vs excluded as a small inset
        ax2.annotate(
            f"{len(valid_final)} valid",
            xy=(0, 0), xytext=(0.5, -0.6),
            fontsize=11, ha="center",
            arrowprops=dict(arrowstyle="->", color="green"),
        )

        out_path = Path(__file__).parent / "filter_overlap_venn.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        print(f"Venn diagram saved to {out_path}")
        plt.show()
    else:
        # text-only fallback
        print("═" * 60)
        print("Set sizes (release_candid level)")
        print("═" * 60)
        print(f"  identifiers           = {len(id_rc)}")
        print(f"  batch.info            = {len(batch_rc)}")
        print(f"  par_visit             = {len(par_candids)}")
        print(f"  excluded              = {len(exc_rc)}")
        print(f"  identifiers ∩ batch   = {len(id_rc & batch_rc)}")
        print(f"  (∩) ∩ par_visit      = {len(id_rc & batch_rc & par_candids)}")
        print(f"  valid (final)         = {len((id_rc & batch_rc & par_candids) - exc_rc)}")


if __name__ == "__main__":
    main()
