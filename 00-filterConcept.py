#!/usr/bin/env python3

"""Explore overlap between data sources using Venn diagrams.

Compares release_candids across:
  - identifiers (release_identifiers CSV)
  - par_visit (completed visit data)
  - HST .fam (de-identified PLINK file from HST_HBCD_Transfer)
"""

import pandas as pd
import numpy as np
from pathlib import Path

from _lib import (
    DATA_DIR,
    load_par_visit_candids,
    load_identifiers,
    load_excluded_release_candids,
)

try:
    HERE = Path(__file__).resolve().parent
except NameError:
    HERE = Path.cwd()

try:
    from matplotlib_venn import venn3
    import matplotlib.pyplot as plt
    HAS_VENN = True
except ImportError:
    HAS_VENN = False
    print("matplotlib-venn not available — skipping Venn diagram.")

# ── data sources ──────────────────────────────────────────────────────────
identifiers = load_identifiers()
par_candids = load_par_visit_candids()
exc_rc = load_excluded_release_candids()

hst_fam = pd.read_csv(
    str(DATA_DIR / ".." / "HST_HBCD_Transfer_July2025" / "HBCD_analysis" / "hbcd.fam"),
    sep=r"\s+",
    header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
hst_fam["pscid"] = (
    hst_fam["FID"].astype(str)
    .str.rsplit("_", n=1).str[-1]
    .str[:-1]                # drop trailing relationship letter
    .pipe(pd.to_numeric, errors="coerce")
)

# ── deduplicate + build key sets ───────────────────────────────────────────
identifiers = identifiers.drop_duplicates(subset="release_candid")

# Enforce integer pscid for joining
identifiers["pscid"] = pd.to_numeric(identifiers["pscid"], errors="coerce")

hst_with_rc = hst_fam.merge(
    identifiers[["pscid", "release_candid"]].drop_duplicates(subset="pscid"),
    on="pscid",
    how="left",
)

id_rc = set(identifiers["release_candid"].unique())
hst_rc = set(hst_with_rc["release_candid"].dropna().unique())

# ── print counts ───────────────────────────────────────────────────────────
print("═" * 60)
print("Set sizes (release_candid level)")
print("═" * 60)
print(f"  identifiers                           : {len(id_rc):>6}")
print(f"  HST .fam (hbcd.fam, excl FID=0)       : {len(hst_rc):>6}")
print(f"  par_visit (completed visits)           : {len(par_candids):>6}")
print(f"  excluded                               : {len(exc_rc):>6}")
print()

print("═" * 60)
print("Pairwise overlaps")
print("═" * 60)
print(f"  identifiers ∩ HST .fam                 : {len(id_rc & hst_rc):>6}")
print(f"  identifiers ∩ par_visit                : {len(id_rc & par_candids):>6}")
print(f"  HST .fam ∩ par_visit                   : {len(hst_rc & par_candids):>6}")
print()

print("═" * 60)
print("Triple overlap & effects of exclusion")
print("═" * 60)
triple = id_rc & hst_rc & par_candids
print(f"  identifiers ∩ HST .fam ∩ par_visit     : {len(triple):>6}")
print(f"    minus excluded                       : {len(triple - exc_rc):>6}")
print()

id_not_hst = id_rc - hst_rc
hst_not_id = hst_rc - id_rc
print(f"  identifiers only (not in HST .fam)     : {len(id_not_hst):>6}")
print(f"  HST .fam only (not in identifiers)     : {len(hst_not_id):>6}")

# ── Venn diagram ───────────────────────────────────────────────────────────
if HAS_VENN:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    subsets = (
        len(id_rc - hst_rc - par_candids),
        len(hst_rc - id_rc - par_candids),
        len((id_rc & hst_rc) - par_candids),
        len(par_candids - id_rc - hst_rc),
        len((id_rc & par_candids) - hst_rc),
        len((hst_rc & par_candids) - id_rc),
        len(id_rc & hst_rc & par_candids),
    )

    colors = ["#1f78b4", "#e31a1c", "#33a02c"]

    v = venn3(subsets, set_labels=("identifiers", "HST .fam", "par_visit"),
              set_colors=colors, ax=ax1)
    ax1.set_title("release_candid overlap", fontsize=12)

    v2 = venn3(subsets, set_labels=("identifiers", "HST .fam", "par_visit"),
               set_colors=colors, ax=ax2)
    ax2.set_title("with excluded highlighted", fontsize=12)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=colors[0], alpha=0.5, label="identifiers"),
        Patch(facecolor=colors[1], alpha=0.5, label="HST .fam"),
        Patch(facecolor=colors[2], alpha=0.5, label="par_visit"),
    ]
    ax1.legend(handles=legend_elements, loc="lower left", fontsize=9)
    ax2.legend(handles=legend_elements, loc="lower left", fontsize=9)

    valid_final = triple - exc_rc
    ax2.text(
        -0.6, -0.7,
        f"Excluded: {len(exc_rc)}\nValid (final): {len(valid_final)}",
        fontsize=10,
        bbox=dict(facecolor="lightcoral", alpha=0.4),
    )

    out_path = HERE / "filter_overlap_venn.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nVenn diagram saved to {out_path}")
    plt.show()

else:
    print("\n═" * 60)
    print("Set sizes (release_candid level)")
    print("═" * 60)
    print(f"  identifiers           = {len(id_rc)}")
    print(f"  HST .fam              = {len(hst_rc)}")
    print(f"  par_visit             = {len(par_candids)}")
    print(f"  excluded              = {len(exc_rc)}")
    print(f"  identifiers ∩ HST     = {len(id_rc & hst_rc)}")
    print(f"  (∩) ∩ par_visit      = {len(id_rc & hst_rc & par_candids)}")
    print(f"  valid (final)         = {len((id_rc & hst_rc & par_candids) - exc_rc)}")
