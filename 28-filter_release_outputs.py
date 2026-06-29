#!/usr/bin/env python3

# step 28: filter all derivative outputs (GRM/PC-AiR/PC-Relate) to release subjects
# Reads keep_list.txt (release IIDs) and writes filtered copies to the
# release data directory, preserving the full-data originals in DATA_DIR.
#
# Filtered outputs (all sourced from HANDOFF_DIR except pc_scores from DATA_DIR):
#   - GCTA binary GRM          →  release_dir/hbcd_gcta_grm.{grm.id,grm.bin,grm.N.bin}
#   - PLINK GRM                →  release_dir/hbcd_plink_grm.{rel,rel.id}
#   - PC-AiR scores            →  release_dir/hbcd_rsid_harmonized_pc_scores.txt
#   - PC-Relate pairs/IBD/self →  release_dir/hbcd_rsid_harmonized_pcrelate_{pairs,ibd,self}.{csv,tsv}
#   - PC-Relate kinmat         →  release_dir/hbcd_rsid_harmonized_pcrelate_kinmat_{long,wide}.{csv,tsv}
#   - PC-Relate GRM            →  release_dir/hbcd_pcrelate_grm.{grm.id,grm.bin,grm.N.bin,gz}
#   - PC-Relate GRM pairwise   →  release_dir/hbcd_pcrelate_grm_pairwise.tsv
#   - PC-Relate relatedness    →  release_dir/hbcd_pcrelate_relatedness.tsv
#   - PC-AiR 32 PCs clean      →  release_dir/hbcd_pcair_32PCs_clean.tsv
#   - PLINK PCA 32 PCs         →  release_dir/hbcd_plink_pca_32PCs.tsv
#   - KING relatedness         →  release_dir/hbcd_king_relatedness.tsv
#
# Usage: python 28-filter_release_outputs.py [--release-dir /path]
# Default release from HBCD_RELEASE env var or br_21p2.

import sys
import shutil
import pandas as pd
import numpy as np
from pathlib import Path
from _lib import DATA_DIR, get_release_dir

HANDOFF_DIR = DATA_DIR.parent / "data_handoff"
NAME = "hbcd_rsid_harmonized"

# Accept --release-dir / -r override, default from env/HBCD_RELEASE
if len(sys.argv) > 1 and sys.argv[1] in ("--release-dir", "-r"):
    RELEASE_DIR = Path(sys.argv[2])
else:
    RELEASE_DIR = get_release_dir()

release_base = RELEASE_DIR.parent

# -- load release IIDs --
keep = pd.read_csv(
    release_base / "keep_list.txt",
    sep=r"\s+",
    header=None,
    names=["FID", "IID"],
)
release_iids = set(keep["IID"].astype(str))
print(f"Release IIDs: {len(release_iids)}")


# ── filter helpers (always write to out_path; leave source untouched) ────


def filter_text(in_path, id_col, out_path):
    """Read a delimited text file, keep rows where *id_col* is in
    ``release_iids``, write filtered result to *out_path*."""
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, sep="\t", dtype=str)
    before = len(df)
    df = df[df[id_col].isin(release_iids)]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_matrix(in_rel, in_id, out_rel, out_id):
    """Filter a square matrix (.rel + .rel.id) to only release subjects."""
    if not in_rel.exists() or not in_id.exists():
        print(f"  SKIP (not found): {in_rel.name} / {in_id.name}")
        return
    ids = pd.read_csv(in_id, sep=r"\s+", header=None, names=["FID", "IID"])
    orig_n = len(ids)
    keep_mask = ids["IID"].astype(str).isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    ids_filtered = ids[keep_mask].reset_index(drop=True)
    out_id.parent.mkdir(parents=True, exist_ok=True)
    ids_filtered.to_csv(out_id, sep=" ", header=False, index=False)

    mat = np.loadtxt(str(in_rel))
    mat_filtered = mat[np.ix_(keep_idx, keep_idx)]
    np.savetxt(str(out_rel), mat_filtered, fmt="%.10g")
    print(f"  {in_rel.name}: {orig_n}x{orig_n} → {len(ids_filtered)}x{len(ids_filtered)}")


def filter_binary_grm(in_prefix, out_prefix, id_ext=".grm.id", bin_ext=".grm.bin", n_ext=".grm.N.bin"):
    """Filter a binary GRM to release subjects.

    Parameters
    ----------
    in_prefix, out_prefix : str or Path
        Prefix before the extension parts.
    id_ext : str
        Extension for the ID file (default ``.grm.id``).
    bin_ext : str
        Extension for the binary GRM file (default ``.grm.bin``).
    n_ext : str
        Extension for the N-obs file (default ``.grm.N.bin``).
    """
    id_in = Path(f"{in_prefix}{id_ext}")
    bin_in = Path(f"{in_prefix}{bin_ext}")
    n_in = Path(f"{in_prefix}{n_ext}")
    id_out = Path(f"{out_prefix}{id_ext}")
    bin_out = Path(f"{out_prefix}{bin_ext}")
    n_out = Path(f"{out_prefix}{n_ext}")

    if not id_in.exists() or not bin_in.exists():
        print(f"  SKIP (not found): {in_prefix}.grm.*")
        return

    ids = pd.read_csv(id_in, sep="\t", header=None, names=["FID", "IID"], dtype=str)
    n = len(ids)
    keep_mask = ids["IID"].isin(release_iids)
    keep_idx = [i for i, k in enumerate(keep_mask) if k]
    n_keep = len(keep_idx)
    print(f"  GRM: {n} → {n_keep} subjects")

    id_out.parent.mkdir(parents=True, exist_ok=True)
    ids_filtered = ids[keep_mask].reset_index(drop=True)
    ids_filtered.to_csv(id_out, sep="\t", header=False, index=False)

    dt = np.dtype("f4")
    grm_flat = np.fromfile(str(bin_in), dtype=dt)
    grm_full = np.zeros((n, n), dtype=dt)
    grm_full[np.tril_indices(n)] = grm_flat
    grm_full += grm_full.T
    np.fill_diagonal(grm_full, np.diag(grm_full) / 2)
    grm_sub = grm_full[np.ix_(keep_idx, keep_idx)]
    grm_sub[np.tril_indices(n_keep)].tofile(str(bin_out))
    print(f"  {bin_in.name}: {n}x{n} → {n_keep}x{n_keep}")

    if n_in.exists():
        n_flat = np.fromfile(str(n_in), dtype=dt)
        n_full = np.zeros((n, n), dtype=dt)
        n_full[np.tril_indices(n)] = n_flat
        n_full += n_full.T
        np.fill_diagonal(n_full, np.diag(n_full) / 2)
        n_sub = n_full[np.ix_(keep_idx, keep_idx)]
        n_sub[np.tril_indices(n_keep)].tofile(str(n_out))
        print(f"  {n_in.name}: {n}x{n} → {n_keep}x{n_keep}")


def filter_csv(in_path, out_path, id_cols, sep=","):
    """Filter a CSV/TSV keeping rows where any *id_cols* matches a release IID."""
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, sep=sep, dtype=str)
    before = len(df)
    mask = df[id_cols].apply(lambda c: c.isin(release_iids)).any(axis=1)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"  {in_path.name}: {before} → {len(df)} rows")


def filter_kinmat_wide(in_path, out_path, sep=","):
    """Filter a wide-format kinship matrix CSV keeping only release columns/rows."""
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    df = pd.read_csv(in_path, sep=sep, dtype=str)
    cols = ["SampleID"] + [c for c in df.columns[1:] if c in release_iids]
    df = df[df["SampleID"].isin(release_iids)][cols]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"  {in_path.name}: {len(df)} samples x {len(cols)-1} cols")


def filter_grm_text_gz(in_path, out_path):
    """Filter a gzipped text GRM (GCTA ``--make-grm-gz`` format).

    Columns assumed: ``FID1 IID1 FID2 IID2 N GRM``, space-separated, no
    header.  Only rows where both **IID1** and **IID2** are in
    ``release_iids`` are kept.
    """
    if not in_path.exists():
        print(f"  SKIP (not found): {in_path}")
        return
    # Read gzipped text with no header — ID columns are str, N/GRM are float
    df = pd.read_csv(
        in_path, sep=r"\s+", header=None,
        names=["FID1", "IID1", "FID2", "IID2", "N", "GRM"],
        dtype={"FID1": str, "IID1": str, "FID2": str, "IID2": str},
    )
    before = len(df)
    mask = df["IID1"].isin(release_iids) & df["IID2"].isin(release_iids)
    df = df[mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep=" ", header=False, index=False, compression="gzip")
    print(f"  {in_path.name}: {before} → {len(df)} rows")


# ===========================================================================
# GCTA binary GRM
# ===========================================================================
print("\n--- GCTA binary GRM ---")
filter_binary_grm(
    str(HANDOFF_DIR / "hbcd_gcta_grm"),
    str(RELEASE_DIR / "hbcd_gcta_grm"),
)

# ===========================================================================
# PLINK GRM
# ===========================================================================
print("\n--- PLINK GRM ---")
filter_matrix(
    HANDOFF_DIR / "hbcd_plink_grm.rel",
    HANDOFF_DIR / "hbcd_plink_grm.rel.id",
    RELEASE_DIR / "hbcd_plink_grm.rel",
    RELEASE_DIR / "hbcd_plink_grm.rel.id",
)

# ===========================================================================
# PC-AiR text exports
# ===========================================================================
print("\n--- PC-AiR ---")
filter_text(
    DATA_DIR / f"{NAME}_pc_scores.txt", "sample.id",
    RELEASE_DIR / f"{NAME}_pc_scores.txt",
)

# ===========================================================================
# PC-AiR eigenvalues (no sample IDs — just copy)
# ===========================================================================
print("\n--- PC-AiR eigenvalues ---")
for fname in [f"{NAME}_pcair_eigenvalues.csv", f"{NAME}_pcair_eigenvalues.tsv"]:
    src = DATA_DIR / fname
    if src.exists():
        shutil.copy2(src, RELEASE_DIR / fname)
        print(f"  Copied: {fname}")

# ===========================================================================
# PC-AiR scores (CSV/TSV from export_pcrelate.R)
# ===========================================================================
print("\n--- PC-AiR scores (export) ---")
for ext in ["csv", "tsv"]:
    filter_csv(
        DATA_DIR / f"{NAME}_pcair_scores.{ext}",
        RELEASE_DIR / f"{NAME}_pcair_scores.{ext}",
        ["SampleID"],
        sep="\t" if ext == "tsv" else ",",
    )

# ===========================================================================
# PC-AiR unrelated/related IDs (export format)
# ===========================================================================
print("\n--- PC-AiR unrelated/related (export) ---")
for fstem in ["pcair_unrelated_ids", "pcair_related_ids"]:
    for ext in ["csv", "tsv"]:
        filter_csv(
            DATA_DIR / f"{NAME}_{fstem}.{ext}",
            RELEASE_DIR / f"{NAME}_{fstem}.{ext}",
            ["SampleID"],
            sep="\t" if ext == "tsv" else ",",
        )

# ===========================================================================
# PC-Relate pairs + IBD + self
# ===========================================================================
print("\n--- PC-Relate pairs ---")
for ext in ["csv", "tsv"]:
    sep = "\t" if ext == "tsv" else ","
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_pairs.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_pairs.{ext}",
        ["ID1", "ID2"], sep=sep,
    )
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_ibd.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_ibd.{ext}",
        ["ID1", "ID2"], sep=sep,
    )
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_self.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_self.{ext}",
        ["ID"], sep=sep,
    )

# ===========================================================================
# PC-Relate kinship matrices
# ===========================================================================
print("\n--- PC-Relate kinmat long ---")
for ext in ["csv", "tsv"]:
    filter_csv(
        DATA_DIR / f"{NAME}_pcrelate_kinmat_long.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_kinmat_long.{ext}",
        ["ID1", "ID2"],
        sep="\t" if ext == "tsv" else ",",
    )

print("\n--- PC-Relate kinmat wide ---")
for ext in ["csv", "tsv"]:
    filter_kinmat_wide(
        DATA_DIR / f"{NAME}_pcrelate_kinmat_wide.{ext}",
        RELEASE_DIR / f"{NAME}_pcrelate_kinmat_wide.{ext}",
        sep="\t" if ext == "tsv" else ",",
    )

# ===========================================================================
# PC-Relate GRM (binary — pcrelate_grm naming, no extra .grm. infix)
# ===========================================================================
print("\n--- PC-Relate GRM (binary) ---")
filter_binary_grm(
    str(HANDOFF_DIR / "hbcd_pcrelate_grm"),
    str(RELEASE_DIR / "hbcd_pcrelate_grm"),
    id_ext=".grm.id", bin_ext=".grm.bin", n_ext=".grm.N.bin",
)

# ===========================================================================
# PC-Relate GRM (gzipped text)
# ===========================================================================
print("\n--- PC-Relate GRM (text gz) ---")
filter_grm_text_gz(
    HANDOFF_DIR / "hbcd_pcrelate_grm.gz",
    RELEASE_DIR / "hbcd_pcrelate_grm.gz",
)

# ===========================================================================
# PC-Relate GRM pairwise
# ===========================================================================
print("\n--- PC-Relate GRM pairwise ---")
filter_csv(
    HANDOFF_DIR / "hbcd_pcrelate_grm_pairwise.tsv",
    RELEASE_DIR / "hbcd_pcrelate_grm_pairwise.tsv",
    ["ID1", "ID2"], sep="\t",
)

# ===========================================================================
# PC-Relate relatedness
# ===========================================================================
print("\n--- PC-Relate relatedness ---")
filter_csv(
    HANDOFF_DIR / "hbcd_pcrelate_relatedness.tsv",
    RELEASE_DIR / "hbcd_pcrelate_relatedness.tsv",
    ["subject_id_1", "subject_id_2"], sep="\t",
)

# ===========================================================================
# PC-AiR 32 PCs clean
# ===========================================================================
print("\n--- PC-AiR 32 PCs clean ---")
filter_text(
    HANDOFF_DIR / "hbcd_pcair_32PCs_clean.tsv",
    "participant_id",
    RELEASE_DIR / "hbcd_pcair_32PCs_clean.tsv",
)

# ===========================================================================
# PLINK PCA 32 PCs
# ===========================================================================
print("\n--- PLINK PCA 32 PCs ---")
filter_text(
    HANDOFF_DIR / "hbcd_plink_pca_32PCs.tsv",
    "participant_id",
    RELEASE_DIR / "hbcd_plink_pca_32PCs.tsv",
)

# ===========================================================================
# KING relatedness
# ===========================================================================
print("\n--- KING relatedness ---")
filter_csv(
    HANDOFF_DIR / "hbcd_king_relatedness.tsv",
    RELEASE_DIR / "hbcd_king_relatedness.tsv",
    ["subject_id_1", "subject_id_2"], sep="\t",
)

print("\nDone.")
