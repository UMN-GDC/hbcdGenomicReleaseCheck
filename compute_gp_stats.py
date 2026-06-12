#!/usr/bin/env python3
"""
Parse bcftools query GP output from stdin, compute per-sample and
ancestry-stratified GP statistics.

Input (stdin): one variant per line:
  CHR\\tPOS\\tGP1\\tGP2\\t...\\tGPn
  GPi = "0.997,0.003,0.000"  (3 comma-separated genotype probabilities)

Output TSV files to --outdir:
  per_sample.tsv   — per-sample aggregates
  hist_bins.tsv    — per-ancestry max-GP histograms
  ridge.tsv        — subsampled ridge observations
  call_rate.tsv    — subsampled per-variant call rates
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ancestry", required=True,
                        help="sample_ancestry.txt (tab: position, ancestry)")
    parser.add_argument("--chr", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--ridge-rate", type=float, default=0.002)
    parser.add_argument("--call-rate-rate", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    anc_df = pd.read_csv(args.ancestry, sep="\t", header=None,
                         names=["pos", "ancestry"])
    ancestry = anc_df["ancestry"].values
    n_samples = len(ancestry)

    unique_anc = list(dict.fromkeys(ancestry))
    anc_to_idx = {a: i for i, a in enumerate(unique_anc)}
    n_anc = len(unique_anc)

    sum_maxgp = np.zeros(n_samples, dtype=np.float64)
    cnt = np.zeros(n_samples, dtype=np.int64)
    conf = np.zeros(n_samples, dtype=np.int64)
    hist = np.zeros((n_anc, 50), dtype=np.int64)

    ridge_records = []
    call_rate_records = []
    n_variants = 0

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        gp_strings = parts[2:]
        n_variants += 1

        anc_tot = {}
        anc_call = {}

        for i, gp_str in enumerate(gp_strings):
            gp = np.fromstring(gp_str, sep=",", dtype=np.float64)
            max_gp = gp.max()
            a = ancestry[i]
            sum_maxgp[i] += max_gp
            cnt[i] += 1
            if max_gp > 0.9:
                conf[i] += 1
            bi = min(int(max_gp * 50), 49)
            hist[anc_to_idx[a], bi] += 1
            if rng.random() < args.ridge_rate:
                ridge_records.append((chrom, a, max_gp))
            anc_tot[a] = anc_tot.get(a, 0) + 1
            if max_gp > 0.9:
                anc_call[a] = anc_call.get(a, 0) + 1

        if rng.random() < args.call_rate_rate:
            for a in anc_tot:
                call_rate_records.append(
                    (chrom, pos, a, anc_call.get(a, 0), anc_tot[a])
                )

    # ── Write outputs ──────────────────────────────────────────────
    per_sample = pd.DataFrame({
        "sample_pos": np.arange(1, n_samples + 1),
        "chr": args.chr,
        "ancestry": ancestry,
        "sum_maxgp": sum_maxgp,
        "n_variants": cnt,
        "n_confident": conf,
    })
    per_sample.to_csv(outdir / "per_sample.tsv", sep="\t", index=False)

    hist_rows = []
    for ai, a in enumerate(unique_anc):
        for bi in range(50):
            if hist[ai, bi] > 0:
                hist_rows.append((a, args.chr, bi + 1, int(hist[ai, bi])))
    hist_df = pd.DataFrame(hist_rows,
                           columns=["ancestry", "chr", "bin", "count"])
    hist_df.to_csv(outdir / "hist_bins.tsv", sep="\t", index=False)

    ridge_df = pd.DataFrame(ridge_records,
                            columns=["chromosome", "ancestry", "maxgp"])
    ridge_df.to_csv(outdir / "ridge.tsv", sep="\t", index=False)

    cr_df = pd.DataFrame(call_rate_records,
                         columns=["chr", "pos", "ancestry", "called", "total"])
    cr_df.to_csv(outdir / "call_rate.tsv", sep="\t", index=False)

    print(f"  Done: {n_variants} variants, {n_samples} samples, "
          f"{len(ridge_records)} ridge, {len(call_rate_records)} call-rate",
          file=sys.stderr)


if __name__ == "__main__":
    main()
