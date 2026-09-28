#!/usr/bin/env python3
"""
Build an analysis-ready, de-duplicated copy of the ported ESM-2 8M baseline.

The ported reference (`data/reference/esm2_8m_pair_distances.csv`) is kept BYTE-FOR-BYTE
as it came from the previous repo, warts and all, for provenance. It has two known issues
(both artifacts of the previous repo's batched, append-mode PLM run — not distance errors):

  1. One corrupted value: `4.547834396362305e-05git` (a stray "git" token appended),
     which forces the whole distance column to dtype=str.
  2. ~18k pairs written twice (append-mode runs with overlapping row ranges). Duplicates
     agree on distance to floating-point noise (max spread ~2.4e-7).

This script produces `*_dedup.csv`: numeric distances, one row per unique ortholog pair.
It does NOT modify the reference file.

IMPORTANT (methodological note carried into provenance):
The previous paper's published percentiles (median 0.000487, 95th 0.008432, 99th 0.039693)
were computed over the FULL row set (68,971 rows, duplicates included) and reproduce exactly
that way. Because the duplicated pairs are disproportionately low-divergence, per-unique-pair
percentiles on THIS dedup file are higher (≈0.000496 / 0.009815 / 0.050938). Use the dedup
file for per-pair analyses; do not expect it to reproduce the published row-wise percentiles.

Usage:
    python analysis/build_baseline_dedup.py \\
        --in data/reference/esm2_8m_pair_distances.csv \\
        --out data/reference/esm2_8m_pair_distances_dedup.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", default="data/reference/esm2_8m_pair_distances.csv")
    ap.add_argument("--out", dest="out", default="data/reference/esm2_8m_pair_distances_dedup.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    n_raw = len(df)

    # (1) coerce distances; the lone corrupted token becomes NaN, then we recover its
    #     numeric value from its byte-identical duplicate row if present.
    df["plm_cosine_distance"] = pd.to_numeric(df["plm_cosine_distance"], errors="coerce")
    n_nan = int(df["plm_cosine_distance"].isna().sum())

    # (2) collapse duplicate pairs. Sort NaN last within each pair so keep="first" prefers
    #     a numeric value when a pair has both a good and a corrupted row.
    df = df.sort_values("plm_cosine_distance", na_position="last")
    ded = df.drop_duplicates(subset=["zebrafish_id", "aescallii_id"], keep="first")
    ded = ded.sort_values(["zebrafish_id", "aescallii_id"]).reset_index(drop=True)

    still_nan = int(ded["plm_cosine_distance"].isna().sum())

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ded.to_csv(out, index=False)

    print(f"in : {n_raw} rows ({n_nan} non-numeric after coercion)")
    print(f"out: {len(ded)} unique pairs -> {out}  ({still_nan} still NaN)")
    d = ded["plm_cosine_distance"].dropna()
    print(f"per-pair distance: median={d.median():.6f} 95th={d.quantile(0.95):.6f} 99th={d.quantile(0.99):.6f}")


if __name__ == "__main__":
    main()
