#!/usr/bin/env python3
"""
Select an ESM-3 re-embedding subset from the ESM-2 8M baseline — with a design that
separates the tail from an unbiased estimate of the whole distribution.

Why not top-N only? Selecting pairs by the outcome variable (ESM-2 distance) and then
reading ESM-3 distance off that same selection invites regression to the mean: the very
pairs picked *because* they scored high on one model will tend to score lower on another,
purely as a selection artifact. To let the analysis distinguish a real cross-model signal
from that artifact, the subset is a labeled UNION of four selection modes, and every pair
carries the mode(s) that put it in:

  1. ``stratified``  — decile-stratified random sample of the 8M distance distribution
                       (``--per-decile`` pairs from each of 10 deciles). An unbiased,
                       distribution-spanning backbone; NOT conditioned on being extreme.
  2. ``top``         — the ``--top-n`` most divergent pairs (tail analysis).
  3. ``random``      — a plain random sample of ``--random-n`` pairs (a control whose
                       ESM-3 values are not selected on 8M distance at all).
  4. ``target``      — the curated proteins (AHR subtypes, Kcnj13 + pigment comparators),
                       always included regardless of rank.

Sampling is seeded (``--seed``) for reproducibility. Input defaults to the de-duplicated
baseline (one row per unique pair); see analysis/build_baseline_dedup.py.

Usage:
    python analysis/select_high_divergence_subset.py \\
        --baseline data/reference/esm2_8m_pair_distances_dedup.csv \\
        --ortholog-table data/reference/ortholog_pairs_rbh.csv \\
        --per-decile 150 --top-n 1000 --random-n 1000 --seed 0 \\
        --out reports/esm3_subset_pairs.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Curated targets keyed by zebrafish (D. rerio) RefSeq id. Sources:
# analysis/ported/ahr_pocket_local_plm.py (AHR pocket subtypes) and
# analysis/ported/kcnj13_local_plm.py (Kcnj13 + pigment comparators). NP_001417877
# (AHR1a full-length) is omitted: not in the RBH table (handled via override FASTA before).
CURATED_TARGETS = {
    "NP_571339": "AHR2",
    "NP_571337": "AHR1b",
    "NP_571338": "AHR1a-like",
    "NP_001039014": "Kcnj13",
    "NP_001038288": "Cx39.4",
    "NP_001030160": "Cx41.8",
    "NP_001265751": "Igsf11",
}

MODES = ["stratified", "top", "random", "target"]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baseline", default="data/reference/esm2_8m_pair_distances_dedup.csv",
                    help="per-pair ESM-2 distances CSV (default: the dedup file).")
    ap.add_argument("--ortholog-table", default="data/reference/ortholog_pairs_rbh.csv",
                    help="RBH pairs CSV, used only to resolve curated targets' partner ids.")
    ap.add_argument("--out", required=True, help="output subset pairs CSV.")
    ap.add_argument("--per-decile", type=int, default=150,
                    help="stratified pairs per decile (10 deciles; default %(default)s -> ~1500).")
    ap.add_argument("--top-n", type=int, default=1000, help="tail: N most divergent pairs.")
    ap.add_argument("--random-n", type=int, default=1000, help="random control sample size.")
    ap.add_argument("--seed", type=int, default=0, help="RNG seed for reproducible sampling.")
    ap.add_argument("--no-targets", action="store_true", help="do not force-include curated targets.")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    base = pd.read_csv(args.baseline)
    if not {"zebrafish_id", "aescallii_id", "plm_cosine_distance"} <= set(base.columns):
        raise SystemExit("Baseline CSV must have zebrafish_id, aescallii_id, plm_cosine_distance columns.")
    base["plm_cosine_distance"] = pd.to_numeric(base["plm_cosine_distance"], errors="coerce")
    n_dup = base.duplicated(subset=["zebrafish_id", "aescallii_id"]).sum()
    if n_dup:
        base = base.drop_duplicates(subset=["zebrafish_id", "aescallii_id"], keep="first")
        print(f"  note: input had {n_dup} duplicate pair rows; collapsed to {len(base)} unique pairs")
    df = base.dropna(subset=["plm_cosine_distance"]).reset_index(drop=True)

    # membership flags per selection mode
    flags = {m: pd.Series(False, index=df.index) for m in MODES}

    # 1. decile-stratified sample
    if args.per_decile > 0:
        # qcut into 10 groups by distance; duplicates="drop" guards against tied edges in
        # the heavily-zero-inflated low tail (fewer than 10 bins is fine).
        deciles = pd.qcut(df["plm_cosine_distance"], 10, labels=False, duplicates="drop")
        for d in sorted(pd.unique(deciles.dropna())):
            idx = df.index[deciles == d].to_numpy()
            k = min(args.per_decile, len(idx))
            pick = rng.choice(idx, size=k, replace=False)
            flags["stratified"].loc[pick] = True

    # 2. tail
    if args.top_n > 0:
        top_idx = df.nlargest(args.top_n, "plm_cosine_distance").index
        flags["top"].loc[top_idx] = True

    # 3. random control
    if args.random_n > 0:
        k = min(args.random_n, len(df))
        pick = rng.choice(df.index.to_numpy(), size=k, replace=False)
        flags["random"].loc[pick] = True

    # 4. curated targets (resolve partner id via the ortholog table so they always appear)
    target_rows = []
    if not args.no_targets:
        flags["target"].loc[df["zebrafish_id"].isin(CURATED_TARGETS)] = True
        present = set(df.loc[flags["target"], "zebrafish_id"])
        missing = [t for t in CURATED_TARGETS if t not in present]
        if missing:
            ortho = pd.read_csv(args.ortholog_table).rename(
                columns={"Danio_rerio": "zebrafish_id", "Danio rerio": "zebrafish_id",
                         "Danio aesculapii": "aescallii_id"})
            om = ortho[ortho["zebrafish_id"].isin(missing)][["zebrafish_id", "aescallii_id"]]
            for z, a in om.itertuples(index=False):
                target_rows.append({"zebrafish_id": z, "aescallii_id": a, "plm_cosine_distance": np.nan})
            still = [t for t in missing if t not in set(om["zebrafish_id"])]
            for t in still:
                print(f"  WARNING: curated target {CURATED_TARGETS[t]} ({t}) not found in baseline or ortholog table")

    keep = df[flags["stratified"] | flags["top"] | flags["random"] | flags["target"]].copy()
    for m in MODES:
        keep[f"in_{m}"] = flags[m].loc[keep.index].to_numpy()

    # append any curated targets that were absent from the baseline
    if target_rows:
        extra = pd.DataFrame(target_rows)
        for m in MODES:
            extra[f"in_{m}"] = (m == "target")
        keep = pd.concat([keep, extra], ignore_index=True)

    keep["target_name"] = keep["zebrafish_id"].map(CURATED_TARGETS).fillna("")
    keep["selection_reasons"] = [";".join(m for m in MODES if row[f"in_{m}"])
                                 for _, row in keep.iterrows()]
    keep = keep.sort_values("plm_cosine_distance", ascending=False, na_position="last").reset_index(drop=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["zebrafish_id", "aescallii_id", "plm_cosine_distance",
            "selection_reasons", "target_name"] + [f"in_{m}" for m in MODES]
    keep[cols].to_csv(out, index=False)

    print(f"Selected {len(keep)} unique pairs -> {out}  (seed={args.seed})")
    for m in MODES:
        print(f"  in_{m:10s}: {int(keep[f'in_{m}'].sum())}")
    d = keep["plm_cosine_distance"].dropna()
    print(f"  ESM-2 distance in subset: min={d.min():.5f} median={d.median():.5f} max={d.max():.5f}")


if __name__ == "__main__":
    main()
