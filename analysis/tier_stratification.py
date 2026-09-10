#!/usr/bin/env python3
"""
Three-tier stratification of local divergence (Prof. Endo's §4):
  structured functional core  /  structured periphery  /  disordered–low-confidence.

Tiers are defined per residue from ESMFold pLDDT plus a functional-core residue set:
  - core        : annotated functional-core residues (AHR2 ligand pocket)
  - structured  : pLDDT >= 70, not core (candidate modulatory periphery)
  - intermediate: 50 <= pLDDT < 70
  - disordered  : pLDDT < 50 (evolvable regulatory periphery candidate)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OI = {"blue": "#0072B2", "green": "#009E73", "orange": "#E69F00", "vermillion": "#D55E00"}
AHR2_POCKET = [281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356, 370, 378, 380, 381]


def read_plddt(pdb):
    d = {}
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 66:
            try:
                d[int(line[22:26])] = float(line[60:66])
            except ValueError:
                pass
    mx = max(d.values()) if d else 1.0
    f = 100.0 if mx <= 1.5 else 1.0
    return {r: v * f for r, v in d.items()}


def tiers(df, pocket):
    def classify(row):
        if row["residue"] in pocket:
            return "core"
        if row["plddt"] >= 70:
            return "structured"
        if row["plddt"] >= 50:
            return "intermediate"
        return "disordered"
    df = df.copy()
    df["tier"] = df.apply(classify, axis=1)
    return df


def main() -> None:
    jobs = [
        ("AHR2", "reports/ahr2_residue_divergence_fixed.csv", "reports/ahr2_esmfold.pdb", AHR2_POCKET),
        ("Kcnj13", "reports/kcnj13_residue_divergence_fixed.csv", "reports/kcnj13_esmfold.pdb", []),
    ]
    order = ["core", "structured", "intermediate", "disordered"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    for ax, (name, dcsv, pdb, pocket) in zip(axes, jobs):
        div = pd.read_csv(dcsv)[["residue", "divergence"]]
        pl = pd.DataFrame(list(read_plddt(pdb).items()), columns=["residue", "plddt"])
        df = tiers(div.merge(pl, on="residue"), set(pocket))
        print(f"\n== {name} ==")
        rows = []
        for t in order:
            sub = df[df["tier"] == t]
            if len(sub):
                rows.append((t, len(sub), sub["divergence"].mean(), sub["divergence"].median()))
                print(f"  {t:12s} n={len(sub):4d}  mean_div={sub['divergence'].mean():.5f}  median={sub['divergence'].median():.5f}")
        r = pd.DataFrame(rows, columns=["tier", "n", "mean", "median"])
        colors = {"core": OI["blue"], "structured": OI["green"], "intermediate": OI["orange"], "disordered": OI["vermillion"]}
        ax.bar(range(len(r)), r["mean"], color=[colors[t] for t in r["tier"]], width=0.62)
        for i, (v, n) in enumerate(zip(r["mean"], r["n"])):
            ax.text(i, v, f"{v:.4f}\n(n={n})", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(range(len(r))); ax.set_xticklabels(r["tier"], fontsize=9)
        ax.set_ylabel("mean local divergence")
        ax.set_title(name, fontsize=12, weight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Divergence by structural tier: conserved core vs evolvable periphery",
                 fontsize=13, weight="bold")
    fig.tight_layout()
    out = Path("reports/figures/tier_stratification.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
