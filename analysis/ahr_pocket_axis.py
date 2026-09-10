#!/usr/bin/env python3
"""
AHR ligand-pocket divergence: ORTHOLOG (cross-species) vs PARALOG (cross-subtype).

Resolves the axis confusion flagged for the paper (§2.2): the previous study's "pocket
divergence" was a subtype (paralog) comparison, whereas the direct D. rerio vs D.
aesculapii ortholog comparison shows the pocket is conserved. This script quantifies both
axes at the AHR ligand-pocket residues (Takeda) and plots the contrast.

Inputs are per-residue divergence CSVs (from analysis/local_residue_divergence.py), each
numbered on its zebrafish_id; we pull that protein's pocket residues.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

POCKET = {
    "AHR2":  [281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356, 370, 378, 380, 381],
    "AHR1a": [278, 280, 286, 288, 300, 304, 308, 311, 319, 321, 332, 337, 345, 353, 367, 375, 377, 378],
}
OI = {"blue": "#0072B2", "vermillion": "#D55E00"}

# (label, csv, pocket_key, axis)
COMPARISONS = [
    ("AHR2\n(rerio vs aesc)",  "reports/ahr2_residue_divergence_fixed.csv",       "AHR2",  "ortholog"),
    ("AHR1a\n(rerio vs aesc)", "reports/ahr1a_ortholog_divergence.csv",         "AHR1a", "ortholog"),
    ("AHR2 vs AHR1a\n(paralog)", "reports/ahr2_vs_ahr1a_paralog_divergence.csv", "AHR2",  "paralog"),
]


def main() -> None:
    rows = []
    for label, csv, pk, axis in COMPARISONS:
        df = pd.read_csv(csv)
        pocket = df[df["residue"].isin(POCKET[pk])]["divergence"]
        whole = df["divergence"]
        rows.append({"label": label, "axis": axis,
                     "pocket_mean": pocket.mean(), "pocket_median": pocket.median(),
                     "whole_mean": whole.mean(), "n_pocket": len(pocket)})
    res = pd.DataFrame(rows)
    print(res.to_string(index=False))

    fig, ax = plt.subplots(figsize=(7.5, 5))
    colors = [OI["blue"] if a == "ortholog" else OI["vermillion"] for a in res["axis"]]
    x = np.arange(len(res))
    ax.bar(x, res["pocket_mean"], color=colors, width=0.62)
    for i, v in enumerate(res["pocket_mean"]):
        ax.text(i, v, f" {v:.4f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(res["label"], fontsize=9)
    ax.set_ylabel("mean pocket local ESM-2 divergence")
    ax.set_title("AHR ligand-pocket divergence:\nconserved across species (ortholog) — divergent across subtypes (paralog)",
                 fontsize=11.5, weight="bold")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=OI["blue"], label="ortholog (cross-species)"),
                       Patch(color=OI["vermillion"], label="paralog (cross-subtype)")],
              frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    out = Path("reports/figures/ahr_pocket_ortholog_vs_paralog.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
