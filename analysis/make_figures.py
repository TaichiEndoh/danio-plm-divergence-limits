#!/usr/bin/env python3
"""
Generate the current result figures for sharing / the paper draft.

Outputs PNGs to reports/figures/:
  fig1_variant_rediscovery_scaling.png  — Stage-2 benchmark: known-variant percentile vs
                                           ESM-2 model size, split into conserved functional
                                           sites (recovered, improves with scale) vs
                                           non-conservation-driven controls (not recovered).
  fig2_baseline_divergence_distribution.png — ESM-2 8M per-pair cosine distance distribution
                                           (heavy right skew; signal in a small tail).
  fig3_ahr2_residue_divergence.png       — AHR2 per-residue local ESM-2 divergence along the
                                           sequence, with the ligand-pocket residues marked.

Colorblind-safe Okabe-Ito palette; recessive axes; direct labels.
NOTE: ESM-3 (1.4B) is not yet run (GPU + HF token pending) — these are the ESM-2 / ESMFold
results computed so far.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Okabe-Ito colorblind-safe palette
OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00",
      "sky": "#56B4E9", "reddish": "#CC79A7", "yellow": "#F0E442", "black": "#000000",
      "grey": "#999999"}

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200, "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.autolayout": True,
})

ROOT = Path(__file__).resolve().parent.parent
FIGDIR = ROOT / "reports" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

SIZE = {"esm2_t6_8M_UR50D": 8, "esm2_t12_35M_UR50D": 35,
        "esm2_t30_150M_UR50D": 150, "esm2_t33_650M_UR50D": 650}


def fig1_rediscovery() -> None:
    df = pd.read_csv(ROOT / "reports" / "variant_rediscovery_benchmark.csv")
    df = df[df["scoring"] == "wt-marginal"].copy()
    df["size"] = df["model"].map(SIZE)
    df = df.dropna(subset=["size"]).sort_values("size")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    groups = {
        "Conserved functional sites\n(recovered — improves with scale)":
            ([("HRAS_HUMAN", "G12V"), ("HRAS_HUMAN", "G12D"), ("TP53_HUMAN", "R175H")],
             [OI["blue"], OI["sky"], OI["green"]]),
        "Non-conservation-driven controls\n(not recovered at any scale)":
            ([("HBB_HUMAN_mature", "E6V"), ("RHO_BOVIN", "E122Q"), ("RHO_BOVIN", "A292S")],
             [OI["vermillion"], OI["orange"], OI["reddish"]]),
    }
    for ax, (title, (members, colors)) in zip(axes, groups.items()):
        end_labels = []  # (y, label, color) for de-collision
        for (prot, var), c in zip(members, colors):
            sub = df[(df["protein"] == prot) & (df["variant"] == var)].sort_values("size")
            if sub.empty:
                continue
            label = f"{prot.split('_')[0]} {var}"
            ax.plot(sub["size"], sub["percentile_most_impactful"], "-o", color=c,
                    lw=2, ms=7, label=label)
            end_labels.append([float(sub.iloc[-1]["percentile_most_impactful"]), label, c])
        # push apart end labels that are within 6 percentile units of each other
        end_labels.sort()
        for i in range(1, len(end_labels)):
            if end_labels[i][0] - end_labels[i - 1][0] < 6:
                end_labels[i][0] = end_labels[i - 1][0] + 6
        for y, label, c in end_labels:
            ax.annotate(label, (650, y), textcoords="offset points", xytext=(10, 0),
                        fontsize=9, color=c, va="center")
        ax.axhline(50, color=OI["grey"], ls="--", lw=1)
        ax.text(8, 52, "random ≈ 50%", color=OI["grey"], fontsize=8.5)
        ax.set_xscale("log")
        ax.set_xticks(list(SIZE.values()))
        ax.set_xticklabels([f"{s}M" for s in SIZE.values()])
        ax.set_xlabel("ESM-2 model size (parameters)")
        ax.set_title(title, fontsize=10.5)
        ax.set_xlim(6, 1500)
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("Known-variant percentile\n(lower = better rediscovery)")
    axes[0].invert_yaxis()
    fig.suptitle("Rediscovering known functional variants with the ESM-2 ladder",
                 fontsize=13, weight="bold")
    fig.savefig(FIGDIR / "fig1_variant_rediscovery_scaling.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig1")


def fig2_distribution() -> None:
    d = pd.read_csv(ROOT / "data/reference/esm2_8m_pair_distances_dedup.csv")
    x = pd.to_numeric(d["plm_cosine_distance"], errors="coerce").dropna().to_numpy()
    x = np.clip(x, 1e-6, None)  # for log axis
    med, p95, p99 = np.median(x), np.percentile(x, 95), np.percentile(x, 99)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    bins = np.logspace(np.log10(x.min()), np.log10(x.max()), 60)
    ax.hist(x, bins=bins, color=OI["blue"], alpha=0.85, edgecolor="white", linewidth=0.3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    for val, name, c in [(med, "median", OI["green"]), (p95, "95th", OI["orange"]),
                         (p99, "99th", OI["vermillion"])]:
        ax.axvline(val, color=c, lw=2)
        ax.text(val, ax.get_ylim()[1] * 0.5, f" {name}\n {val:.4f}", color=c, fontsize=9,
                rotation=0, va="top")
    ax.set_xlabel("ESM-2 8M per-pair cosine distance (log scale)")
    ax.set_ylabel("number of ortholog pairs (log)")
    ax.set_title(f"Ortholog divergence is concentrated in a small tail\n"
                 f"{len(x):,} D. rerio–D. aesculapii pairs; heavy right-skew",
                 fontsize=12, weight="bold")
    fig.savefig(FIGDIR / "fig2_baseline_divergence_distribution.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig2")


def fig3_ahr2() -> None:
    d = pd.read_csv(ROOT / "reports" / "ahr2_residue_divergence_fixed.csv")
    pocket = [281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356,
              370, 378, 380, 381]  # AHR2 ligand-pocket residues (Takeda)

    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.plot(d["residue"], d["divergence"], color=OI["blue"], lw=0.9)
    ax.axvspan(min(pocket), max(pocket), color=OI["orange"], alpha=0.15,
               label="ligand-pocket region (281–381)")
    pk = d.loc[d["divergence"].idxmax()]
    ax.annotate(f"peak (res {int(pk['residue'])})", (pk["residue"], pk["divergence"]),
                textcoords="offset points", xytext=(10, -2), fontsize=9, color=OI["vermillion"])
    ax.scatter([pk["residue"]], [pk["divergence"]], color=OI["vermillion"], zorder=5, s=30)
    pv = d[d["residue"].isin(pocket)]
    ax.scatter(pv["residue"], pv["divergence"], color=OI["orange"], s=22, zorder=4,
               label="pocket residues")
    ax.set_xlabel("AHR2 residue (D. rerio NP_571339, 1027 aa)")
    ax.set_ylabel("local ESM-2 divergence\n(±10 aa window cosine dist.)")
    ax.set_title("Per-residue local divergence along AHR2 (D. rerio vs D. aesculapii)",
                 fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_xlim(0, len(d) + 20)
    fig.savefig(FIGDIR / "fig3_ahr2_residue_divergence.png", bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3")


if __name__ == "__main__":
    fig1_rediscovery()
    fig2_distribution()
    fig3_ahr2()
    print(f"figures -> {FIGDIR}")
