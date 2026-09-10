#!/usr/bin/env python3
"""
Fig 1 — concept diagram: the three-tier hierarchy of the study.

    Sequence conservation  ->  PLM representation divergence  ->  Structural grounding

Drawn with matplotlib primitives (no data); colorblind-safe Okabe-Ito palette.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "vermillion": "#D55E00", "grey": "#BBBBBB", "sky": "#56B4E9", "ink": "#222222"}


def stage_box(ax, x, w, title, color):
    ax.add_patch(FancyBboxPatch((x, 0.30), w, 0.52, boxstyle="round,pad=0.012,rounding_size=0.02",
                                linewidth=1.6, edgecolor=color, facecolor=color + "22"))
    ax.text(x + w / 2, 0.78, title, ha="center", va="center", fontsize=12.5, weight="bold", color=OI["ink"])


def main() -> None:
    fig, ax = plt.subplots(figsize=(12.5, 5))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.95, "Localized ortholog divergence, from sequence to structure",
            ha="center", va="center", fontsize=15, weight="bold", color=OI["ink"])
    ax.text(0.5, 0.89, "Danio rerio  vs  Danio aesculapii", ha="center", va="center",
            fontsize=11, style="italic", color=OI["grey"])

    # ---- Stage 1: Sequence ----
    x1, w = 0.03, 0.28
    stage_box(ax, x1, w, "1. Sequence", OI["blue"])
    # two aligned rows of residue cells; a few divergent columns highlighted
    rng = np.random.default_rng(0)
    ncol = 16
    div_cols = {4, 5, 11}
    for row, yy in [(0, 0.58), (1, 0.50)]:
        for c in range(ncol):
            col = OI["orange"] if c in div_cols else OI["grey"]
            ax.add_patch(plt.Rectangle((x1 + 0.02 + c * (w - 0.04) / ncol, yy),
                                       (w - 0.04) / ncol * 0.86, 0.06,
                                       facecolor=col, edgecolor="white", linewidth=0.4))
    ax.text(x1 + w / 2, 0.65, "globally conserved,\na few local differences",
            ha="center", va="bottom", fontsize=9.5, color=OI["ink"])
    ax.text(x1 + w / 2, 0.40, "ortholog pair\n(RBH; 50,921 pairs)", ha="center", va="center",
            fontsize=9, color=OI["ink"])

    # ---- Stage 2: PLM representation ----
    x2 = 0.36
    stage_box(ax, x2, w, "2. PLM divergence", OI["green"])
    # per-residue divergence profile with a peak
    xs = np.linspace(x2 + 0.03, x2 + w - 0.03, 100)
    prof = np.exp(-((np.linspace(-3, 3, 100) - (-1.2)) ** 2) / 0.25) * 0.11 \
        + np.exp(-((np.linspace(-3, 3, 100) - 1.8) ** 2) / 1.2) * 0.05
    ax.plot(xs, 0.50 + prof, color=OI["green"], lw=2)
    ax.hlines(0.50, x2 + 0.03, x2 + w - 0.03, color=OI["grey"], lw=0.8)
    ax.text(x2 + w / 2, 0.65, "localized divergence\n(per-residue / per-pair)",
            ha="center", va="bottom", fontsize=9.5, color=OI["ink"])
    ax.text(x2 + w / 2, 0.405, "ESM-2 (8M) → ESM-3 (1.4B)\nrank preserved (Spearman ρ = 0.88)",
            ha="center", va="center", fontsize=9, color=OI["ink"])

    # ---- Stage 3: Structure ----
    x3 = 0.69
    stage_box(ax, x3, w, "3. Structural grounding", OI["vermillion"])
    # simple backbone squiggle with a highlighted conserved pocket + a divergent tail
    t = np.linspace(0, 3 * np.pi, 200)
    bx = x3 + 0.04 + (np.cos(t) * 0.5 + 0.5) * (w - 0.08)
    by = 0.50 + np.sin(t * 1.3) * 0.05
    ax.plot(bx, by, color=OI["grey"], lw=2.2, solid_capstyle="round")
    ax.scatter(bx[80:110], by[80:110], s=10, color=OI["blue"], zorder=3)      # conserved core/pocket
    ax.scatter(bx[:14], by[:14], s=12, color=OI["orange"], zorder=3)          # divergent terminal
    ax.text(x3 + w / 2, 0.65, "map divergence onto\nESMFold / AlphaFold structure",
            ha="center", va="bottom", fontsize=9.5, color=OI["ink"])
    ax.text(x3 + w / 2, 0.40, "pocket conserved (core) ·\ndivergence at flexible / low-pLDDT regions",
            ha="center", va="center", fontsize=9, color=OI["ink"])

    # ---- arrows ----
    for xa, lbl in [(x1 + w + 0.005, "embed"), (x2 + w + 0.005, "fold + map")]:
        ax.annotate("", xy=(xa + 0.028, 0.56), xytext=(xa, 0.56),
                    arrowprops=dict(arrowstyle="-|>", lw=2.2, color=OI["ink"]))
        ax.text(xa + 0.014, 0.60, lbl, ha="center", va="bottom", fontsize=8.5, color=OI["ink"])

    out = Path("reports/figures/fig1_concept.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
