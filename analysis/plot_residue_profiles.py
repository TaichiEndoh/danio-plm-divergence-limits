#!/usr/bin/env python3
"""
Per-residue local divergence profiles: ESM-2 (8M, corrected windows) vs ESM-3 (1.4B),
for AHR2 and Kcnj13.

For Kcnj13 the ESM-3 profile is contaminated at the N-terminus by the four-residue offset
between the two RefSeq entries (Results §3.7); that region is shaded and labelled rather than
silently dropped, and the reported correlation is given both over the whole profile and
restricted to residues beyond the affected window.

Writes reports/figures/esm2_vs_esm3_divergence.png and esm2_ladder_to_esm3.png.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OI = {"blue": "#0072B2", "vermillion": "#D55E00", "grey": "#999999", "orange": "#E69F00",
      "green": "#009E73"}
FIG = Path("reports/figures")
FLANK = 10  # window half-width used for the ESM-2 profiles


def col(df: pd.DataFrame) -> str:
    return [c for c in df.columns if "diver" in c][0]


def load(path: str) -> pd.DataFrame:
    d = pd.read_csv(path)
    return d[["residue", col(d)]].rename(columns={col(d): "v"})


def panel(ax, e2, e3, name, show_artifact=False, restrict_above=None, marks=()):
    """Plot both profiles over their own ranges.

    ESM-2 has no values within `flank` of either terminus (those windows would be clipped),
    so its line starts later; ESM-3 has values there and, where the two entries are offset at
    the N-terminus, spikes. Both facts are shown rather than hidden.
    """
    m = e2.merge(e3, on="residue", suffixes=("_2", "_3"))
    ax.plot(e3.residue, e3.v, lw=0.9, color=OI["vermillion"], alpha=0.85, label="ESM-3 (1.4B)")
    ax.plot(e2.residue, e2.v, lw=0.9, color=OI["blue"], label="ESM-2 (8M), corrected windows")

    # y-limit set by the biology, not by the artifact
    top = max(m.v_2.max(), m.v_3.max()) * 1.35
    ax.set_ylim(0, top)

    rho_all = spearmanr(m.v_2, m.v_3).statistic
    stats = f"Spearman ρ = {rho_all:.3f}  (n = {len(m)}, residues {m.residue.min()}–{m.residue.max()})"
    if restrict_above:
        sub = m[m.residue > restrict_above]
        stats += f"\nρ = {spearmanr(sub.v_2, sub.v_3).statistic:.3f} for residues > {restrict_above}"

    if show_artifact:
        lo, hi = e3.residue.min(), e2.residue.min() - 1
        ax.axvspan(lo - 0.5, hi + 0.5, color=OI["grey"], alpha=0.25, lw=0)
        peak = e3.loc[e3.v.idxmax()]
        ax.annotate(
            f"ESM-3 reaches {peak.v:.3f} at residue {int(peak.residue)}\n"
            f"(off scale; {peak.v / m.v_2.max():.0f}× the corrected ESM-2 maximum).\n"
            f"Terminal-offset artifact — the two RefSeq entries differ\n"
            f"by four residues at the N-terminus. ESM-2 has no window\n"
            f"here by construction, so residues {lo}–{hi} are excluded.",
            xy=(int(peak.residue), top * 0.97), xytext=(0.30, 0.93),
            textcoords="axes fraction", fontsize=7.6, va="top", color="#333",
            arrowprops=dict(arrowstyle="->", color="#777", lw=1,
                            connectionstyle="arc3,rad=-0.2"),
            bbox=dict(fc="#fffbf0", ec="#b8860b", lw=.6, alpha=.95))

    for pos, lab in marks:
        ax.annotate(lab, xy=(pos, 0), xytext=(pos, top * 0.46), ha="center",
                    fontsize=8.5, color=OI["green"], weight="bold",
                    arrowprops=dict(arrowstyle="-", color=OI["green"], lw=1))

    ax.set_title(f"{name}: per-residue local divergence", fontsize=11, weight="bold")
    ax.set_xlabel("residue (D. rerio numbering)")
    ax.set_ylabel("local cosine distance")
    ax.text(0.995, 0.97, stats, transform=ax.transAxes, ha="right", va="top", fontsize=8.2,
            bbox=dict(fc="white", ec="#bbb", lw=.5, alpha=.92))
    ax.legend(frameon=False, fontsize=8.2, loc="upper left", bbox_to_anchor=(0.0, 1.0))
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11, 8.0))
    panel(axes[0], load("reports/ahr2_residue_divergence_fixed.csv"),
          load("reports/ahr2_residue_divergence_esm3.csv"), "AHR2")
    panel(axes[1], load("reports/kcnj13_residue_divergence_fixed.csv"),
          load("reports/kcnj13_residue_divergence_esm3.csv"), "Kcnj13",
          show_artifact=True, restrict_above=12,
          marks=((19, "Q19L"), (176, "D176G")))
    fig.tight_layout(h_pad=2.2)
    fig.savefig(FIG / "esm2_vs_esm3_divergence.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", FIG / "esm2_vs_esm3_divergence.png")

    # ladder agreement with ESM-3, per protein
    ladder = {"8M": "fixed", "35M": "esm2_t12_35M_UR50D",
              "150M": "esm2_t30_150M_UR50D", "650M": "esm2_t33_650M_UR50D"}
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    width, xs = 0.36, np.arange(len(ladder))
    for k, (prot, e3path) in enumerate({"AHR2": "reports/ahr2_residue_divergence_esm3.csv",
                                        "Kcnj13": "reports/kcnj13_residue_divergence_esm3.csv"}.items()):
        e3 = load(e3path)
        vals = []
        for size, suffix in ladder.items():
            p = f"reports/{prot.lower()}_residue_divergence_{suffix}.csv"
            if not Path(p).is_file():
                vals.append(np.nan); continue
            m = load(p).merge(e3, on="residue", suffixes=("_2", "_3"))
            vals.append(spearmanr(m.v_2, m.v_3).statistic)
        ax.bar(xs + (k - 0.5) * width, vals, width,
               color=[OI["blue"], OI["orange"]][k], label=prot)
        for x, v in zip(xs + (k - 0.5) * width, vals):
            if not np.isnan(v):
                ax.text(x, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xs); ax.set_xticklabels([f"ESM-2 {s}" for s in ladder])
    ax.set_ylabel("Spearman ρ against ESM-3 (1.4B)")
    ax.set_title("Residue-profile agreement with ESM-3 across the ESM-2 ladder",
                 fontsize=11, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "esm2_ladder_to_esm3.png", dpi=200, bbox_inches="tight")
    print("Wrote", FIG / "esm2_ladder_to_esm3.png")


if __name__ == "__main__":
    main()
