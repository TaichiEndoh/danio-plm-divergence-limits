#!/usr/bin/env python3
"""
Fig 2 / Fig 3: subset-scale ESM-2 vs ESM-3 per-pair divergence.

Fig 2 — scatter of ESM-2 (8M) vs ESM-3 (1.4B) per-pair cosine distance (+ Spearman).
Fig 3 — candidate-prioritization reproducibility: overlap of the top-1% / top-5% most
        divergent pairs between each ESM-2 ladder size and ESM-3 (1.4B).

The paper's point (§2.1): absolute distances need NOT match across embedding spaces, but
divergence RANKING / hotspot prioritization is preserved.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00"}


def load(path, col):
    d = pd.read_csv(path)
    d["z"] = d["zebrafish_id"].astype(str).str.split(".").str[0]
    d["a"] = d["aescallii_id"].astype(str).str.split(".").str[0]
    return d[["z", "a", col]].rename(columns={col: "d"})


def topk_overlap(df, c1, c2, frac):
    k = max(1, int(round(frac * len(df))))
    s1 = set(map(tuple, df.nlargest(k, c1)[["z", "a"]].values))
    s2 = set(map(tuple, df.nlargest(k, c2)[["z", "a"]].values))
    return len(s1 & s2) / k, k


def main() -> None:
    esm3 = load("reports/esm3_subset_esm3_1p4b.csv", "esm3_cosine_distance").rename(columns={"d": "esm3"})
    ladder_all = {"8M": "reports/esm3_subset_esm2_8m.csv",
                  "35M": "reports/esm3_subset_esm2_35m.csv",
                  "150M": "reports/esm3_subset_esm2_150m.csv",
                  "650M": "reports/esm3_subset_esm2_650m.csv"}
    ladder = {k: v for k, v in ladder_all.items() if Path(v).is_file()}
    m = esm3.copy()
    for name, p in ladder.items():
        m = m.merge(load(p, "esm2_cosine_distance").rename(columns={"d": name}), on=["z", "a"])
    print(f"merged pairs: {len(m)}")

    # ---- Fig 2: ESM-2 8M vs ESM-3 scatter ----
    rho, pval = spearmanr(m["8M"], m["esm3"])
    fig, ax = plt.subplots(figsize=(7, 6.2))
    x = np.clip(m["8M"], 1e-5, None); y = np.clip(m["esm3"], 1e-5, None)
    ax.scatter(x, y, s=10, c=OI["blue"], alpha=0.4, edgecolors="none")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("ESM-2 (8M) per-pair cosine distance")
    ax.set_ylabel("ESM-3 (1.4B) per-pair cosine distance")
    ax.set_title(f"Subset per-pair divergence: ESM-2 (8M) vs ESM-3 (1.4B)\n"
                 f"Spearman ρ = {rho:.2f} (p = {pval:.1e}); n = {len(m)}", fontsize=11.5, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig("reports/figures/fig2_esm2_vs_esm3_subset.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ---- Fig 3: top-1% / top-5% overlap with ESM-3 across the ladder ----
    print("\nSpearman(size vs ESM-3) and top-k overlap:")
    sizes = list(ladder.keys())
    ov1 = []; ov5 = []; rhos = []
    for s in sizes:
        r, _ = spearmanr(m[s], m["esm3"]); rhos.append(r)
        o1, k1 = topk_overlap(m, s, "esm3", 0.01)
        o5, k5 = topk_overlap(m, s, "esm3", 0.05)
        ov1.append(o1 * 100); ov5.append(o5 * 100)
        print(f"  {s:5s} vs ESM-3: ρ={r:.3f}  top1%({k1})={o1*100:.0f}%  top5%({k5})={o5*100:.0f}%")

    fig, ax = plt.subplots(figsize=(7.5, 5))
    xi = np.arange(len(sizes)); w = 0.38
    ax.bar(xi - w/2, ov1, w, label="top-1% overlap", color=OI["vermillion"])
    ax.bar(xi + w/2, ov5, w, label="top-5% overlap", color=OI["orange"])
    for i, (a, b) in enumerate(zip(ov1, ov5)):
        ax.text(i - w/2, a, f"{a:.0f}%", ha="center", va="bottom", fontsize=9)
        ax.text(i + w/2, b, f"{b:.0f}%", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(xi); ax.set_xticklabels([f"ESM-2 {s}" for s in sizes])
    ax.set_ylabel("overlap of most-divergent pairs with ESM-3 (1.4B)")
    ax.set_ylim(0, 100)
    ax.set_title("Divergence-ranking reproducibility vs ESM-3 (1.4B)\n"
                 "(top-1% / top-5% most-divergent pairs)", fontsize=11.5, weight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig("reports/figures/fig3_ranking_overlap_esm3.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("\nWrote reports/figures/fig2_esm2_vs_esm3_subset.png and fig3_ranking_overlap_esm3.png")


if __name__ == "__main__":
    main()
