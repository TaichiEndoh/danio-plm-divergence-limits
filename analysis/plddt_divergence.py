#!/usr/bin/env python3
"""
Quantify local embedding divergence vs ESMFold confidence (pLDDT) per residue.

Answers the reviewer's question for the structure figure: are the high-divergence
residues merely low-confidence / disordered (pLDDT noise), or structured functional
changes? Reads per-residue pLDDT from an ESMFold PDB (B-factor of CA, 0-1 -> x100) and
per-residue local divergence, correlates them (Spearman), stratifies residues by pLDDT
band, and reports where the divergence hotspots fall. Writes a scatter figure.

Usage:
    python analysis/plddt_divergence.py \\
        --pdb reports/ahr2_esmfold.pdb \\
        --divergence reports/ahr2_residue_divergence_mac.csv \\
        --out reports/figures/ahr2_plddt_vs_divergence.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OI = {"blue": "#0072B2", "orange": "#E69F00", "vermillion": "#D55E00", "grey": "#999999"}


def read_plddt(pdb: Path) -> pd.DataFrame:
    rows = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 66:
            try:
                rows.append((int(line[22:26]), float(line[60:66])))
            except ValueError:
                pass
    df = pd.DataFrame(rows, columns=["residue", "plddt"])
    if df["plddt"].max() <= 1.5:      # ESMFold transformers writes 0-1; normalize to 0-100
        df["plddt"] *= 100.0
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True, help="ESMFold PDB (pLDDT in B-factor).")
    ap.add_argument("--divergence", required=True, help="CSV: residue,divergence,...")
    ap.add_argument("--out", required=True)
    ap.add_argument("--protein", default="protein", help="protein name for the figure title.")
    ap.add_argument("--highlight", default="", help="comma-separated residues to ring (e.g. pocket).")
    ap.add_argument("--highlight-label", default="highlighted residues")
    args = ap.parse_args()

    highlight = {int(v) for v in args.highlight.split(",") if v.strip()}
    plddt = read_plddt(Path(args.pdb))
    div = pd.read_csv(args.divergence)[["residue", "divergence"]]
    df = div.merge(plddt, on="residue", how="inner").dropna()
    df["pocket"] = df["residue"].isin(highlight)

    rho, p = spearmanr(df["divergence"], df["plddt"])
    print(f"n residues: {len(df)}")
    print(f"Spearman(divergence, pLDDT): rho={rho:.3f}  p={p:.2e}")

    # stratify by pLDDT band
    bands = [("disordered (<50)", df["plddt"] < 50),
             ("intermediate (50-70)", (df["plddt"] >= 50) & (df["plddt"] < 70),),
             ("confident (>=70)", df["plddt"] >= 70)]
    print("\npLDDT band            n     mean_div   median_div")
    for name, m in bands:
        s = df.loc[m, "divergence"]
        print(f"  {name:20s} {len(s):4d}   {s.mean():.5f}   {s.median():.5f}")

    # where do the top-divergence hotspots sit?
    top = df.nlargest(int(np.ceil(0.10 * len(df))), "divergence")
    print(f"\ntop-10% divergent residues (n={len(top)}): "
          f"mean pLDDT={top['plddt'].mean():.1f}  "
          f"({(top['plddt']<50).mean()*100:.0f}% in disordered band)")
    print(f"rest: mean pLDDT={df.loc[~df.index.isin(top.index),'plddt'].mean():.1f}")

    # figure
    fig, ax = plt.subplots(figsize=(8.2, 6))
    nonp = df[~df["pocket"]]
    ax.scatter(nonp["plddt"], nonp["divergence"], s=16, c=OI["blue"], alpha=0.55,
               edgecolors="none", label="residue")
    pk = df[df["pocket"]]
    if len(pk):
        ax.scatter(pk["plddt"], pk["divergence"], s=45, facecolors="none", edgecolors=OI["vermillion"],
                   linewidths=1.4, label=args.highlight_label)
    ax.axvline(50, color=OI["grey"], ls="--", lw=1)
    ax.text(50, ax.get_ylim()[1], " pLDDT 50", color=OI["grey"], fontsize=8, va="top")
    ax.set_xlabel("ESMFold pLDDT (per residue; higher = more confident/structured)")
    ax.set_ylabel("local ESM-2 divergence (±10 aa window)")
    ax.set_title(f"{args.protein}: local divergence vs structural confidence\n"
                 f"Spearman ρ = {rho:.2f} (p = {p:.1e}); n = {len(df)}", fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
