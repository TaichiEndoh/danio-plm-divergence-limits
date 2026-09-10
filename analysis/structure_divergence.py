#!/usr/bin/env python3
"""
Per-residue structural divergence between two species' ESMFold structures, correlated
with local embedding divergence:  D_PLM(i)  <->  D_structure(i).

For an ortholog pair, superpose the two predicted structures on their shared confident
core (both pLDDT > --plddt-core) via Kabsch, then measure the per-(zebrafish-)residue
Cα displacement as local structural divergence, and correlate it with the per-residue
local ESM divergence. Tests whether PLM divergence tracks predicted structural change.

Residue correspondence (zebrafish i <-> aesculapii j) comes from the divergence CSV's
`aescallii_residue` column (the NW alignment used to compute local divergence).

Usage:
    python analysis/structure_divergence.py \\
        --pdb-a reports/kcnj13_esmfold.pdb --pdb-b reports/kcnj13_aesc_esmfold.pdb \\
        --divergence reports/kcnj13_residue_divergence.csv --protein Kcnj13 \\
        --out reports/figures/kcnj13_plm_vs_structure.png
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


def read_ca(pdb: Path):
    d = {}
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 66:
            try:
                r = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                d[r] = (xyz, float(line[60:66]))
            except ValueError:
                pass
    mx = max((b for _, b in d.values()), default=1.0)
    f = 100.0 if mx <= 1.5 else 1.0            # normalize pLDDT 0-1 -> 0-100
    return {r: (xyz, b * f) for r, (xyz, b) in d.items()}


def kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rotation R (3x3) minimising |P@R.T - Q|; P,Q are centered (N,3)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb-a", required=True, help="zebrafish (D. rerio) ESMFold PDB")
    ap.add_argument("--pdb-b", required=True, help="aesculapii ESMFold PDB")
    ap.add_argument("--divergence", required=True, help="CSV: residue, divergence, aescallii_residue")
    ap.add_argument("--plddt-core", type=float, default=70.0, help="superpose on residues with both pLDDT > this")
    ap.add_argument("--protein", default="protein")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    A, B = read_ca(Path(args.pdb_a)), read_ca(Path(args.pdb_b))
    div = pd.read_csv(args.divergence)

    rows = []
    for _, row in div.iterrows():
        ri, aj = int(row["residue"]), int(row["aescallii_residue"])
        if aj == 0 or ri not in A or aj not in B:
            continue
        rows.append((ri, aj, float(row["divergence"]), A[ri][1], B[aj][1]))
    pairs = pd.DataFrame(rows, columns=["ri", "aj", "divergence", "plddt_a", "plddt_b"])
    if len(pairs) < 10:
        raise SystemExit("too few matched residues")

    Pa = np.array([A[ri][0] for ri in pairs["ri"]])
    Pb = np.array([B[aj][0] for aj in pairs["aj"]])
    core = ((pairs["plddt_a"] > args.plddt_core) & (pairs["plddt_b"] > args.plddt_core)).to_numpy()
    if core.sum() < 10:
        core = np.ones(len(pairs), bool)      # fallback: use all

    mu_a, mu_b = Pa[core].mean(0), Pb[core].mean(0)
    R = kabsch(Pb[core] - mu_b, Pa[core] - mu_a)
    Pb_sup = (Pb - mu_b) @ R.T + mu_a
    pairs["ca_dev"] = np.linalg.norm(Pa - Pb_sup, axis=1)   # Å

    rho, p = spearmanr(pairs["divergence"], pairs["ca_dev"])
    if core.sum() > 10:
        crho, cp = spearmanr(pairs.loc[core, "divergence"], pairs.loc[core, "ca_dev"])
    else:
        crho, cp = float("nan"), float("nan")
    core_rmsd = float(np.sqrt((pairs.loc[core, "ca_dev"] ** 2).mean()))
    print(f"n matched residues: {len(pairs)}  core(both pLDDT>{args.plddt_core}): {int(core.sum())}")
    print(f"core superposition RMSD: {core_rmsd:.2f} Å")
    print(f"Spearman(D_PLM, Ca-dev)  all : rho={rho:.3f}  p={p:.2e}")
    print(f"Spearman(D_PLM, Ca-dev)  core: rho={crho:.3f}  p={cp:.2e}")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(pairs.loc[~core, "ca_dev"], pairs.loc[~core, "divergence"], s=14,
               c="#999999", alpha=0.5, edgecolors="none", label="low-confidence residue")
    ax.scatter(pairs.loc[core, "ca_dev"], pairs.loc[core, "divergence"], s=16,
               c="#0072B2", alpha=0.6, edgecolors="none", label=f"core (pLDDT>{int(args.plddt_core)})")
    ax.set_xlabel("per-residue Cα deviation between species (Å)")
    ax.set_ylabel("local ESM divergence (±10 aa window)")
    ax.set_title(f"{args.protein}: PLM divergence vs predicted structural deviation\n"
                 f"Spearman ρ = {rho:.2f} (p={p:.1e}); core ρ = {crho:.2f}; n = {len(pairs)}",
                 fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    pairs.to_csv(out.with_suffix(".csv"), index=False)
    print(f"Wrote {out}  and  {out.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
