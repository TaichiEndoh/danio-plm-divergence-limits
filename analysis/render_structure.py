#!/usr/bin/env python3
"""
Render a B-factor-colored PDB to a PNG — a structural figure without PyMOL/ChimeraX.

Reads the Cα trace of a PDB and colors each residue by its B-factor column (which
analysis/alphafold_local_mapping.py fills with local divergence). Produces a shareable
3D image directly, so no desktop molecular viewer is needed.

Usage:
    python analysis/render_structure.py --pdb reports/ahr2_colored.pdb \\
        --out reports/figures/ahr2_structure_divergence.png \\
        --title "AHR2 ESMFold — colored by local divergence" --label "divergence x1000"
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_ca(pdb: Path):
    xs, ys, zs, b, res = [], [], [], [], []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 66:
            try:
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
                b.append(float(line[60:66])); res.append(int(line[22:26]))
            except ValueError:
                pass
    return (np.array(xs), np.array(ys), np.array(zs), np.array(b), np.array(res))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="structure colored by B-factor")
    ap.add_argument("--label", default="B-factor")
    ap.add_argument("--cmap", default="viridis")
    ap.add_argument("--highlight", default="", help="comma-separated residue numbers to ring (e.g. pocket).")
    args = ap.parse_args()

    x, y, z, b, res = read_ca(Path(args.pdb))
    if len(x) == 0:
        raise SystemExit("no CA atoms found in PDB")

    fig = plt.figure(figsize=(9, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, color="0.6", lw=1.0, alpha=0.7, zorder=1)  # backbone trace
    p = ax.scatter(x, y, z, c=b, cmap=args.cmap, s=22, zorder=2, depthshade=False)
    if args.highlight:
        hl = {int(v) for v in args.highlight.split(",") if v.strip()}
        m = np.isin(res, list(hl))
        if m.any():
            ax.scatter(x[m], y[m], z[m], facecolors="none", edgecolors="red", s=70,
                       linewidths=1.2, zorder=3, label="highlighted")
            ax.legend(loc="upper left", frameon=False)
    cb = fig.colorbar(p, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label(args.label)
    ax.set_title(args.title, fontsize=12, weight="bold")
    ax.set_axis_off()
    try:
        ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
    except Exception:
        pass
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}  ({len(x)} residues; B-factor range {b.min():.2f}..{b.max():.2f})")


if __name__ == "__main__":
    main()
