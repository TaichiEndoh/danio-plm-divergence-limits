#!/usr/bin/env python3
"""
Map local embedding divergence onto a predicted structure.

Idea: for a target protein (e.g., AHR2, Kcnj13), take per-residue (or per-window) local
divergence values and write them into the B-factor column of a predicted PDB, so the
structure can be colored by divergence in PyMOL/ChimeraX. This is the "structural
validation" bridge: do high-divergence residues cluster in pockets / functional regions?

Structure source is flexible: pass a local `--pdb` (e.g. an ESMFold prediction from
analysis/esmfold_predict.py — the project's chosen route, since full-length AHR2 has no
AlphaFold DB model) or `--uniprot` to fetch from the AlphaFold DB.

Divergence values are small (~0-0.07 cosine distance) while the PDB B-factor column keeps
only 2 decimals, so raw values collapse to a few colors. Use `--scale` (e.g. 1000) to
spread them across a usable range before writing. A `--demo` mode writes synthetic values
so the flow can be tested offline.

Usage:
    # color an ESMFold prediction by a divergence table (scale for visible contrast):
    python analysis/alphafold_local_mapping.py --pdb reports/ahr2_esmfold.pdb \\
        --divergence reports/ahr2_residue_divergence_esm2_8m.csv --scale 1000 \\
        --out reports/ahr2_colored.pdb

    # or fetch from the AlphaFold DB by UniProt accession:
    python analysis/alphafold_local_mapping.py --uniprot Q1LWM8 \\
        --divergence div.csv --scale 1000 --out reports/af_colored.pdb
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import numpy as np

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
AFDB_FILE = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{ver}.pdb"


def fetch_alphafold_pdb(uniprot: str, dest: Path) -> Path:
    """Download an AlphaFold DB predicted structure by UniProt accession.

    The AF DB bumps its model version over time (v4 -> v6 as of 2026), so the file URL
    is resolved from the prediction API rather than hard-coded; if the API is unreachable,
    fall back to trying the newest known version numbers.
    """
    import requests  # local import so --demo works without the dependency installed
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = None
    try:
        meta = requests.get(AFDB_API.format(acc=uniprot), timeout=60).json()
        if meta:
            url = meta[0].get("pdbUrl")
    except Exception:
        url = None
    if url:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        dest.write_bytes(r.content)
        return dest
    for ver in (6, 5, 4):  # fallback if the API is unavailable
        r = requests.get(AFDB_FILE.format(acc=uniprot, ver=ver), timeout=60)
        if r.ok:
            dest.write_bytes(r.content)
            return dest
    raise SystemExit(f"No AlphaFold model found for {uniprot} (API and v6/v5/v4 file URLs all failed).")


def load_divergence(path: Path) -> Dict[int, float]:
    """CSV with columns: residue (1-based int), divergence (float)."""
    import pandas as pd
    df = pd.read_csv(path)
    return {int(r.residue): float(r.divergence) for r in df.itertuples()}


def write_bfactor(pdb_in: Path, pdb_out: Path, per_residue: Dict[int, float]) -> None:
    """Rewrite the B-factor column (cols 61-66) of ATOM records with per-residue values.

    Kept dependency-free (plain PDB text munging) so it runs anywhere. For CIF or complex
    structures, switch to Biopython's PDBParser/MMCIFParser.
    """
    default = 0.0
    lines_out = []
    for line in pdb_in.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                resseq = int(line[22:26])
            except ValueError:
                lines_out.append(line)
                continue
            val = per_residue.get(resseq, default)
            lines_out.append(f"{line[:60]}{val:6.2f}{line[66:]}")
        else:
            lines_out.append(line)
    pdb_out.parent.mkdir(parents=True, exist_ok=True)
    pdb_out.write_text("\n".join(lines_out) + "\n")


def demo_divergence(pdb: Path) -> Dict[int, float]:
    """Synthetic per-residue divergence (seeded) for offline testing."""
    resids = set()
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and len(line) >= 26:
            try:
                resids.add(int(line[22:26]))
            except ValueError:
                pass
    rng = np.random.default_rng(0)
    return {r: float(abs(rng.standard_normal()) * 0.05) for r in sorted(resids)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uniprot", help="UniProt accession to fetch from AlphaFold DB")
    ap.add_argument("--pdb", help="local PDB path (instead of fetching)")
    ap.add_argument("--divergence", help="CSV: residue,divergence")
    ap.add_argument("--demo", action="store_true", help="use synthetic divergence (offline)")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply divergence before writing (e.g. 1000; B-factor keeps 2 decimals).")
    ap.add_argument("--out", required=True, help="output PDB with divergence in B-factor column")
    args = ap.parse_args()

    if args.pdb:
        pdb_in = Path(args.pdb)
    elif args.uniprot:
        pdb_in = fetch_alphafold_pdb(args.uniprot, Path("data/structures") / f"{args.uniprot}.pdb")
    else:
        ap.error("provide --pdb or --uniprot")

    if args.demo:
        per_residue = demo_divergence(pdb_in)
    elif args.divergence:
        per_residue = load_divergence(Path(args.divergence))
    else:
        ap.error("provide --divergence or --demo")

    if args.scale != 1.0:
        per_residue = {r: v * args.scale for r, v in per_residue.items()}

    write_bfactor(pdb_in, Path(args.out), per_residue)
    if per_residue:
        lo, hi = min(per_residue.values()), max(per_residue.values())
        print(f"Wrote {args.out}  (B-factor holds divergence*{args.scale:g}, range {lo:.2f}..{hi:.2f})")
    else:
        print(f"Wrote {args.out}")
    print("  open in PyMOL/ChimeraX and color by B-factor — PyMOL: spectrum b, blue_white_red")


if __name__ == "__main__":
    main()
