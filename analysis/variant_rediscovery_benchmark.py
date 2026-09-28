#!/usr/bin/env python3
"""
Stage-2 "rediscovery" benchmark: can ESM-2 rank a KNOWN functional variant above the bulk
of possible substitutions in the same protein? And does that improve with model scale?

For each control protein we build the full single-substitution set (L x 19) with the ESM-2
variant score (log P(mut) - log P(wt); see analysis/esm2_variant_score.py) and report where
each known pathogenic/functional variant falls: its percentile among all substitutions
(lower = more impactful = better rediscovery). A random substitution sits near 50% by
construction, so a known variant landing in the low single digits is the signal.

Runs across the ESM-2 ladder (8M/35M/150M/650M) so the result doubles as a scaling readout:
the whole paper asks whether small models already capture what large ones do.

wt-marginal scoring (one forward pass per sequence) is the default so the full ladder is cheap;
masked-marginal (one pass per position) is available but L x slower.

Usage:
    python analysis/variant_rediscovery_benchmark.py \\
        --fasta data/reference/sequences/positive_controls.fasta \\
        --controls "HBB_HUMAN_mature:E6V;RHO_BOVIN:E122Q,A292S" \\
        --models esm2_t6_8M_UR50D,esm2_t12_35M_UR50D,esm2_t30_150M_UR50D,esm2_t33_650M_UR50D \\
        --out reports/variant_rediscovery_benchmark.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur, buf = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0].split(".")[0], []
        elif cur:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def parse_controls(spec: str) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for part in spec.split(";"):
        part = part.strip()
        if not part:
            continue
        pid, variants = part.split(":")
        out[pid.strip()] = [v.strip().upper() for v in variants.split(",") if v.strip()]
    return out


def parse_variant(v: str) -> Tuple[str, int, str]:
    return v[0], int(v[1:-1]), v[-1]


def dms_scores(seq: str, model_name: str, scoring: str, device: str) -> pd.DataFrame:
    """Full single-substitution table with ESM-2 variant scores."""
    import torch
    import esm
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    bc = alphabet.get_batch_converter()
    _, _, toks = bc([("wt", seq)])
    toks = toks.to(device)
    aa_idx = {a: alphabet.get_idx(a) for a in AA}
    L = len(seq)

    lp_by_pos: Dict[int, np.ndarray] = {}
    with torch.no_grad():
        if scoring == "wt-marginal":
            lp = torch.log_softmax(model(toks)["logits"], dim=-1)[0].cpu().numpy()
            for p in range(1, L + 1):
                lp_by_pos[p] = lp[p]
        else:
            for p in range(1, L + 1):
                t = toks.clone()
                t[0, p] = alphabet.mask_idx
                lp_by_pos[p] = torch.log_softmax(model(t)["logits"], dim=-1)[0, p].cpu().numpy()

    rows = []
    for p in range(1, L + 1):
        wt = seq[p - 1]
        if wt not in aa_idx:
            continue
        for mut in AA:
            if mut == wt:
                continue
            rows.append({"variant": f"{wt}{p}{mut}", "position": p, "wt": wt, "mut": mut,
                         "score": float(lp_by_pos[p][aa_idx[mut]] - lp_by_pos[p][aa_idx[wt]])})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--controls", required=True, help='e.g. "HBB_HUMAN_mature:E6V;RHO_BOVIN:E122Q,A292S"')
    ap.add_argument("--models", default="esm2_t6_8M_UR50D,esm2_t12_35M_UR50D,esm2_t30_150M_UR50D,esm2_t33_650M_UR50D")
    ap.add_argument("--scoring", choices=["wt-marginal", "masked-marginal"], default="wt-marginal")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(Path(args.fasta))
    controls = parse_controls(args.controls)
    models = [m.strip() for m in args.models.split(",") if m.strip()]

    results = []
    for model_name in models:
        for pid, variants in controls.items():
            seq = seqs.get(pid)
            if not seq:
                print(f"  WARNING: {pid} not in fasta; skipping")
                continue
            df = dms_scores(seq, model_name, args.scoring, args.device)
            n = len(df)
            for v in variants:
                wt, pos, mut = parse_variant(v)
                if pos < 1 or pos > len(seq) or seq[pos - 1] != wt:
                    print(f"  WARNING: {pid} {v}: wt {wt} != residue at {pos}; skipping")
                    continue
                row = df[df["variant"] == v]
                score = float(row["score"].iloc[0])
                pct = 100.0 * (df["score"] <= score).sum() / n
                results.append({"protein": pid, "model": model_name, "scoring": args.scoring,
                                "variant": v, "known_score": round(score, 4),
                                "n_substitutions": n, "percentile_most_impactful": round(pct, 2)})
                print(f"  {model_name:22s} {pid:18s} {v:7s} score={score:7.3f} "
                      f"-> top {pct:5.2f}% of {n}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(results).to_csv(out, index=False)
    print(f"\nWrote {len(results)} rows -> {out}")


if __name__ == "__main__":
    main()
