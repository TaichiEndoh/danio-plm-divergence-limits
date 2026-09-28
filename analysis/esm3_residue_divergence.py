#!/usr/bin/env python3
"""
Per-residue local divergence with ESM-3 (1.4B, EvolutionaryScale) for one ortholog pair.

Robustness check for the paper (§2.1): does the localized-divergence signal found with
ESM-2 (8M) persist with a much larger, different-architecture model (ESM-3 1.4B)?

Method: embed each full sequence once with ESM-3 (per-residue 1536-d representations),
then for each zebrafish residue i take the cosine distance to its aligned aesculapii
residue as the local divergence. The residue<->residue mapping is taken from the
matching ESM-2 divergence CSV (its `aescallii_residue` column, from the NW alignment),
so ESM-2 and ESM-3 profiles are directly comparable per residue.

Runs in the env-esm3 environment (esm>=3, separate from fair-esm). CPU works; a 1000-aa
forward on the 1.4B model is a few minutes on Apple Silicon CPU.

Usage:
    ~/env-esm3/bin/python analysis/esm3_residue_divergence.py \\
        --zebrafish-id NP_001039014 --aescallii-id XP_056330517 \\
        --fasta data/reference/sequences/kcnj13_pair.fasta \\
        --mapping-csv reports/kcnj13_residue_divergence.csv \\
        --out reports/kcnj13_residue_divergence_esm3.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


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


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    denom = (np.linalg.norm(a) * np.linalg.norm(b)) or 1e-12
    return float(1.0 - float(np.dot(a, b)) / denom)


def embed_per_residue(model, seq: str) -> np.ndarray:
    """ESM-3 per-residue embeddings [L, d], excluding BOS/EOS."""
    import torch
    from esm.sdk.api import ESMProtein, LogitsConfig
    with torch.no_grad():
        t = model.encode(ESMProtein(sequence=seq))
        out = model.logits(t, LogitsConfig(sequence=True, return_embeddings=True))
    return out.embeddings[0, 1:len(seq) + 1].float().cpu().numpy()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zebrafish-id", required=True)
    ap.add_argument("--aescallii-id", required=True)
    ap.add_argument("--fasta", required=True, help="FASTA containing both sequences.")
    ap.add_argument("--mapping-csv", required=True,
                    help="ESM-2 divergence CSV with columns residue, aescallii_residue (NW mapping).")
    ap.add_argument("--device", default="cpu", help="cpu | mps (cpu recommended for stability).")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(Path(args.fasta))
    zseq, aseq = seqs[args.zebrafish_id], seqs[args.aescallii_id]
    mapping = pd.read_csv(args.mapping_csv)[["residue", "aescallii_residue"]]

    from esm.pretrained import ESM3_sm_open_v0
    model = ESM3_sm_open_v0(args.device)
    print(f"loaded ESM3 on {args.device}; embedding {args.zebrafish_id} ({len(zseq)} aa) "
          f"and {args.aescallii_id} ({len(aseq)} aa)...", flush=True)
    z_emb = embed_per_residue(model, zseq)
    a_emb = embed_per_residue(model, aseq)

    rows = []
    for _, r in mapping.iterrows():
        i, j = int(r["residue"]), int(r["aescallii_residue"])
        if j == 0 or i < 1 or i > len(zseq) or j < 1 or j > len(aseq):
            continue
        rows.append({"residue": i,
                     "divergence_esm3": cosine_distance(z_emb[i - 1], a_emb[j - 1]),
                     "aescallii_residue": j})
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    d = df["divergence_esm3"]
    print(f"Wrote {len(df)} residue ESM-3 divergences -> {out}")
    print(f"  ESM-3 divergence: mean={d.mean():.5f} median={d.median():.5f} max={d.max():.5f} "
          f"(peak at residue {int(df.loc[d.idxmax(),'residue'])})")


if __name__ == "__main__":
    main()
