#!/usr/bin/env python3
"""
ESM-3 (1.4B) sequence embedding + pairwise cosine distance — starter skeleton.

Goal: re-embed ortholog pairs (and local windows) with ESM-3 to compare against the
previous ESM-2 8M results. This is a SKELETON: the heavy embedding step requires the
`esm` package, accepted model weights, and (practically) a GPU. A deterministic
`--demo` mode runs the full data flow without the model so the shape can be tested.

Usage:
    # real run (needs GPU + esm-3 weights):
    python analysis/esm3_embed.py --pairs data/orthologs/pairs.csv \\
        --fasta-a data/sequences/rerio.fasta --fasta-b data/sequences/aesculapii.fasta \\
        --out reports/esm3_pair_distances.csv

    # demo (no model, deterministic fake embeddings):
    python analysis/esm3_embed.py --pairs data/orthologs/pairs.csv \\
        --fasta-a data/sequences/rerio.fasta --fasta-b data/sequences/aesculapii.fasta \\
        --out /tmp/demo.csv --demo
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

MODEL_ID = "esm3-open"  # 1.4B open weights; record the exact id you actually load


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


# ----- embedding backends -------------------------------------------------

def embed_demo(seq: str, dim: int = 320) -> np.ndarray:
    """Deterministic pseudo-embedding (seeded by sequence). NOT science — demo only."""
    seed = int(hashlib.sha256(seq.encode()).hexdigest()[:8], 16)
    return np.random.default_rng(seed).standard_normal(dim).astype(np.float32)


def embed_esm3(sequences: Dict[str, str]) -> Dict[str, np.ndarray]:
    """
    TODO: real ESM-3 embedding.
      - from esm.models.esm3 import ESM3 (or the current API — CHECK the installed version)
      - load esm3-open (1.4B), move to CUDA, eval()
      - for each sequence: tokenize, forward, take per-residue representations,
        mean-pool over residues excluding special tokens (mirror the previous repo's pooling).
    Returns {id: 1D embedding}.
    """
    raise NotImplementedError(
        "Real ESM-3 embedding not implemented yet. Install `esm`, accept weights, "
        "run on GPU, and mirror the previous repo's mean-pooling. Use --demo to test the flow."
    )


# ----- main ---------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", required=True, help="CSV with columns: zebrafish_id, aescallii_id")
    ap.add_argument("--fasta-a", required=True, help="D. rerio protein FASTA")
    ap.add_argument("--fasta-b", required=True, help="D. aesculapii protein FASTA")
    ap.add_argument("--out", required=True, help="output CSV of per-pair ESM-3 cosine distance")
    ap.add_argument("--demo", action="store_true", help="use deterministic fake embeddings (no model)")
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs)
    a_seqs = read_fasta(Path(args.fasta_a))
    b_seqs = read_fasta(Path(args.fasta_b))

    # collect the sequences we actually need
    needed_a = {r.zebrafish_id: a_seqs[r.zebrafish_id] for r in pairs.itertuples()
                if r.zebrafish_id in a_seqs}
    needed_b = {r.aescallii_id: b_seqs[r.aescallii_id] for r in pairs.itertuples()
                if r.aescallii_id in b_seqs}

    if args.demo:
        emb_a = {k: embed_demo(v) for k, v in needed_a.items()}
        emb_b = {k: embed_demo(v) for k, v in needed_b.items()}
    else:
        emb_a = embed_esm3(needed_a)
        emb_b = embed_esm3(needed_b)

    rows = []
    for r in pairs.itertuples():
        va, vb = emb_a.get(r.zebrafish_id), emb_b.get(r.aescallii_id)
        if va is None or vb is None:
            continue
        rows.append({
            "zebrafish_id": r.zebrafish_id,
            "aescallii_id": r.aescallii_id,
            "esm3_cosine_distance": cosine_distance(va, vb),
            "model": "demo" if args.demo else MODEL_ID,
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote {len(rows)} pair distances to {out}  (model={'demo' if args.demo else MODEL_ID})")


if __name__ == "__main__":
    main()
