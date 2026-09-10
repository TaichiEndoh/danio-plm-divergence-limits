#!/usr/bin/env python3
"""
Local PLM distance around pocket residues for AHR subtypes.

For each zebrafish AHR protein, map pocket residues to aescallii candidates
via global alignment, extract local windows (±10 aa), embed with ESM-2,
and compute mean/median PLM distances per pair.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from numba import njit

import esm


POCKET_RESIDUES = {
    "NP_571339": [
        281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356,
        370, 378, 380, 381,
    ],
    "NP_571337": [
        278, 280, 286, 288, 300, 304, 308, 311, 319, 321, 332, 337, 345, 353,
        367, 375, 377, 378,
    ],
    # AHR1a full-length RefSeq (NP_001417877.1); same residue indices as NP_571337
    "NP_001417877": [
        278, 280, 286, 288, 300, 304, 308, 311, 319, 321, 332, 337, 345, 353,
        367, 375, 377, 378,
    ],
    "NP_571338": [
        281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356,
        370, 378, 380, 381,
    ],
}


def load_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur_id = None
    cur: List[str] = []
    with path.open() as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur)
                cur_id = line[1:].strip().split()[0].split(".")[0]
                cur = []
            else:
                if cur_id:
                    cur.append(line.strip())
        if cur_id:
            seqs[cur_id] = "".join(cur)
    return seqs


@njit
def _nw_traceback(a, b):
    n = a.shape[0]
    m = b.shape[0]
    dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    for j in range(m + 1):
        dp[0, j] = -j
    for i in range(n + 1):
        dp[i, 0] = -i
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = 1 if a[i - 1] == b[j - 1] else 0
            diag = dp[i - 1, j - 1] + score
            up = dp[i - 1, j] - 1
            left = dp[i, j - 1] - 1
            if diag >= up and diag >= left:
                dp[i, j] = diag
            elif up >= left:
                dp[i, j] = up
            else:
                dp[i, j] = left

    # map positions (1-based) in seq A to seq B
    i, j = n, m
    a_pos = n
    b_pos = m
    mapping = np.zeros(n + 1, dtype=np.int32)  # 0 = gap
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = 1 if a[i - 1] == b[j - 1] else 0
            if dp[i, j] == dp[i - 1, j - 1] + score:
                mapping[a_pos] = b_pos
                i -= 1
                j -= 1
                a_pos -= 1
                b_pos -= 1
                continue
        if i > 0 and dp[i, j] == dp[i - 1, j] - 1:
            mapping[a_pos] = 0
            i -= 1
            a_pos -= 1
            continue
        if j > 0 and dp[i, j] == dp[i, j - 1] - 1:
            j -= 1
            b_pos -= 1
            continue
        if i > 0 and j > 0:
            mapping[a_pos] = b_pos
            i -= 1
            j -= 1
            a_pos -= 1
            b_pos -= 1
        elif i > 0:
            mapping[a_pos] = 0
            i -= 1
            a_pos -= 1
        else:
            j -= 1
            b_pos -= 1
    return mapping


def map_positions(seq_a: str, seq_b: str) -> Dict[int, int]:
    a = np.frombuffer(seq_a.encode("utf-8"), dtype=np.uint8)
    b = np.frombuffer(seq_b.encode("utf-8"), dtype=np.uint8)
    mapping = _nw_traceback(a, b)
    return {i: int(mapping[i]) for i in range(1, len(mapping))}


def extract_window(seq: str, center: int, flank: int) -> Optional[str]:
    if center <= 0 or center > len(seq):
        return None
    start = max(1, center - flank)
    end = min(len(seq), center + flank)
    return seq[start - 1 : end]


def embed_sequences(sequences: Dict[str, str], model, alphabet, device: str) -> Dict[str, torch.Tensor]:
    batch_converter = alphabet.get_batch_converter()
    items = list(sequences.items())
    embeddings: Dict[str, torch.Tensor] = {}
    for start in range(0, len(items), 1):
        pid, seq = items[start]
        labels, strs, tokens = batch_converter([(pid, seq)])
        tokens = tokens.to(device)
        with torch.no_grad():
            results = model(tokens, repr_layers=[model.num_layers])
        reps = results["representations"][model.num_layers].cpu()
        token_repr = reps[0, 1 : len(seq) + 1].mean(0)
        embeddings[pid] = token_repr
    return embeddings


def cosine_distance(vec_a: torch.Tensor, vec_b: torch.Tensor) -> float:
    a = vec_a.unsqueeze(0)
    b = vec_b.unsqueeze(0)
    sim = torch.nn.functional.cosine_similarity(a, b)
    return float(1 - sim.item())


def local_identity(seq_a: str, seq_b: str) -> float:
    if not seq_a or not seq_b:
        return 0.0
    n = min(len(seq_a), len(seq_b))
    if n == 0:
        return 0.0
    matches = sum(1 for i in range(n) if seq_a[i] == seq_b[i])
    return matches / n


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs-csv", required=True)
    parser.add_argument("--zebrafish-fasta", required=True)
    parser.add_argument("--aescallii-fasta", required=True)
    parser.add_argument("--override-fasta", default="")
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--model-name", default="esm2_t6_8M_UR50D")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--flank", type=int, default=10)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs_csv)
    zseqs = load_fasta(Path(args.zebrafish_fasta))
    if args.override_fasta:
        override = load_fasta(Path(args.override_fasta))
        zseqs.update(override)
    aseqs = load_fasta(Path(args.aescallii_fasta))

    model, alphabet = esm.pretrained.load_model_and_alphabet(args.model_name)
    model = model.to(args.device)
    model.eval()

    rows = []
    for row in pairs.itertuples(index=False):
        z_id = row.zebrafish_id
        a_id = row.aescallii_id
        if z_id not in POCKET_RESIDUES:
            continue
        zseq = zseqs.get(z_id)
        aseq = aseqs.get(a_id)
        if not zseq or not aseq:
            continue

        mapping = map_positions(zseq, aseq)
        window_seqs = {}
        window_pairs = []
        for pos in POCKET_RESIDUES[z_id]:
            a_pos = mapping.get(pos, 0)
            if a_pos == 0:
                continue
            z_win = extract_window(zseq, pos, args.flank)
            a_win = extract_window(aseq, a_pos, args.flank)
            if not z_win or not a_win:
                continue
            z_key = f"{z_id}:{pos}"
            a_key = f"{a_id}:{a_pos}"
            window_seqs[z_key] = z_win
            window_seqs[a_key] = a_win
            window_pairs.append((z_key, a_key))

        if not window_pairs:
            continue

        embeddings = embed_sequences(window_seqs, model, alphabet, args.device)
        distances = []
        identities = []
        for z_key, a_key in window_pairs:
            dist = cosine_distance(embeddings[z_key], embeddings[a_key])
            distances.append(dist)
            identities.append(local_identity(window_seqs[z_key], window_seqs[a_key]) * 100.0)

        rows.append(
            {
                "zebrafish_id": z_id,
                "aescallii_id": a_id,
                "mean_local_plm": float(np.mean(distances)),
                "median_local_plm": float(np.median(distances)),
                "mean_local_identity": float(np.mean(identities)),
                "median_local_identity": float(np.median(identities)),
                "n_windows": len(distances),
            }
        )

    out_df = pd.DataFrame(rows)
    out_df = out_df.drop_duplicates(subset=["zebrafish_id", "aescallii_id"])
    out_df.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
