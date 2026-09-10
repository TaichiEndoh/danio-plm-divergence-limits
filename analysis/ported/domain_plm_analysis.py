#!/usr/bin/env python3
"""
Domain-level PLM analysis using approximate coordinates.

This script slices domains from full-length sequences, computes
domain-level identity (global alignment), and domain-level PLM distances.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from numba import njit

import esm


@dataclass
class DomainSpec:
    name: str
    start: int
    end: int


PASB = DomainSpec("PAS_B", 230, 350)
KINASE = DomainSpec("Pkinase_Ty", 600, 900)


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
def _global_identity_numba(a, b):
    n = a.shape[0]
    m = b.shape[0]
    if n == 0 or m == 0:
        return 0.0

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

    i, j = n, m
    matches = 0
    length = 0
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = 1 if a[i - 1] == b[j - 1] else 0
            if dp[i, j] == dp[i - 1, j - 1] + score:
                length += 1
                if score == 1:
                    matches += 1
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i, j] == dp[i - 1, j] - 1:
            length += 1
            i -= 1
            continue
        if j > 0 and dp[i, j] == dp[i, j - 1] - 1:
            length += 1
            j -= 1
            continue
        if i > 0 and j > 0:
            length += 1
            if a[i - 1] == b[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif i > 0:
            length += 1
            i -= 1
        else:
            length += 1
            j -= 1

    return matches / length if length else 0.0


def global_identity(seq1: str, seq2: str) -> float:
    if not seq1 or not seq2:
        return 0.0
    if seq1 == seq2:
        return 1.0
    a = np.frombuffer(seq1.encode("utf-8"), dtype=np.uint8)
    b = np.frombuffer(seq2.encode("utf-8"), dtype=np.uint8)
    return _global_identity_numba(a, b)


def slice_domain(seq: str, spec: DomainSpec) -> Tuple[Optional[str], str, int, int]:
    if not seq:
        return None, "missing sequence", spec.start, spec.end
    if spec.start > len(seq):
        return None, "start beyond sequence length", spec.start, spec.end
    end = min(spec.end, len(seq))
    note = "" if end == spec.end else "truncated to sequence end"
    return seq[spec.start - 1 : end], note, spec.start, end


def cosine_distance(vec_a: torch.Tensor, vec_b: torch.Tensor) -> Optional[float]:
    if vec_a is None or vec_b is None:
        return None
    a = vec_a.unsqueeze(0)
    b = vec_b.unsqueeze(0)
    sim = torch.nn.functional.cosine_similarity(a, b)
    return float(1 - sim.item())


def embed_sequences(
    sequences: Dict[str, str],
    model,
    alphabet,
    device: str,
    batch_size: int,
) -> Dict[str, torch.Tensor]:
    batch_converter = alphabet.get_batch_converter()
    items = list(sequences.items())
    embeddings: Dict[str, torch.Tensor] = {}
    for start in range(0, len(items), batch_size):
        chunk = items[start : start + batch_size]
        batch = [(pid, seq) for pid, seq in chunk]
        labels, strs, tokens = batch_converter(batch)
        tokens = tokens.to(device)
        with torch.no_grad():
            results = model(tokens, repr_layers=[model.num_layers])
        reps = results["representations"][model.num_layers].cpu()
        for idx, (pid, seq) in enumerate(chunk):
            seq_len = len(seq)
            token_repr = reps[idx, 1 : seq_len + 1].mean(0)
            embeddings[pid] = token_repr
    return embeddings


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict-csv", required=True)
    parser.add_argument("--zebrafish-fasta", required=True)
    parser.add_argument("--aescallii-fasta", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--output-plot", required=True)
    parser.add_argument("--model-name", default="esm2_t6_8M_UR50D")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=1)
    args = parser.parse_args()

    df = pd.read_csv(args.strict_csv)
    df = df[df["note"].isna()]

    pasb_symbols = {"AHR1B", "AHR2", "AHR1A", "AHRRA", "AHRRB", "NR1I2"}
    kinase_symbols = {"LTK", "CSF1RA", "KITA", "EDNRB1A", "EDNRB1B"}

    zseqs = load_fasta(Path(args.zebrafish_fasta))
    aseqs = load_fasta(Path(args.aescallii_fasta))

    rows = []
    domain_sequences = {}

    for row in df.itertuples(index=False):
        symbol = row.symbol
        if symbol in pasb_symbols:
            spec = PASB
        elif symbol in kinase_symbols:
            spec = KINASE
        else:
            # skip symbols without requested domain definition (e.g., CYP1A)
            continue

        z_id = row.zebrafish_id
        a_id = row.aescallii_id

        zseq = zseqs.get(z_id)
        aseq = aseqs.get(a_id)

        z_dom, z_note, z_start, z_end = slice_domain(zseq, spec)
        a_dom, a_note, a_start, a_end = slice_domain(aseq, spec)

        note = "; ".join([n for n in [z_note, a_note] if n])
        if not z_dom or not a_dom:
            rows.append({
                "group": row.group,
                "symbol": symbol,
                "domain": spec.name,
                "zebrafish_id": z_id,
                "aescallii_id": a_id,
                "z_start": z_start,
                "z_end": z_end,
                "a_start": a_start,
                "a_end": a_end,
                "identity_percent_domain": None,
                "plm_distance_domain": None,
                "note": note or "missing domain sequence",
            })
            continue

        identity = global_identity(z_dom, a_dom) * 100.0

        z_key = f"{z_id}|{spec.name}|{z_start}-{z_end}"
        a_key = f"{a_id}|{spec.name}|{a_start}-{a_end}"
        domain_sequences[z_key] = z_dom
        domain_sequences[a_key] = a_dom

        rows.append({
            "group": row.group,
            "symbol": symbol,
            "domain": spec.name,
            "zebrafish_id": z_id,
            "aescallii_id": a_id,
            "z_start": z_start,
            "z_end": z_end,
            "a_start": a_start,
            "a_end": a_end,
            "identity_percent_domain": identity,
            "plm_distance_domain": None,
            "note": note,
        })

    if not domain_sequences:
        raise SystemExit("No domain sequences to embed.")

    model, alphabet = esm.pretrained.load_model_and_alphabet(args.model_name)
    model = model.to(args.device)
    model.eval()

    embeddings = embed_sequences(domain_sequences, model, alphabet, args.device, args.batch_size)

    for rec in rows:
        if rec["identity_percent_domain"] is None:
            continue
        z_key = f"{rec['zebrafish_id']}|{rec['domain']}|{rec['z_start']}-{rec['z_end']}"
        a_key = f"{rec['aescallii_id']}|{rec['domain']}|{rec['a_start']}-{rec['a_end']}"
        vec_a = embeddings.get(z_key)
        vec_b = embeddings.get(a_key)
        rec["plm_distance_domain"] = cosine_distance(vec_a, vec_b)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output_csv, index=False)

    # scatter plot
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise SystemExit("matplotlib is required for plotting.")

    plot_df = out_df.dropna(subset=["identity_percent_domain", "plm_distance_domain"])
    colors = {"conserved": "#1f77b4", "variant": "#d62728"}
    fig, ax = plt.subplots(figsize=(6, 4))
    for group, gdf in plot_df.groupby("group"):
        ax.scatter(
            gdf["identity_percent_domain"],
            gdf["plm_distance_domain"],
            label=group,
            c=colors.get(group, "#333333"),
            alpha=0.8,
        )
    ax.set_xlabel("Domain identity (%)")
    ax.set_ylabel("Domain PLM distance (cosine)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.output_plot, dpi=200)


if __name__ == "__main__":
    main()
