#!/usr/bin/env python3
"""
Per-residue local embedding divergence for one ortholog pair — the bridge from PLM
distance to structure. For every zebrafish residue that aligns to an aesculapii residue,
take a local window (+/- flank), embed both windows, and record the cosine distance as
that residue's local divergence. The output (residue, divergence, ...) feeds directly
into analysis/alphafold_local_mapping.py, which writes it into the B-factor column of the
AlphaFold model so the structure can be colored by divergence.

This generalizes the ported pocket-only analysis (analysis/ported/ahr_pocket_local_plm.py)
to a sliding window over the whole protein, using the same Needleman-Wunsch position map
and the same ESM pooling. Runs on CPU with ESM-2 (the AHR2 pair is ~1k residues -> ~2k
short-window embeddings, seconds on 8M).

Usage:
    python analysis/local_residue_divergence.py \\
        --zebrafish-id NP_571339 --aescallii-id XP_056303610 \\
        --fasta-a data/reference/sequences/target_domain_sequences.fasta \\
        --fasta-b data/reference/sequences/target_domain_sequences.fasta \\
        --model-name esm2_t6_8M_UR50D --flank 10 \\
        --out reports/ahr2_residue_divergence.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur, buf = None, []
    with Path(path).open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0].split(".")[0], []
            elif cur:
                buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def nw_map(seq_a: str, seq_b: str) -> Dict[int, int]:
    """Needleman-Wunsch; return {1-based pos in A -> 1-based pos in B} (0 = gap).

    Unit match/mismatch scoring with a -1 gap, matching the ported pocket analysis.
    """
    n, m = len(seq_a), len(seq_b)
    dp = np.zeros((n + 1, m + 1), dtype=np.int32)
    dp[0, :] = -np.arange(m + 1)
    dp[:, 0] = -np.arange(n + 1)
    for i in range(1, n + 1):
        ai = seq_a[i - 1]
        row, prev = dp[i], dp[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (1 if ai == seq_b[j - 1] else 0)
            row[j] = max(diag, prev[j] - 1, row[j - 1] - 1)
    mapping: Dict[int, int] = {}
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i, j] == dp[i - 1, j - 1] + (1 if seq_a[i - 1] == seq_b[j - 1] else 0):
            mapping[i] = j
            i -= 1
            j -= 1
        elif i > 0 and dp[i, j] == dp[i - 1, j] - 1:
            mapping[i] = 0
            i -= 1
        else:
            j -= 1
    return mapping


def window(seq: str, center: int, flank: int) -> Optional[str]:
    """Full-width (2*flank+1) window centred on `center` (1-based), or None near a terminus.

    Returning a clipped window here would be a bug: the two species' termini are not
    generally aligned (an N-terminal indel shifts one sequence relative to the other), so a
    clipped window on one side would be compared against a differently-sized, differently-framed
    window on the other, and the resulting cosine distance would reflect window framing rather
    than sequence divergence. Residues within `flank` of either terminus are therefore skipped.
    """
    if center <= 0 or center > len(seq):
        return None
    start = center - 1 - flank
    end = center + flank
    if start < 0 or end > len(seq):
        return None
    return seq[start:end]


def embed_windows(win: Dict[str, str], model_name: str, device: str, batch_size: int) -> Dict[str, np.ndarray]:
    import torch
    import esm
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device).eval()
    bc = alphabet.get_batch_converter()
    last = model.num_layers
    items = list(win.items())
    out: Dict[str, np.ndarray] = {}
    for s in range(0, len(items), batch_size):
        chunk = items[s:s + batch_size]
        _, _, toks = bc(chunk)
        toks = toks.to(device)
        with torch.no_grad():
            rep = model(toks, repr_layers=[last])["representations"][last].cpu()
        for i, (key, seq) in enumerate(chunk):
            out[key] = rep[i, 1:len(seq) + 1].mean(0).numpy().astype(np.float32)
    return out


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    denom = (np.linalg.norm(a) * np.linalg.norm(b)) or 1e-12
    return float(1.0 - float(np.dot(a, b)) / denom)


def local_identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    return 100.0 * sum(1 for i in range(n) if a[i] == b[i]) / n if n else 0.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zebrafish-id", required=True)
    ap.add_argument("--aescallii-id", required=True)
    ap.add_argument("--fasta-a", required=True, help="FASTA containing the zebrafish sequence.")
    ap.add_argument("--fasta-b", required=True, help="FASTA containing the aesculapii sequence.")
    ap.add_argument("--out", required=True, help="output CSV: residue,divergence,local_identity,aescallii_residue.")
    ap.add_argument("--model-name", default="esm2_t6_8M_UR50D")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--flank", type=int, default=10, help="half-window size in residues.")
    ap.add_argument("--batch-size", type=int, default=8)
    args = ap.parse_args()

    a_seqs = read_fasta(Path(args.fasta_a))
    b_seqs = read_fasta(Path(args.fasta_b))
    zseq = a_seqs.get(args.zebrafish_id)
    aseq = b_seqs.get(args.aescallii_id)
    if not zseq:
        raise SystemExit(f"{args.zebrafish_id} not found in {args.fasta_a}")
    if not aseq:
        raise SystemExit(f"{args.aescallii_id} not found in {args.fasta_b}")

    mapping = nw_map(zseq, aseq)
    win: Dict[str, str] = {}
    plan: List[tuple] = []  # (zebrafish_residue, aescallii_residue, zkey, akey)
    for zpos in range(1, len(zseq) + 1):
        apos = mapping.get(zpos, 0)
        if apos == 0:
            continue
        zw, aw = window(zseq, zpos, args.flank), window(aseq, apos, args.flank)
        if not zw or not aw:
            continue
        zk, ak = f"z{zpos}", f"a{apos}"
        win[zk], win[ak] = zw, aw
        plan.append((zpos, apos, zk, ak))

    emb = embed_windows(win, args.model_name, args.device, args.batch_size)

    rows = []
    for zpos, apos, zk, ak in plan:
        rows.append({
            "residue": zpos,
            "divergence": cosine_distance(emb[zk], emb[ak]),
            "local_identity": local_identity(win[zk], win[ak]),
            "aescallii_residue": apos,
        })
    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote {len(df)} residue divergences ({args.zebrafish_id} vs {args.aescallii_id}) -> {out}")
    if len(df):
        d = df["divergence"]
        print(f"  divergence: mean={d.mean():.5f} median={d.median():.5f} max={d.max():.5f} "
              f"(peak at residue {int(df.loc[d.idxmax(), 'residue'])})")
        print(f"  {int((mapping and len(zseq)) - len(df))} residues unmapped/edge-skipped of {len(zseq)}")


if __name__ == "__main__":
    main()
