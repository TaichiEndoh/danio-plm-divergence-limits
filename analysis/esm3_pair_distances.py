#!/usr/bin/env python3
"""
Per-pair ESM-3 (1.4B) cosine distance over a set of ortholog pairs — the ESM-3 axis of
paper Fig 2/3 (ESM-2 vs ESM-3 subset comparison).

Mean-pools ESM-3 per-residue embeddings (1536-d, excluding BOS/EOS) per protein, caches
each protein embedding to disk (resumable across restarts), then computes cosine distance
per pair. CPU only (ESM-3 rejects the mps device); ~5-6 s per ~1000-aa protein.

Usage (env-esm3):
    ~/env-esm3/bin/python analysis/esm3_pair_distances.py \\
        --pairs reports/esm3_subset_pairs.csv \\
        --fasta-a data/reference/sequences/subset/subset_rerio.fasta \\
        --fasta-b data/reference/sequences/subset/subset_aesculapii.fasta \\
        --cache-dir data/cache/esm3_1p4b --max-seq-len 2046 \\
        --out reports/esm3_subset_esm3_1p4b.csv
"""
from __future__ import annotations

import argparse
import time
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--fasta-a", required=True, help="D. rerio FASTA (zebrafish_id).")
    ap.add_argument("--fasta-b", required=True, help="D. aesculapii FASTA (aescallii_id).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache-dir", required=True, help="per-protein embedding cache (*.npy); resumable.")
    ap.add_argument("--max-seq-len", type=int, default=2046)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--num-threads", type=int, default=0, help="torch CPU threads (0 = default).")
    ap.add_argument("--log-every", type=int, default=25)
    args = ap.parse_args()

    import torch
    if args.num_threads:
        torch.set_num_threads(args.num_threads)
    from esm.pretrained import ESM3_sm_open_v0
    from esm.sdk.api import ESMProtein, LogitsConfig

    pairs = pd.read_csv(args.pairs)
    pairs["z"] = pairs["zebrafish_id"].astype(str).str.split(".").str[0]
    pairs["a"] = pairs["aescallii_id"].astype(str).str.split(".").str[0]
    a_seqs, b_seqs = read_fasta(Path(args.fasta_a)), read_fasta(Path(args.fasta_b))

    cache = Path(args.cache_dir)
    cache.mkdir(parents=True, exist_ok=True)

    # unique proteins needed, tagged by which fasta they come from
    need = []  # (pid, seq)
    seen = set()
    for pid in pairs["z"]:
        if pid not in seen and pid in a_seqs:
            seen.add(pid); need.append((pid, a_seqs[pid]))
    for pid in pairs["a"]:
        if pid not in seen and pid in b_seqs:
            seen.add(pid); need.append((pid, b_seqs[pid]))
    todo = [(p, s) for p, s in need if not (cache / f"{p}.npy").is_file()]
    print(f"proteins needed: {len(need)} | already cached: {len(need)-len(todo)} | to embed: {len(todo)}",
          flush=True)

    model = ESM3_sm_open_v0(args.device)
    print(f"loaded ESM3 on {args.device}", flush=True)
    t0 = time.time()
    for i, (pid, seq) in enumerate(todo, 1):
        s = seq[:args.max_seq_len] if args.max_seq_len and len(seq) > args.max_seq_len else seq
        try:
            with torch.no_grad():
                t = model.encode(ESMProtein(sequence=s))
                out = model.logits(t, LogitsConfig(sequence=True, return_embeddings=True))
            emb = out.embeddings[0, 1:len(s) + 1].float().mean(0).cpu().numpy().astype(np.float32)
            np.save(cache / f"{pid}.npy", emb)
        except Exception as e:
            print(f"  ERROR embedding {pid} (len {len(s)}): {type(e).__name__}: {str(e)[:100]}", flush=True)
        if i % args.log_every == 0 or i == len(todo):
            el = time.time() - t0
            rate = el / i
            eta = rate * (len(todo) - i)
            print(f"  {i}/{len(todo)} embedded  ({rate:.1f}s/prot, ETA {eta/60:.0f} min)", flush=True)

    # compute pair distances from cached embeddings
    emb_cache: Dict[str, np.ndarray] = {}
    def get(pid):
        if pid not in emb_cache:
            f = cache / f"{pid}.npy"
            emb_cache[pid] = np.load(f) if f.is_file() else None
        return emb_cache[pid]

    rows = []
    for r in pairs.itertuples():
        va, vb = get(r.z), get(r.a)
        rows.append({"zebrafish_id": r.z, "aescallii_id": r.a,
                     "esm3_cosine_distance": cosine_distance(va, vb) if va is not None and vb is not None else np.nan,
                     "model": "esm3-sm-open-v1"})
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    n_ok = int(df["esm3_cosine_distance"].notna().sum())
    print(f"ESM3_PAIRS_DONE  wrote {len(df)} pairs ({n_ok} with distance) -> {out}", flush=True)


if __name__ == "__main__":
    main()
