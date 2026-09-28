#!/usr/bin/env python3
"""
ESM-2 ladder: per-pair cosine distance at 8M / 35M / 150M / 650M — the paper's MAIN
scaling axis. Running the same ortholog pairs across four ESM-2 sizes isolates the
"is 8M enough (scale)?" question within a single model family; ESM-3 (a different family)
is then the architecture-robustness confirmation, not the primary evidence.

Pooling is IDENTICAL to the previous repo (final layer, exclude the BOS/EOS special
tokens, mean over residues) so the 8M run reproduces the ported baseline. Run one size
per invocation (clean GPU-memory profile, independently resumable); loop over the four
names for the full ladder:

    esm2_t6_8M_UR50D  esm2_t12_35M_UR50D  esm2_t30_150M_UR50D  esm2_t33_650M_UR50D

Usage (real run, GPU recommended for 150M/650M):
    python analysis/esm2_ladder_embed.py \\
        --pairs data/reference/esm2_8m_pair_distances_dedup.csv \\
        --fasta-a data/sequences/rerio.fasta --fasta-b data/sequences/aesculapii.fasta \\
        --model-name esm2_t33_650M_UR50D --cache-dir data/cache/esm2_650m \\
        --out reports/esm2_650m_pair_distances.csv

    # flow test without weights (deterministic fake embeddings):
    python analysis/esm2_ladder_embed.py --pairs P.csv \\
        --fasta-a A.fasta --fasta-b B.fasta --out /tmp/demo.csv --demo
"""
from __future__ import annotations

import argparse
import hashlib
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# short tag -> fair-esm checkpoint name (for logging / output labelling)
LADDER = {
    "8M": "esm2_t6_8M_UR50D",
    "35M": "esm2_t12_35M_UR50D",
    "150M": "esm2_t30_150M_UR50D",
    "650M": "esm2_t33_650M_UR50D",
}


def read_fasta(path: Path) -> Dict[str, str]:
    """{id_without_version: sequence}. Mirrors the ported normalize_id (drop '.N')."""
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


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    denom = (np.linalg.norm(a) * np.linalg.norm(b)) or 1e-12
    return float(1.0 - float(np.dot(a, b)) / denom)


def embed_demo(seq: str, dim: int = 320) -> np.ndarray:
    seed = int(hashlib.sha256(seq.encode()).hexdigest()[:8], 16)
    return np.random.default_rng(seed).standard_normal(dim).astype(np.float32)


def embed_esm2(ids: List[str], seqs: Dict[str, str], model_name: str, device: str,
               batch_size: int, cache_dir: Optional[Path], max_seq_len: Optional[int]) -> Dict[str, np.ndarray]:
    """Mean-pool final-layer per-residue reps, excluding special tokens. == previous repo."""
    import torch
    import esm

    log = logging.getLogger("esm2_ladder")
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)

    out: Dict[str, np.ndarray] = {}
    to_encode: List[Tuple[str, str]] = []
    for pid in ids:
        if cache_dir and (cache_dir / f"{pid}.npy").is_file():
            out[pid] = np.load(cache_dir / f"{pid}.npy")
            continue
        seq = seqs.get(pid)
        if not seq:
            log.warning("sequence for %s missing; pair distance will be NaN", pid)
            continue
        if max_seq_len and len(seq) > max_seq_len:
            seq = seq[:max_seq_len]
        to_encode.append((pid, seq))

    if to_encode:
        if device == "cuda" and not torch.cuda.is_available():
            raise SystemExit("CUDA requested but not available; re-run with --device cpu.")
        log.info("loading %s on %s; encoding %d sequences", model_name, device, len(to_encode))
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model = model.to(device).eval()
        bc = alphabet.get_batch_converter()
        last = model.num_layers
        for s in range(0, len(to_encode), batch_size):
            chunk = to_encode[s:s + batch_size]
            _, _, toks = bc(chunk)
            toks = toks.to(device)
            with torch.no_grad():
                rep = model(toks, repr_layers=[last])["representations"][last].cpu()
            for i, (pid, seq) in enumerate(chunk):
                emb = rep[i, 1:len(seq) + 1].mean(0).numpy().astype(np.float32)
                out[pid] = emb
                if cache_dir:
                    np.save(cache_dir / f"{pid}.npy", emb)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True, help="CSV with zebrafish_id, aescallii_id columns.")
    ap.add_argument("--fasta-a", required=True, help="D. rerio protein FASTA (full-length).")
    ap.add_argument("--fasta-b", required=True, help="D. aesculapii protein FASTA (full-length).")
    ap.add_argument("--out", required=True, help="output per-pair distance CSV.")
    ap.add_argument("--model-name", default="esm2_t6_8M_UR50D", help="fair-esm checkpoint (see LADDER).")
    ap.add_argument("--device", default="auto", help="cuda | cpu | mps | auto.")
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--cache-dir", default=None, help="optional per-protein embedding cache (*.npy).")
    ap.add_argument("--max-seq-len", type=int, default=None)
    ap.add_argument("--num-threads", type=int, default=1)
    ap.add_argument("--row-offset", type=int, default=0, help="start pair row (resumable runs).")
    ap.add_argument("--row-limit", type=int, default=None, help="max pair rows this run.")
    ap.add_argument("--demo", action="store_true", help="deterministic fake embeddings (no model).")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    device = args.device
    if device == "auto":
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else (
                "mps" if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available() else "cpu")
        except Exception:
            device = "cpu"

    if not args.demo:
        try:
            import torch
            torch.set_num_threads(max(1, args.num_threads))
        except Exception:
            pass

    pairs = pd.read_csv(args.pairs)
    for col in ("zebrafish_id", "aescallii_id"):
        if col not in pairs.columns:
            raise SystemExit(f"pairs CSV must contain '{col}'")
    pairs["zebrafish_id"] = pairs["zebrafish_id"].astype(str).str.split(".").str[0]
    pairs["aescallii_id"] = pairs["aescallii_id"].astype(str).str.split(".").str[0]
    if args.row_offset or args.row_limit is not None:
        end = args.row_offset + args.row_limit if args.row_limit is not None else None
        pairs = pairs.iloc[args.row_offset:end].reset_index(drop=True)

    a_seqs = read_fasta(Path(args.fasta_a))
    b_seqs = read_fasta(Path(args.fasta_b))
    need_a = sorted(set(pairs["zebrafish_id"]) & set(a_seqs))
    need_b = sorted(set(pairs["aescallii_id"]) & set(b_seqs))

    # Combine both species into one embedding pass so the model loads once (matters for
    # 650M over 50k pairs). NP_/XP_ id prefixes never collide.
    combined = {**{k: a_seqs[k] for k in need_a}, **{k: b_seqs[k] for k in need_b}}
    if args.demo:
        emb = {k: embed_demo(v) for k, v in combined.items()}
        model_label = "demo"
    else:
        cache = Path(args.cache_dir) if args.cache_dir else None
        emb = embed_esm2(sorted(combined), combined, args.model_name, device,
                         args.batch_size, cache, args.max_seq_len)
        model_label = args.model_name

    rows = []
    for r in pairs.itertuples():
        va, vb = emb.get(r.zebrafish_id), emb.get(r.aescallii_id)
        rows.append({
            "zebrafish_id": r.zebrafish_id,
            "aescallii_id": r.aescallii_id,
            "esm2_cosine_distance": cosine_distance(va, vb) if va is not None and vb is not None else np.nan,
            "model": model_label,
        })
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    n_ok = int(df["esm2_cosine_distance"].notna().sum())
    print(f"Wrote {len(df)} pairs ({n_ok} with distance) to {out}  (model={model_label}, device={device})")


if __name__ == "__main__":
    main()
