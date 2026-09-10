#!/usr/bin/env python3
"""
ESM-2 variant-effect scoring (masked-marginal Δlog-likelihood) — the core primitive for
the "rediscover known functional variants" stage of the project.

For a substitution wt->mut at position p, the score is
    score = log P(mut | context) - log P(wt | context)
read from the ESM-2 masked-language-model head (Meier et al. 2021, ESM-1v). More negative =
the model finds the mutant less plausible = candidate functional impact. Two schemes:

  --scoring wt-marginal      one forward pass over the wild-type sequence; read all positions'
                             logits. Fast; good for a whole deep-mutational scan.
  --scoring masked-marginal  mask each scored position and run a forward per position. Slower
                             (one pass per position) but the stronger scheme in the ESM papers.

Two jobs:
  1. --variants A111T,S102A   score specific substitutions (1-based positions).
  2. --scan                   score every single substitution (L x 19); with --known VARIANT it
                              also reports the known variant's percentile rank among all
                              substitutions in that protein — the Stage-2 "rediscovery" test.

Positive controls to validate the method (see docs/PAPER_DESIGN_ja.md): SLC24A5 A111T,
rhodopsin E122Q/A292S, ChATa S102, HBB E6V, ...

Usage:
    python analysis/esm2_variant_score.py --fasta seqs.fasta --id SLC24A5 \\
        --scan --known A111T --model-name esm2_t33_650M_UR50D \\
        --scoring masked-marginal --out reports/slc24a5_dms.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


def parse_variant(v: str) -> Tuple[str, int, str]:
    """'A111T' -> ('A', 111, 'T') (1-based position)."""
    v = v.strip().upper()
    wt, mut, pos = v[0], v[-1], int(v[1:-1])
    if wt not in AA or mut not in AA:
        raise SystemExit(f"bad variant {v}")
    return wt, pos, mut


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence")
    src.add_argument("--fasta")
    ap.add_argument("--id", help="sequence id within --fasta.")
    ap.add_argument("--variants", help="comma-separated, e.g. A111T,S102A.")
    ap.add_argument("--scan", action="store_true", help="score every single substitution (L x 19).")
    ap.add_argument("--known", help="with --scan: report this variant's percentile rank.")
    ap.add_argument("--model-name", default="esm2_t33_650M_UR50D")
    ap.add_argument("--scoring", choices=["wt-marginal", "masked-marginal"], default="masked-marginal")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if not args.variants and not args.scan:
        ap.error("provide --variants and/or --scan")

    if args.sequence:
        seq = args.sequence.strip().upper()
    else:
        if not args.id:
            ap.error("--fasta requires --id")
        seq = read_fasta(Path(args.fasta))[args.id]
    L = len(seq)

    import torch
    import esm
    model, alphabet = esm.pretrained.load_model_and_alphabet(args.model_name)
    model = model.to(args.device).eval()
    bc = alphabet.get_batch_converter()
    _, _, base_tokens = bc([("wt", seq)])
    base_tokens = base_tokens.to(args.device)
    aa_idx = {a: alphabet.get_idx(a) for a in AA}
    mask_idx = alphabet.mask_idx

    # which positions do we need logits for?
    if args.scan:
        positions = list(range(1, L + 1))
    else:
        positions = sorted({parse_variant(v)[1] for v in args.variants.split(",")})

    # log-softmax logits per needed position (token index = pos, because of the leading BOS)
    logprobs: Dict[int, np.ndarray] = {}
    with torch.no_grad():
        if args.scoring == "wt-marginal":
            lp = torch.log_softmax(model(base_tokens)["logits"], dim=-1)[0].cpu().numpy()
            for p in positions:
                logprobs[p] = lp[p]
        else:  # masked-marginal: one forward per masked position
            for p in positions:
                t = base_tokens.clone()
                t[0, p] = mask_idx
                lp = torch.log_softmax(model(t)["logits"], dim=-1)[0, p].cpu().numpy()
                logprobs[p] = lp

    def score(wt: str, pos: int, mut: str) -> float:
        lp = logprobs[pos]
        return float(lp[aa_idx[mut]] - lp[aa_idx[wt]])

    rows: List[dict] = []
    if args.scan:
        for pos in positions:
            wt = seq[pos - 1]
            for mut in AA:
                if mut == wt:
                    continue
                rows.append({"variant": f"{wt}{pos}{mut}", "position": pos,
                             "wt": wt, "mut": mut, "score": score(wt, pos, mut)})
    if args.variants:
        for v in args.variants.split(","):
            wt, pos, mut = parse_variant(v)
            if seq[pos - 1] != wt:
                print(f"  WARNING: {v} wt {wt} != sequence residue {seq[pos-1]} at {pos}")
            rows.append({"variant": f"{wt}{pos}{mut}", "position": pos,
                         "wt": wt, "mut": mut, "score": score(wt, pos, mut)})

    df = pd.DataFrame(rows).drop_duplicates("variant").sort_values("score").reset_index(drop=True)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote {len(df)} variant scores ({args.id or 'seq'}, {args.model_name}, {args.scoring}) -> {out}")

    if args.scan and args.known:
        wt, pos, mut = parse_variant(args.known)
        kv = f"{wt}{pos}{mut}"
        sub = df[df["variant"] == kv]
        if sub.empty:
            print(f"  known {kv} not in scan (check position/wt).")
        else:
            s = float(sub["score"].iloc[0])
            # rank among all substitutions: lower (more negative) score = more impactful
            pct = 100.0 * (df["score"] <= s).sum() / len(df)
            print(f"  known variant {kv}: score={s:.3f}, ranks in the most-impactful {pct:.1f}% "
                  f"of {len(df)} substitutions (lower % = better rediscovery).")


if __name__ == "__main__":
    main()
