#!/usr/bin/env python3
"""
Predict a protein structure with ESMFold — the structure route for this project.

Why ESMFold (not AlphaFold DB): the full-length AHR2 (1027 aa, NP_571339) has NO AlphaFold
DB model, and AlphaFold needs an MSA. ESMFold predicts a 3D structure from a SINGLE sequence
(same ESM family as the ESM-2 ladder and ESM-3 used elsewhere here), so we can fold the full
length ourselves with residue numbering that lines up 1:1 with our per-residue divergence.
The predicted PDB carries per-residue pLDDT confidence in the B-factor column, letting us
ignore low-confidence regions when reading structure/pocket (interaction) sites.

Two backends:
  --backend api    ESM Atlas public API (no GPU; sequences up to ~400 aa) — great for the
                   Kcnj13-size targets and for smoke tests.
  --backend local  fair-esm ESMFold (esm.pretrained.esmfold_v1); needs a GPU for long
                   proteins like full-length AHR2. Install: pip install "fair-esm[esmfold]".

Pipe the resulting PDB into analysis/alphafold_local_mapping.py to colour it by divergence
(that step overwrites the B-factor, so keep this PDB if you also want the pLDDT).

Usage:
    # short target / no GPU (Kcnj13, fragments):
    python analysis/esmfold_predict.py --fasta seqs.fasta --id NP_001039014 \\
        --backend api --out reports/kcnj13_esmfold.pdb

    # full-length AHR2 on a local GPU:
    python analysis/esmfold_predict.py --fasta seqs.fasta --id NP_571339 \\
        --backend local --device cuda --chunk-size 128 --out reports/ahr2_esmfold.pdb
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
API_MAX_LEN = 400  # ESM Atlas public limit (approx)


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


def fold_api(seq: str) -> str:
    import requests
    if len(seq) > API_MAX_LEN:
        raise SystemExit(f"sequence is {len(seq)} aa; the ESM Atlas API caps at ~{API_MAX_LEN}. "
                         f"Use --backend local on a GPU for long proteins.")
    r = requests.post(API_URL, data=seq, timeout=300)
    r.raise_for_status()
    if not r.text.lstrip().startswith(("HEADER", "ATOM", "REMARK")):
        raise SystemExit(f"unexpected API response: {r.text[:200]}")
    return r.text


def fold_local(seq: str, device: str, chunk_size: int | None) -> str:
    import torch
    import esm
    model = esm.pretrained.esmfold_v1().eval()
    if device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA requested but not available; try --device cpu (slow) or a GPU box.")
    model = model.to(device)
    if chunk_size:
        model.set_chunk_size(chunk_size)  # trade speed for memory on long sequences
    with torch.no_grad():
        return model.infer_pdb(seq)


def fold_transformers(seq: str, device: str, chunk_size: int | None) -> str:
    """ESMFold via HuggingFace transformers (facebook/esmfold_v1).

    Unlike fair-esm's esmfold, this path does NOT need openfold's compiled CUDA kernels,
    so it installs cleanly in the Docker image (pip install transformers accelerate) and
    runs full-length proteins on the RTX A5000. First run downloads ~2.8 GB of weights.
    """
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    if device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA requested but not available; try --device cpu (slow) or a GPU box.")
    tok = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).eval()
    model = model.to(device)
    if device == "cuda":
        model.esm = model.esm.half()          # half-precision language tower saves VRAM
        torch.backends.cuda.matmul.allow_tf32 = True
    if chunk_size:
        model.trunk.set_chunk_size(chunk_size)  # lower = less VRAM on long sequences
    ids = tok([seq], return_tensors="pt", add_special_tokens=False)["input_ids"].to(device)
    with torch.no_grad():
        out = model(ids)
    return model.output_to_pdb(out)[0]


def mean_plddt(pdb_text: str) -> float:
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and len(line) >= 66:
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return sum(vals) / len(vals) if vals else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence", help="raw amino-acid sequence.")
    src.add_argument("--fasta", help="FASTA file (use with --id).")
    ap.add_argument("--id", help="sequence id within --fasta.")
    ap.add_argument("--out", required=True, help="output PDB path.")
    ap.add_argument("--backend", choices=["api", "local", "transformers"], default="api",
                    help="api=ESM Atlas (≤400aa, no GPU); transformers=HF ESMFold (GPU, no openfold, "
                         "recommended for full-length); local=fair-esm ESMFold (needs openfold).")
    ap.add_argument("--device", default="cuda", help="local/transformers backend: cuda | cpu.")
    ap.add_argument("--chunk-size", type=int, default=None, help="local/transformers: memory/speed tradeoff.")
    args = ap.parse_args()

    if args.sequence:
        seq = args.sequence.strip().upper()
    else:
        if not args.id:
            ap.error("--fasta requires --id")
        seqs = read_fasta(Path(args.fasta))
        if args.id not in seqs:
            raise SystemExit(f"{args.id} not found in {args.fasta}")
        seq = seqs[args.id]

    if args.backend == "api":
        pdb = fold_api(seq)
    elif args.backend == "transformers":
        pdb = fold_transformers(seq, args.device, args.chunk_size)
    else:
        pdb = fold_local(seq, args.device, args.chunk_size)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(pdb)
    print(f"Wrote ESMFold structure ({len(seq)} aa, backend={args.backend}) -> {out}")
    print(f"  mean pLDDT (B-factor of CA): {mean_plddt(pdb):.1f}")
    print("  next: color by divergence with analysis/alphafold_local_mapping.py "
          "(that overwrites B-factor; keep this file for pLDDT).")


if __name__ == "__main__":
    main()
