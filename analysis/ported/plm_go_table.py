#!/usr/bin/env python3
"""
Generate a PLM (ESM-2 embedding) distance table merged with GO annotations for
Danio rerio ↔ Danio aesculapii ortholog pairs.
"""

from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

try:
    import torch
except ImportError as exc:  # pragma: no cover - torch is required at runtime
    raise SystemExit("PyTorch is required. Install it via `pip install torch`.") from exc

try:
    import esm
except ImportError as exc:  # pragma: no cover - esm is required at runtime
    raise SystemExit(
        "facebookresearch/esm is required. Install it via `pip install fair-esm`."
    ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute PLM distances for ortholog pairs and merge with GO annotations."
    )
    parser.add_argument(
        "--ortholog-table",
        required=True,
        help="Path to RBH table (CSV) with Danio rerio / Danio aesculapii IDs.",
    )
    parser.add_argument(
        "--zebrafish-fasta",
        required=True,
        help="FASTA file containing zebrafish (Danio rerio) protein sequences.",
    )
    parser.add_argument(
        "--aescallii-fasta",
        required=True,
        help="FASTA file containing Danio aesculapii protein sequences.",
    )
    parser.add_argument(
        "--zebrafish-go",
        required=True,
        help="TSV from the GO pipeline (e.g., results/zebrafish_protein_go_annotated.tsv).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path for the combined GO + PLM table.",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to the output file instead of overwriting (no header will be written).",
    )
    parser.add_argument(
        "--model-name",
        default="esm2_t6_8M_UR50D",
        help="ESM-2 checkpoint to load (default: %(default)s).",
    )
    parser.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="Device to run inference on (default: auto-detected).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=4,
        help="Number of sequences per PLM batch.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Optional directory to cache per-protein embeddings (*.pt).",
    )
    parser.add_argument(
        "--max-seq-len",
        type=int,
        default=None,
        help="Optional max amino-acid length; longer sequences are truncated to reduce memory.",
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=1,
        help="Limit PyTorch CPU threads to reduce memory pressure.",
    )
    parser.add_argument(
        "--row-offset",
        type=int,
        default=0,
        help="Start processing from this zero-based row index of the ortholog table.",
    )
    parser.add_argument(
        "--row-limit",
        type=int,
        default=None,
        help="Maximum number of ortholog rows to process in this run.",
    )
    return parser.parse_args()


def normalize_id(raw_id: str) -> str:
    """Drop version suffixes such as '.2' to match ortholog table entries."""
    return raw_id.split(".")[0]


def load_fasta_sequences(fasta_path: Path, target_ids: Optional[Sequence[str]] = None) -> Dict[str, str]:
    """Return {protein_id_without_version: sequence} for a FASTA file."""
    sequences: Dict[str, str] = {}
    target_set = set(target_ids) if target_ids else None
    seq_id: Optional[str] = None
    chunks: List[str] = []
    with fasta_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if seq_id and chunks:
                    base_id = normalize_id(seq_id)
                    if target_set is None or base_id in target_set:
                        sequences[base_id] = "".join(chunks)
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if seq_id and chunks:
            base_id = normalize_id(seq_id)
            if target_set is None or base_id in target_set:
                sequences[base_id] = "".join(chunks)
    return sequences


class EmbeddingCache:
    """Simple filesystem cache for per-protein embeddings."""

    def __init__(self, cache_dir: Optional[Path]):
        self.cache_dir = cache_dir
        if cache_dir:
            cache_dir.mkdir(parents=True, exist_ok=True)

    def load(self, protein_id: str) -> Optional[torch.Tensor]:
        if not self.cache_dir:
            return None
        path = self.cache_dir / f"{protein_id}.pt"
        if not path.is_file():
            return None
        tensor = torch.load(path, map_location="cpu")
        return tensor if isinstance(tensor, torch.Tensor) else torch.tensor(tensor)

    def save(self, protein_id: str, tensor: torch.Tensor) -> None:
        if not self.cache_dir:
            return
        path = self.cache_dir / f"{protein_id}.pt"
        torch.save(tensor.cpu(), path)


def init_logger() -> logging.Logger:
    logger = logging.getLogger("plm_go")
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(handler)
    return logger


def load_model(model_name: str, device: str):
    if device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA requested but not available. Re-run with --device cpu.")
    logger = logging.getLogger("plm_go")
    logger.info("Loading ESM-2 model %s on %s", model_name, device)
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.to(device)
    model.eval()
    return model, alphabet


def generate_embeddings(
    protein_ids: Sequence[str],
    sequences: Dict[str, str],
    model,
    alphabet,
    device: str,
    batch_size: int,
    cache: EmbeddingCache,
    max_seq_len: Optional[int],
) -> Dict[str, torch.Tensor]:
    logger = logging.getLogger("plm_go")
    embeddings: Dict[str, torch.Tensor] = {}
    batch_converter = alphabet.get_batch_converter()
    to_encode: List[Tuple[str, str]] = []
    for pid in protein_ids:
        cached = cache.load(pid)
        if cached is not None:
            embeddings[pid] = cached
            continue
        seq = sequences.get(pid)
        if not seq:
            logger.warning("Sequence for %s not found; distances will be NaN", pid)
            continue
        if max_seq_len and len(seq) > max_seq_len:
            seq = seq[:max_seq_len]
        to_encode.append((pid, seq))

    logger.info("Encoding %d sequences (batch size=%d)", len(to_encode), batch_size)
    for start in range(0, len(to_encode), batch_size):
        chunk = to_encode[start : start + batch_size]
        batch = [(pid, seq) for pid, seq in chunk]
        batch_labels, batch_strs, batch_tokens = batch_converter(batch)
        batch_tokens = batch_tokens.to(device)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[model.num_layers])
        representations = results["representations"][model.num_layers].cpu()
        for idx, (pid, seq) in enumerate(chunk):
            seq_len = len(seq)
            token_repr = representations[idx, 1 : seq_len + 1].mean(0)
            embeddings[pid] = token_repr
            cache.save(pid, token_repr)
    return embeddings


def cosine_distance(vec_a: torch.Tensor, vec_b: torch.Tensor) -> Optional[float]:
    if vec_a is None or vec_b is None:
        return None
    a = vec_a.unsqueeze(0)
    b = vec_b.unsqueeze(0)
    sim = torch.nn.functional.cosine_similarity(a, b)
    return float(1 - sim.item())


def load_go_annotations(
    go_path: Path, target_ids: Optional[Sequence[str]] = None
) -> Dict[str, Dict[str, List[str]]]:
    go_lookup: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    target_set = set(target_ids) if target_ids else None
    go_df = pd.read_csv(go_path, sep="\t")
    for row in go_df.itertuples(index=False):
        protein_full = getattr(row, "protein_id_full")
        aspect = getattr(row, "aspect")
        go_id = getattr(row, "go_id")
        go_name = getattr(row, "go_name")
        base_id = normalize_id(str(protein_full))
        if target_set is not None and base_id not in target_set:
            continue
        go_lookup[base_id][aspect].append(f"{go_id}:{go_name}")
    # Deduplicate lists
    for base_id in go_lookup:
        for aspect in go_lookup[base_id]:
            unique = sorted(set(go_lookup[base_id][aspect]))
            go_lookup[base_id][aspect] = unique
    return go_lookup


def collapse_go(go_lookup: Dict[str, Dict[str, List[str]]], protein_id: str, aspect: str) -> Tuple[str, int]:
    terms = go_lookup.get(protein_id, {}).get(aspect, [])
    return (";".join(terms), len(terms))


def main() -> None:
    args = parse_args()
    logger = init_logger()
    torch.set_num_threads(max(1, args.num_threads))
    torch.set_num_interop_threads(max(1, args.num_threads))
    ortholog_path = Path(args.ortholog_table)
    zebra_fasta = Path(args.zebrafish_fasta)
    aesc_fasta = Path(args.aescallii_fasta)
    go_path = Path(args.zebrafish_go)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Loading ortholog table %s", ortholog_path)
    ortho_df = pd.read_csv(ortholog_path)
    ortho_df = ortho_df.rename(
        columns={
            "Danio_rerio": "zebrafish_id",
            "Danio rerio": "zebrafish_id",
            "Danio aesculapii": "aescallii_id",
        }
    )
    if "zebrafish_id" not in ortho_df.columns or "aescallii_id" not in ortho_df.columns:
        raise SystemExit("Ortholog table must contain Danio_rerio and Danio aesculapii columns.")

    ortho_df["zebrafish_base"] = ortho_df["zebrafish_id"].astype(str).apply(normalize_id)
    ortho_df["aescallii_base"] = ortho_df["aescallii_id"].astype(str).apply(normalize_id)

    if args.row_offset or args.row_limit is not None:
        start = max(0, args.row_offset)
        end = start + args.row_limit if args.row_limit is not None else None
        logger.info("Subsetting ortholog rows: start=%s end=%s", start, end)
        ortho_df = ortho_df.iloc[start:end].reset_index(drop=True)

    required_zebra = sorted(set(ortho_df["zebrafish_base"]))
    required_aesc = sorted(set(ortho_df["aescallii_base"]))
    logger.info("Loading FASTA sequences (filtered to %d + %d IDs)", len(required_zebra), len(required_aesc))
    zebra_sequences = load_fasta_sequences(zebra_fasta, required_zebra)
    aesc_sequences = load_fasta_sequences(aesc_fasta, required_aesc)
    cache = EmbeddingCache(Path(args.cache_dir) if args.cache_dir else None)

    model, alphabet = load_model(args.model_name, args.device)

    zebra_embeddings = generate_embeddings(
        required_zebra,
        zebra_sequences,
        model,
        alphabet,
        args.device,
        args.batch_size,
        cache,
        args.max_seq_len,
    )
    aesc_embeddings = generate_embeddings(
        required_aesc,
        aesc_sequences,
        model,
        alphabet,
        args.device,
        args.batch_size,
        cache,
        args.max_seq_len,
    )

    logger.info("Loading GO annotations from %s (filtered)", go_path)
    go_lookup = load_go_annotations(go_path, required_zebra)

    records = []
    for row in ortho_df.itertuples(index=False, name=None):
        row_dict = dict(zip(ortho_df.columns, row))
        zebra_id = row_dict["zebrafish_id"]
        zebra_base = row_dict["zebrafish_base"]
        aesc_id = row_dict["aescallii_id"]
        aesc_base = row_dict["aescallii_base"]
        vec_a = zebra_embeddings.get(zebra_base)
        vec_b = aesc_embeddings.get(aesc_base)
        distance = cosine_distance(vec_a, vec_b) if vec_a is not None and vec_b is not None else None

        bp_terms, bp_count = collapse_go(go_lookup, zebra_base, "P")
        mf_terms, mf_count = collapse_go(go_lookup, zebra_base, "F")
        cc_terms, cc_count = collapse_go(go_lookup, zebra_base, "C")

        record = {
            "zebrafish_id": zebra_id,
            "aescallii_id": aesc_id,
            "zebrafish_base": zebra_base,
            "aescallii_base": aesc_base,
            "plm_cosine_distance": distance,
            "go_terms_bp": bp_terms,
            "go_terms_mf": mf_terms,
            "go_terms_cc": cc_terms,
            "go_count_bp": bp_count,
            "go_count_mf": mf_count,
            "go_count_cc": cc_count,
        }
        # Preserve other RBH columns (bitscores, coverage, flags, etc.).
        for col, value in row_dict.items():
            if col not in record:
                record[col] = value
        records.append(record)

    result_df = pd.DataFrame(records)
    mode = "a" if args.append else "w"
    header = not args.append or not output_path.exists()
    logger.info(
        "Writing %d rows to %s (mode=%s, header=%s)",
        len(result_df),
        output_path,
        mode,
        header,
    )
    result_df.to_csv(output_path, sep="\t", index=False, mode=mode, header=header)


if __name__ == "__main__":
    main()
