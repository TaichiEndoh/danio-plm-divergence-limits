#!/usr/bin/env python3
import argparse
import itertools
import math
import random
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import torch

try:
    import esm
except ImportError:
    esm = None


def parse_fasta(path: Path):
    seqs = {}
    current_id = None
    chunks = []
    if path.suffix == ".gz":
        import gzip
        fh = gzip.open(path, "rt")
    else:
        fh = path.open()
    with fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                raw_id = line[1:].split()[0]
                current_id = raw_id.split(".")[0]
                chunks = []
            else:
                chunks.append(line)
        if current_id is not None:
            seqs[current_id] = "".join(chunks)
    return seqs


def load_go_terms(bp_str):
    if not isinstance(bp_str, str) or not bp_str:
        return []
    terms = []
    for item in bp_str.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(":")
        if len(parts) >= 2 and parts[0] == "GO":
            go_id = f"{parts[0]}:{parts[1]}"
            terms.append(go_id)
        elif item.startswith("GO:"):
            go_id = item.split()[0]
            terms.append(go_id)
    return terms


def limit_go_terms(terms, max_terms=5):
    if not terms:
        return []
    uniq = sorted(set(terms))
    return uniq[:max_terms]


def build_go_graph(obo_path: Path):
    graph = {}
    current_id = None
    with obo_path.open() as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                current_id = None
                continue
            if line.startswith("id:"):
                current_id = line.split("id:")[1].strip()
                graph.setdefault(current_id, set())
                continue
            if current_id is None:
                continue
            if line.startswith("is_a:"):
                parent = line.split("is_a:")[1].split("!")[0].strip()
                graph.setdefault(current_id, set()).add(parent)
                graph.setdefault(parent, set()).add(current_id)
            if line.startswith("relationship: part_of"):
                parent = line.split("part_of")[1].split("!")[0].strip()
                graph.setdefault(current_id, set()).add(parent)
                graph.setdefault(parent, set()).add(current_id)
    return graph


def shortest_path(graph, start, goal):
    if start == goal:
        return 0
    if start not in graph or goal not in graph:
        return None
    visited = {start}
    queue = deque([(start, 0)])
    while queue:
        node, dist = queue.popleft()
        for nbr in graph.get(node, []):
            if nbr == goal:
                return dist + 1
            if nbr not in visited:
                visited.add(nbr)
                queue.append((nbr, dist + 1))
    return None


def min_go_distance(graph, terms_a, terms_b, cache):
    min_dist = None
    for ta in terms_a:
        for tb in terms_b:
            key = (ta, tb) if ta <= tb else (tb, ta)
            if key in cache:
                dist = cache[key]
            else:
                dist = shortest_path(graph, ta, tb)
                cache[key] = dist
            if dist is None:
                continue
            if min_dist is None or dist < min_dist:
                min_dist = dist
                if min_dist == 0:
                    return 0
    return min_dist


def embed_sequences(model, alphabet, seqs, device="cpu", batch_size=8):
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    embeds = {}
    with torch.no_grad():
        items = list(seqs.items())
        for i in range(0, len(items), batch_size):
            batch = items[i : i + batch_size]
            data = [(sid, s) for sid, s in batch]
            _, _, tokens = batch_converter(data)
            tokens = tokens.to(device)
            out = model(tokens, repr_layers=[model.num_layers])
            reps = out["representations"][model.num_layers]
            for j, (seq_id, seq) in enumerate(batch):
                emb = reps[j, 1 : len(seq) + 1].mean(0)
                embeds[seq_id] = emb.cpu()
    return embeds


def cosine_distance(v1, v2):
    num = torch.dot(v1, v2).item()
    den = (torch.norm(v1) * torch.norm(v2)).item()
    if den == 0:
        return None
    return 1.0 - num / den


def spearman_corr(x, y):
    if len(x) < 3:
        return None
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    mx = rx.mean()
    my = ry.mean()
    num = ((rx - mx) * (ry - my)).sum()
    den = math.sqrt(((rx - mx) ** 2).sum() * ((ry - my) ** 2).sum())
    if den == 0:
        return None
    return num / den


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plm-table", required=True)
    parser.add_argument("--obo", required=True)
    parser.add_argument("--zebrafish-fasta", required=True)
    parser.add_argument("--aescallii-fasta", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--out-fig", required=True)
    parser.add_argument("--embedding-cache", default="/Users/d/AI/plm_go_analysis/embeddings_cache")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--per-go", type=int, default=30)
    args = parser.parse_args()

    random.seed(args.seed)

    go_targets = [
        ("GO:0006955", "immune response"),
        ("GO:0002376", "immune system process"),
        ("GO:0045087", "innate immune response"),
        ("GO:0006470", "protein dephosphorylation"),
        ("GO:0016311", "dephosphorylation"),
        ("GO:0008150", "biological_process"),
        ("GO:0006281", "DNA repair"),
        ("GO:0006631", "fatty acid metabolic process"),
        ("GO:0006508", "proteolysis"),
        ("GO:0032259", "methylation"),
    ]

    df = pd.read_csv(args.plm_table, sep="\t")
    df["go_terms_bp_list"] = df["go_terms_bp"].apply(load_go_terms)

    go_graph = build_go_graph(Path(args.obo))

    # Use embedding cache presence to filter available proteins (faster than parsing full FASTA)
    cache_dir = Path(args.embedding_cache)
    cached_ids = {p.stem for p in cache_dir.glob("*.pt")}

    out_rows = []
    cache = {}

    # collect proteins per GO term (labelled by the GO term used for selection)
    # store full GO term lists per protein for distance calculation
    z_pool = {}
    a_pool = {}
    for go_id, go_name in go_targets:
        subset = df[df["go_terms_bp"].fillna("").str.contains(go_id)]
        if subset.empty:
            continue
        sample = subset.sample(n=min(args.per_go, len(subset)), random_state=args.seed)
        for z_id, a_id, z_terms in zip(
            sample["zebrafish_id"],
            sample["aescallii_id"],
            sample["go_terms_bp_list"],
        ):
            z_terms = limit_go_terms(z_terms)
            if z_id in cached_ids and z_id not in z_pool:
                z_pool[z_id] = (go_id, go_name, z_terms)
            if a_id in cached_ids and a_id not in a_pool:
                a_pool[a_id] = (go_id, go_name, z_terms)

    z_seqs = {pid: pid for pid in z_pool.keys()}
    a_seqs = {pid: pid for pid in a_pool.keys()}

    def load_cache(ids):
        cached = {}
        missing = []
        for pid in ids:
            path = cache_dir / f"{pid}.pt"
            if path.exists():
                cached[pid] = torch.load(path, map_location="cpu")
            else:
                missing.append(pid)
        return cached, missing

    z_cached, z_missing = load_cache(z_seqs.keys())
    a_cached, a_missing = load_cache(a_seqs.keys())

    if z_missing or a_missing:
        print("Warning: missing embeddings detected; skipping those proteins to avoid model download.")

    z_emb = z_cached
    a_emb = a_cached

    for species, emb, pool in [
        ("zebrafish", z_emb, z_pool),
        ("aescallii", a_emb, a_pool),
    ]:
        ids = list(emb.keys())
        for id_a, id_b in itertools.combinations(ids, 2):
            go_a_label, _, go_a_terms = pool[id_a]
            go_b_label, _, go_b_terms = pool[id_b]
            go_dist = min_go_distance(go_graph, go_a_terms, go_b_terms, cache)
            if go_dist is None:
                continue
            esm_dist = cosine_distance(emb[id_a], emb[id_b])
            if esm_dist is None:
                continue
            out_rows.append({
                "species": species,
                "go_id": f"{go_a_label}|{go_b_label}",
                "go_name": None,
                "protein_a": id_a,
                "protein_b": id_b,
                "go_graph_distance": int(go_dist),
                "esm_distance": esm_dist,
            })

    out_df = pd.DataFrame(out_rows)
    out_path = Path(args.out_tsv)
    out_df.to_csv(out_path, sep="\t", index=False)

    fig_path = Path(args.out_fig)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    dist_levels = sorted(out_df["go_graph_distance"].dropna().unique())
    data = [out_df[out_df["go_graph_distance"] == d]["esm_distance"].values for d in dist_levels]
    ax.boxplot(data, labels=[str(d) for d in dist_levels], showfliers=False)
    ax.set_xlabel("GO graph distance (shortest path, integer)")
    ax.set_ylabel("ESM-2 cosine distance")
    ax.set_title("GO graph distance vs ESM-2 distance")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200)

    # print correlations
    for species, group in out_df.groupby("species"):
        corr = spearman_corr(group["go_graph_distance"].to_numpy(), group["esm_distance"].to_numpy())
        print(f"{species}\tN={len(group)}\tSpearman={corr}")


if __name__ == "__main__":
    main()
