#!/usr/bin/env python3
"""
Fetch the ESM-3 subset's protein sequences from NCBI (batched efetch), so the ESM-2 /
ESM-3 per-pair embedding runs (paper Fig 2/3) have their inputs ready.

Reads reports/esm3_subset_pairs.csv, collects the unique zebrafish (D. rerio) and
aesculapii RefSeq accessions, and writes two FASTAs. Version suffixes are stripped from
the FASTA ids so they match the pairs table.

Usage:
    python analysis/fetch_subset_sequences.py \\
        --pairs reports/esm3_subset_pairs.csv \\
        --out-a data/reference/sequences/subset/subset_rerio.fasta \\
        --out-b data/reference/sequences/subset/subset_aesculapii.fasta
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch_batch(ids, tries=5):
    data = urlencode({"db": "protein", "id": ",".join(ids),
                      "rettype": "fasta", "retmode": "text"}).encode()
    for t in range(1, tries + 1):
        try:
            return urlopen(EFETCH, data=data, timeout=120).read().decode()
        except Exception as e:
            print(f"    batch retry {t}: {type(e).__name__}: {str(e)[:80]}", flush=True)
            time.sleep(2 * t)
    raise RuntimeError("efetch batch failed after retries")


def normalize(fasta_text: str) -> dict:
    """Return {id_without_version: fasta_record_text}."""
    out = {}
    cur_id, buf = None, []
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if cur_id:
                out[cur_id] = "".join(buf)
            cur_id = line[1:].split()[0].split(".")[0]
            buf = [line + "\n"]
        elif cur_id:
            buf.append(line + "\n")
    if cur_id:
        out[cur_id] = "".join(buf)
    return out


def fetch_all(ids, batch=200):
    got = {}
    ids = list(ids)
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        got.update(normalize(fetch_batch(chunk)))
        print(f"  fetched {min(i+batch,len(ids))}/{len(ids)} (have {len(got)})", flush=True)
        time.sleep(0.5)  # NCBI politeness (<=3 req/s without an API key)
    return got


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default="reports/esm3_subset_pairs.csv")
    ap.add_argument("--out-a", default="data/reference/sequences/subset/subset_rerio.fasta")
    ap.add_argument("--out-b", default="data/reference/sequences/subset/subset_aesculapii.fasta")
    args = ap.parse_args()

    d = pd.read_csv(args.pairs)
    z = sorted(d["zebrafish_id"].astype(str).str.split(".").str[0].unique())
    a = sorted(d["aescallii_id"].astype(str).str.split(".").str[0].unique())
    print(f"fetching {len(z)} zebrafish + {len(a)} aesculapii sequences...")

    for label, ids, out in [("zebrafish", z, args.out_a), ("aesculapii", a, args.out_b)]:
        print(f"[{label}]")
        recs = fetch_all(ids)
        missing = [i for i in ids if i not in recs]
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        Path(out).write_text("".join(recs[i] for i in ids if i in recs))
        print(f"  wrote {len(recs)}/{len(ids)} -> {out}  (missing {len(missing)})")
        if missing:
            print("  missing ids:", missing[:20], "..." if len(missing) > 20 else "")
    print("FETCH_DONE")


if __name__ == "__main__":
    main()
