#!/usr/bin/env python3
"""
Sequence-level identity at the AHR2 ligand pocket and at the annotated domain intervals.

"Conserved at the pocket" can mean two different things: identical amino acids, or a small
embedding distance. The manuscript needs the first stated explicitly, from a committed output
rather than an ad-hoc check, so that the second is not read as more than it is.

The point matters: at these positions the embedding result is not a model inference at all —
the residues are simply identical, which anyone can verify by counting.

Writes reports/ahr2_pocket_sequence_identity.csv
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import Align

FASTA = Path("data/reference/sequences/target_domain_sequences.fasta")
ZEB, AESC = "NP_571339", "XP_056303610"

# Takeda ligand-pocket residues, D. rerio AHR2 numbering (as used throughout the paper)
POCKET = [281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356, 370, 378, 380, 381]

# Domain intervals as reported in Results §3.5 (InterPro, transferred to our numbering)
INTERVALS = [
    ("bHLH, DNA binding", 25, 86),
    ("PAS-A", 119, 231),
    ("PAS domain (wider)", 277, 385),
    ("PAS fold 3, ligand-binding", 308, 385),
    ("C-terminal transactivation domain", 386, 1027),
]


def read_fasta(path: Path) -> dict[str, str]:
    seqs, cur, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0].split(".")[0], []
        elif cur:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def main() -> None:
    seqs = read_fasta(FASTA)
    z, a = seqs[ZEB], seqs[AESC]

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 1, 0
    aligner.open_gap_score, aligner.extend_gap_score = -2, -0.5
    aln = aligner.align(z, a)[0]
    top, bottom = str(aln[0]), str(aln[1])

    # map D. rerio residue number -> (rerio aa, aesculapii aa)
    pos, pairs = 0, {}
    for x, y in zip(top, bottom):
        if x != "-":
            pos += 1
            pairs[pos] = (x, y)

    whole = 100.0 * sum(1 for x, y in zip(top, bottom) if x == y) / len(top)
    print(f"AHR2 {ZEB} vs {AESC}: whole-protein identity {whole:.1f}% (aligned length {len(top)})\n")

    rows = []
    hits = [p for p in POCKET if p in pairs and pairs[p][0] == pairs[p][1]]
    rows.append(("curated ligand pocket (Takeda, 18 residues)", "", "",
                 len(POCKET), len(hits), 100.0 * len(hits) / len(POCKET)))
    print(f"  curated ligand pocket: {len(hits)}/{len(POCKET)} identical "
          f"({100*len(hits)/len(POCKET):.1f}%)")
    mismatched = [(p, pairs[p]) for p in POCKET if p in pairs and pairs[p][0] != pairs[p][1]]
    if mismatched:
        print(f"    mismatches: {mismatched}")

    for label, lo, hi in INTERVALS:
        span = [p for p in range(lo, hi + 1) if p in pairs]
        same = sum(1 for p in span if pairs[p][0] == pairs[p][1])
        pct = 100.0 * same / len(span) if span else float("nan")
        rows.append((label, lo, hi, len(span), same, pct))
        print(f"  {label} ({lo}-{hi}): {same}/{len(span)} identical ({pct:.1f}%)")

    out = Path("reports/ahr2_pocket_sequence_identity.csv")
    df = pd.DataFrame(rows, columns=["region", "start", "end", "n_residues",
                                     "n_identical", "percent_identical"])
    df.attrs["whole_protein_identity"] = whole
    df.to_csv(out, index=False)
    with out.open("a") as fh:
        fh.write(f"# whole-protein identity,{whole:.4f}\n")
        fh.write(f"# zebrafish={ZEB}, aesculapii={AESC}, aligned_length={len(top)}\n")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
