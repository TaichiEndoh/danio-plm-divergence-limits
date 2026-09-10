#!/usr/bin/env python3
"""
Strengthen the ranking-reproducibility result (paper Fig 2/3, Prof. Endo's §4):
bootstrap CI for Spearman, top-k overlap vs a random baseline, stratification by
sequence length and identity, ESM distance vs plain sequence identity, and the
ESM-2/ESM-3 discordant pairs.

Sequence identity per pair is computed by global alignment (Biopython PairwiseAligner,
sequences capped at 2046 to match the embedding runs) and cached to
reports/subset_identity.csv.
"""
from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OI = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00", "grey": "#999999"}
CAP = 2046


def read_fasta(p):
    s, cur, buf = {}, None, []
    for line in Path(p).read_text().splitlines():
        if line.startswith(">"):
            if cur: s[cur] = "".join(buf)
            cur, buf = line[1:].split()[0].split(".")[0], []
        elif cur: buf.append(line.strip())
    if cur: s[cur] = "".join(buf)
    return s


def load(path, col):
    d = pd.read_csv(path)
    d["z"] = d["zebrafish_id"].astype(str).str.split(".").str[0]
    d["a"] = d["aescallii_id"].astype(str).str.split(".").str[0]
    return d[["z", "a", col]]


def compute_identity(pairs, zseq, aseq, out):
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 1, 0
    aligner.open_gap_score, aligner.extend_gap_score = -1, -0.5
    rows, t0 = [], time.time()
    for i, r in enumerate(pairs.itertuples(index=False), 1):
        s1, s2 = zseq.get(r.z), aseq.get(r.a)
        if not s1 or not s2:
            rows.append((r.z, r.a, np.nan, np.nan)); continue
        s1, s2 = s1[:CAP], s2[:CAP]
        aln = aligner.align(s1, s2)[0]
        c = aln.counts()
        ident = 100.0 * c.identities / aln.length if aln.length else np.nan
        rows.append((r.z, r.a, ident, max(len(s1), len(s2))))
        if i % 400 == 0:
            print(f"  identity {i}/{len(pairs)} ({(time.time()-t0)/i:.3f}s/pair)", flush=True)
    df = pd.DataFrame(rows, columns=["z", "a", "identity", "maxlen"])
    df.to_csv(out, index=False)
    return df


def boot_spearman(x, y, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    idx = np.arange(len(x))
    rs = [spearmanr(x[s], y[s]).statistic for s in (rng.choice(idx, len(idx), replace=True) for _ in range(n))]
    return float(np.percentile(rs, 2.5)), float(np.percentile(rs, 97.5))


def topk_overlap(df, c1, c2, frac):
    k = max(1, int(round(frac * len(df))))
    s1 = set(map(tuple, df.nlargest(k, c1)[["z", "a"]].values))
    s2 = set(map(tuple, df.nlargest(k, c2)[["z", "a"]].values))
    return len(s1 & s2) / k, k, frac  # observed, k, random-expectation=frac


def main() -> None:
    esm3 = load("reports/esm3_subset_esm3_1p4b.csv", "esm3_cosine_distance").rename(columns={"esm3_cosine_distance": "esm3"})
    m = esm3.copy()
    for name, p in [("d8M", "reports/esm3_subset_esm2_8m.csv"), ("d35M", "reports/esm3_subset_esm2_35m.csv"),
                    ("d150M", "reports/esm3_subset_esm2_150m.csv"), ("d650M", "reports/esm3_subset_esm2_650m.csv")]:
        if Path(p).is_file():
            m = m.merge(load(p, "esm2_cosine_distance").rename(columns={"esm2_cosine_distance": name}), on=["z", "a"])

    # sequence identity (cached)
    idf_path = Path("reports/subset_identity.csv")
    if idf_path.is_file():
        idf = pd.read_csv(idf_path)
    else:
        print("computing per-pair sequence identity (global alignment)...", flush=True)
        idf = compute_identity(m[["z", "a"]], read_fasta("data/reference/sequences/subset/subset_rerio.fasta"),
                               read_fasta("data/reference/sequences/subset/subset_aesculapii.fasta"), idf_path)
    m = m.merge(idf, on=["z", "a"], how="left")
    m = m.dropna(subset=["d8M", "esm3"]).reset_index(drop=True)
    print(f"\nn pairs: {len(m)}")

    # 1) bootstrap Spearman CI (8M vs ESM-3)
    rho, p = spearmanr(m["d8M"], m["esm3"])
    lo, hi = boot_spearman(m["d8M"].to_numpy(), m["esm3"].to_numpy())
    print(f"\nSpearman(ESM-2 8M, ESM-3) = {rho:.3f}  95%CI [{lo:.3f}, {hi:.3f}]  p={p:.1e}")

    # 2) top-k overlap vs random baseline
    print("\ntop-k overlap with ESM-3 (observed vs random expectation):")
    for frac in (0.01, 0.05, 0.10):
        for size in [c for c in ["d8M", "d35M", "d150M", "d650M"] if c in m]:
            o, k, rnd = topk_overlap(m, size, "esm3", frac)
            print(f"  {size:5s} top-{int(frac*100):>2d}% (k={k}): {o*100:5.1f}%  (random ~{rnd*100:.1f}%, {o/rnd:.1f}x)")

    # 3) ESM distance vs sequence identity
    m["seqdiv"] = 100.0 - m["identity"]
    rho_id, _ = spearmanr(m["d8M"], m["seqdiv"])
    rho_id3, _ = spearmanr(m["esm3"], m["seqdiv"])
    print(f"\nSpearman(ESM-2 8M, sequence divergence=100-identity) = {rho_id:.3f}")
    print(f"Spearman(ESM-3,    sequence divergence)               = {rho_id3:.3f}")
    print(f"  (median identity = {m['identity'].median():.1f}%)")

    # 4) stratify Spearman(8M vs ESM3) by length and identity terciles
    for var, lbl in [("maxlen", "sequence length"), ("identity", "sequence identity")]:
        q = pd.qcut(m[var], 3, labels=["low", "mid", "high"], duplicates="drop")
        print(f"\nSpearman(8M vs ESM-3) by {lbl} tercile:")
        for g in q.cat.categories:
            sub = m[q == g]
            r, _ = spearmanr(sub["d8M"], sub["esm3"])
            print(f"  {g:4s} (n={len(sub)}, {var} {sub[var].min():.0f}-{sub[var].max():.0f}): rho={r:.3f}")

    # 5) discordant pairs (8M vs ESM-3 rank)
    m["r8"] = m["d8M"].rank(); m["r3"] = m["esm3"].rank()
    m["rankdiff"] = (m["r8"] - m["r3"]).abs()
    disc = m.nlargest(10, "rankdiff")[["z", "a", "d8M", "esm3", "identity", "maxlen", "rankdiff"]]
    print("\ntop-10 most ESM-2/ESM-3 discordant pairs:")
    print(disc.to_string(index=False))
    print(f"\ndiscordant set median identity {disc['identity'].median():.1f}% vs overall {m['identity'].median():.1f}%; "
          f"median len {disc['maxlen'].median():.0f} vs {m['maxlen'].median():.0f}")

    # figure: top-k overlap vs random + ESM-vs-identity
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    sizes = [c for c in ["d8M", "d35M", "d150M", "d650M"] if c in m]
    labels = {"d8M": "8M", "d35M": "35M", "d150M": "150M", "d650M": "650M"}
    xi = np.arange(len(sizes)); w = 0.25
    for j, frac in enumerate((0.01, 0.05, 0.10)):
        vals = [topk_overlap(m, s, "esm3", frac)[0] * 100 for s in sizes]
        ax[0].bar(xi + (j - 1) * w, vals, w, label=f"top-{int(frac*100)}%",
                  color=[OI["vermillion"], OI["orange"], OI["green"]][j])
    # Two baselines: chance (dotted) and, far more demanding, ranking by sequence identity
    # alone (solid horizontal lines). The identity baseline is the one a practitioner would
    # actually have to beat.
    for frac, c in zip((0.01, 0.05, 0.10), [OI["vermillion"], OI["orange"], OI["green"]]):
        ax[0].axhline(frac * 100, ls=":", lw=0.9, color=c, alpha=0.75)
        idv = topk_overlap(m, "seqdiv", "esm3", frac)[0] * 100
        ax[0].axhline(idv, ls="--", lw=1.6, color=c)
        # label sits in the right margin, clear of the bars
        ax[0].text(len(sizes) - 0.30, idv, f"identity {idv:.0f}%", fontsize=7.8, color=c,
                   va="center", ha="left",
                   bbox=dict(fc="white", ec="none", alpha=.85, pad=1.2))
    ax[0].set_xticks(xi); ax[0].set_xticklabels([f"ESM-2 {labels[s]}" for s in sizes])
    ax[0].set_ylabel("overlap with ESM-3 top-k (%)")
    ax[0].set_title("Top-k overlap against two baselines\n"
                    "dashed = ranking by sequence identity;  dotted = chance",
                    fontsize=10.5, weight="bold", pad=26)
    ax[0].set_xlim(-0.55, len(sizes) + 0.45)
    ax[0].set_ylim(0, 95)
    # legend above the axes so it cannot collide with bars or baseline labels
    ax[0].legend(frameon=False, fontsize=9, ncol=3, loc="lower left",
                 bbox_to_anchor=(0.0, 1.005))
    ax[0].spines[["top", "right"]].set_visible(False)

    ax[1].scatter(m["seqdiv"], np.clip(m["d8M"], 1e-5, None), s=8, alpha=0.35, c=OI["blue"], edgecolors="none")
    ax[1].set_yscale("log")
    ax[1].set_xlabel("sequence divergence (100 − % identity)")
    ax[1].set_ylabel("ESM-2 (8M) cosine distance (log)")
    ax[1].set_title(f"ESM distance vs sequence identity (ρ={rho_id:.2f})", fontsize=11, weight="bold")
    ax[1].spines[["top", "right"]].set_visible(False)
    out = Path("reports/figures/fig_ranking_stats.png")
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
