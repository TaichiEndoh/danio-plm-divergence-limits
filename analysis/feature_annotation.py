#!/usr/bin/env python3
"""
Priority 3 & 4 (Prof. Endo): map local ortholog divergence onto biological features.

AHR2  — domain architecture from InterPro (bHLH / PAS-A / PAS-B ligand-binding /
        C-terminal transactivation domain) + a curated short-linear-motif (SLiM) scan
        of the divergent C-terminal region, comparing D. rerio vs D. aesculapii to flag
        motif gain/loss.
Kcnj13 — UniProt features (transmembrane helices, pore/selectivity region, cytoplasmic
        domains) with the divergence hotspots mapped onto them.

Feature coordinates are fetched live (UniProt REST / InterPro API) and aligned to our
zebrafish reference sequence, so residue numbering matches the divergence CSVs. Live
coordinates fall back to the values fetched on 2026-08-23 if the network is unavailable.
"""
from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align

OI = {"blue": "#0072B2", "green": "#009E73", "orange": "#E69F00", "vermillion": "#D55E00", "grey": "#999999"}

# --- fallback feature coordinates (UniProt/InterPro numbering), fetched 2026-08-23 ---
AHR2_DOMAINS_FALLBACK = [
    ("bHLH (DNA binding)", 25, 86, "core"),
    ("PAS-A", 115, 185, "core"),
    ("PAS-B (ligand pocket)", 277, 385, "core"),
    ("PAC motif", 353, 394, "core"),
    ("C-terminal TAD", 395, 1027, "periphery"),
]
KCNJ13_FEATURES_FALLBACK = [
    ("TM1 (outer helix)", 66, 89, "structured"),
    ("Pore + selectivity filter", 90, 145, "structured"),
    ("TM2 (inner helix)", 146, 170, "structured"),
    ("Kir C-terminal cytoplasmic", 182, 319, "structured"),
]
AHR2_POCKET = [281, 283, 289, 291, 303, 307, 311, 314, 322, 324, 335, 340, 348, 356, 370, 378, 380, 381]


def fetch(url, timeout=40):
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return json.load(r)
    except Exception as e:  # noqa: BLE001
        print(f"  [warn] fetch failed ({e}); using fallback")
        return None


def read_fasta(p):
    s, cur, buf = {}, None, []
    for line in Path(p).read_text().splitlines():
        if line.startswith(">"):
            if cur:
                s[cur] = "".join(buf)
            cur, buf = line[1:].split()[0].split(".")[0], []
        elif cur:
            buf.append(line.strip())
    if cur:
        s[cur] = "".join(buf)
    return s


def uniprot_seq(acc):
    d = fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.json?fields=sequence")
    return d["sequence"]["value"] if d else None


def map_ranges_to_ref(ref_seq, feat_seq, ranges):
    """Map (label,start,end,cls) ranges in feat_seq numbering onto ref_seq numbering
    via global alignment. Returns dict label -> set(ref residue numbers, 1-based)."""
    # Key by (label, start, end): UniProt returns several features with identical labels
    # (e.g. two "Transmembrane: Helical" helices), and keying by label alone silently drops
    # all but the last of them.
    if ref_seq == feat_seq:
        return {(lab, s, e): set(range(s, e + 1)) for lab, s, e, _ in ranges}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 1, 0
    aligner.open_gap_score, aligner.extend_gap_score = -1, -0.5
    aln = aligner.align(ref_seq, feat_seq)[0]
    # build feat_pos(1-based) -> ref_pos(1-based)
    f2r = {}
    for (rs, re_), (fs, fe) in zip(aln.aligned[0], aln.aligned[1]):
        for off in range(re_ - rs):
            f2r[fs + off + 1] = rs + off + 1
    out = {}
    for lab, s, e, _ in ranges:
        out[(lab, s, e)] = {f2r[p] for p in range(s, e + 1) if p in f2r}
    return out


def summarize(name, div, feat_map, order):
    print(f"\n== {name}: mean divergence per feature ==")
    rows = []
    for key in order:
        lab = key[0] if isinstance(key, tuple) else key
        res = feat_map.get(key, set())
        sub = div[div["residue"].isin(res)]
        if len(sub):
            rows.append((f"{lab} ({key[1]}-{key[2]})" if isinstance(key, tuple) else lab, len(sub), sub["divergence"].mean(), sub["divergence"].median()))
            print(f"  {lab:28s} n={len(sub):4d}  mean={sub['divergence'].mean():.5f}  median={sub['divergence'].median():.5f}")
    # residues not in any annotated feature
    annotated = set().union(*feat_map.values()) if feat_map else set()
    unann = div[~div["residue"].isin(annotated)]
    if len(unann):
        rows.append(("(unannotated / linker)", len(unann), unann["divergence"].mean(), unann["divergence"].median()))
        print(f"  {'(unannotated / linker)':28s} n={len(unann):4d}  mean={unann['divergence'].mean():.5f}  median={unann['divergence'].median():.5f}")
    return pd.DataFrame(rows, columns=["feature", "n", "mean", "median"])


# ---- curated SLiM classes (ELM-style regex; exploratory) ----
SLIMS = {
    "NLS (basic cluster)": r"[KR]{4,}",
    "bipartite NLS": r"[KR]{2}.{9,12}[KR]{3,}",
    "CK2 phospho": r"[ST]..[DE]",
    "MAPK phospho (S/T-P)": r"[ST]P",
    "SUMO (VKxE)": r"[VILMAFP]K.E",
    "phosphodegron (DSGxxS)": r"DSG..S",
}


def slim_scan_gain_loss(name, zseq, aseq, div, hi_thresh_pct=90):
    """Within high-divergence windows, report SLiM motifs present in one ortholog but
    disrupted in the other (motif gain/loss), aligning rerio<->aesc by the divergence CSV."""
    thr = np.percentile(div["divergence"], hi_thresh_pct)
    hi = set(div[div["divergence"] >= thr]["residue"])
    # build rerio_pos -> aesc_pos from the divergence CSV
    r2a = dict(zip(div["residue"], div["aescallii_residue"]))
    print(f"\n== {name}: SLiM gain/loss in high-divergence regions (top {100-hi_thresh_pct}% , thr={thr:.4f}) ==")
    events = []
    for lab, pat in SLIMS.items():
        rx = re.compile(pat)
        zhits = {m.start() + 1 for m in rx.finditer(zseq)}  # 1-based motif start
        # a motif "in a hotspot" if its span overlaps a high-divergence residue
        zspan = [(m.start() + 1, m.end()) for m in rx.finditer(zseq)]
        for s, e in zspan:
            if not any(p in hi for p in range(s, e + 1)):
                continue
            # corresponding aesc window
            a_positions = [int(r2a[p]) for p in range(s, e + 1) if p in r2a and not pd.isna(r2a.get(p))]
            if not a_positions:
                continue
            a_lo, a_hi = min(a_positions) - 2, max(a_positions) + 2
            awin = aseq[max(0, a_lo - 1):a_hi]
            present_in_aesc = bool(rx.search(awin))
            zwin = zseq[s - 1:e]
            if not present_in_aesc:
                events.append((lab, s, e, zwin, awin, "LOST in aesc"))
    if events:
        for lab, s, e, zw, aw, tag in events[:40]:
            print(f"  {lab:22s} rerio {s}-{e} '{zw}'  aesc~'{aw}'  -> {tag}")
    else:
        print("  (no motif gain/loss events in high-divergence windows)")
    return events


def bar(ax, r, title, palette):
    colors = [palette.get(f, OI["grey"]) for f in r["feature"]]
    ax.bar(range(len(r)), r["mean"], color=colors, width=0.66)
    for i, (v, n) in enumerate(zip(r["mean"], r["n"])):
        ax.text(i, v, f"{v:.4f}\n(n={n})", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(range(len(r)))
    ax.set_xticklabels(r["feature"], rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("mean local divergence")
    ax.set_title(title, fontsize=11, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    # ---------- AHR2 (priority 3) ----------
    ref = read_fasta("data/reference/sequences/target_domain_sequences.fasta")
    ahr2_z = ref["NP_571339"]
    ahr2_a = ref["XP_056303610"]
    up = uniprot_seq("A0A8M1PB00") or ahr2_z
    # domains from InterPro (live), fallback to constants
    doms = []
    ip = fetch("https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/A0A8M1PB00/?page_size=100")
    if ip:
        want = {"IPR011598": "bHLH (DNA binding)", "IPR013767": "PAS-A", "IPR013655": "PAS-B (ligand pocket)"}
        for rr in ip.get("results", []):
            acc = rr["metadata"]["accession"]
            if acc in want:
                for pr in rr.get("proteins", []):
                    for l in pr.get("entry_protein_locations", []):
                        for f in l.get("fragments", []):
                            doms.append((want[acc], f["start"], f["end"], "core"))
        # add C-terminal TAD as everything past the last PAS/PAC end
        last = max((e for _, _, e, _ in doms), default=394)
        doms.append(("C-terminal TAD", last + 1, len(up), "periphery"))
    if not doms:
        doms = AHR2_DOMAINS_FALLBACK
    ahr2_div = pd.read_csv("reports/ahr2_residue_divergence_fixed.csv")
    fmap = map_ranges_to_ref(ahr2_z, up, doms)
    order_a = [(lab, s, e) for lab, s, e, c in sorted(doms, key=lambda x: x[1])]
    r_ahr2 = summarize("AHR2", ahr2_div, fmap, order_a)
    slim_events = slim_scan_gain_loss("AHR2", ahr2_z, ahr2_a, ahr2_div)

    # ---------- Kcnj13 (priority 4) ----------
    kp = read_fasta("data/reference/sequences/kcnj13_pair.fasta")
    kz = kp["NP_001039014"]
    feats = []
    kj = fetch("https://rest.uniprot.org/uniprotkb/A0AC58HB01.json?fields=ft_transmem,ft_domain,ft_intramem,sequence")
    kup = kj["sequence"]["value"] if kj else kz
    if kj and kj.get("features"):
        for f in kj["features"]:
            lab = f"{f['type']}: {f.get('description','')}".strip(": ")
            feats.append((lab, f["location"]["start"]["value"], f["location"]["end"]["value"], "structured"))
    if not feats:
        feats = KCNJ13_FEATURES_FALLBACK
        kup = kz
    kdiv = pd.read_csv("reports/kcnj13_residue_divergence_fixed.csv")
    kmap = map_ranges_to_ref(kz, kup, feats)
    order_k = [(lab, s, e) for lab, s, e, c in sorted(feats, key=lambda x: x[1])]
    r_kcnj = summarize("Kcnj13", kdiv, kmap, order_k)
    # where are the top hotspots?
    print("\n== Kcnj13: top-10 divergence hotspots and their feature ==")
    res2feat = {}
    for lab, res in kmap.items():
        for p in res:
            res2feat.setdefault(p, []).append(lab)
    for row in kdiv.nlargest(10, "divergence").itertuples(index=False):
        print(f"  res {int(row.residue):3d}  div={row.divergence:.4f}  -> {res2feat.get(int(row.residue), ['(unannotated)'])}")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    bar(axes[0], r_ahr2, "AHR2: divergence by domain (bHLH/PAS conserved, TAD divergent)",
        {"bHLH (DNA binding)": OI["blue"], "PAS-A": OI["blue"], "PAS-B (ligand pocket)": OI["blue"],
         "PAC motif": OI["green"], "C-terminal TAD": OI["vermillion"], "(unannotated / linker)": OI["grey"]})
    bar(axes[1], r_kcnj, "Kcnj13: divergence by membrane-channel feature",
        {f[0]: OI["green"] for f in feats} | {"(unannotated / linker)": OI["grey"]})
    fig.suptitle("Divergence mapped onto biological features (UniProt / InterPro)", fontsize=13, weight="bold")
    fig.tight_layout()
    out = Path("reports/figures/feature_annotation.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nWrote {out}")

    # ---------- machine-readable report ----------
    r_ahr2.insert(0, "protein", "AHR2")
    r_kcnj.insert(0, "protein", "Kcnj13")
    pd.concat([r_ahr2, r_kcnj]).to_csv("reports/feature_divergence.csv", index=False)
    print("Wrote reports/feature_divergence.csv")


if __name__ == "__main__":
    main()
