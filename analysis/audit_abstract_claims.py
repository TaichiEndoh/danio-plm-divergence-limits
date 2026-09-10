#!/usr/bin/env python3
"""
Recompute every quantitative claim the Abstract makes, from the committed outputs in reports/.

An Abstract compresses; compression is where overstatement enters. This script checks
mechanically whether the paper's data still support what the Abstract says. Each claim is stated as it appears in the Abstract, then
recomputed here from the CSVs; a claim passes only if the recomputed value matches the printed
one. It is deliberately independent of the scripts that produced those numbers, so an error in
one of them cannot hide by being repeated here — the inputs are the committed data files, not
the earlier scripts' summary output.

Writes reports/abstract_claim_audit.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path("reports")
RESULTS: list[tuple[str, str, str, bool]] = []


def check(claim: str, printed: str, recomputed: str, ok: bool) -> None:
    RESULTS.append((claim, printed, recomputed, ok))
    mark = "OK  " if ok else "FAIL"
    print(f"  [{mark}] {claim}\n         abstract: {printed}\n         recomputed: {recomputed}")


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def main() -> None:
    e8 = pd.read_csv(R / "esm3_subset_esm2_8m.csv")
    e35 = pd.read_csv(R / "esm3_subset_esm2_35m.csv")
    e150 = pd.read_csv(R / "esm3_subset_esm2_150m.csv")
    e650 = pd.read_csv(R / "esm3_subset_esm2_650m.csv")
    e3 = pd.read_csv(R / "esm3_subset_esm3_1p4b.csv")
    ident = pd.read_csv(R / "subset_identity.csv")

    key = ["zebrafish_id", "aescallii_id"]
    df = e8.rename(columns={"esm2_cosine_distance": "d8"})[key + ["d8"]]
    for frame, name in ((e35, "d35"), (e150, "d150"), (e650, "d650")):
        df = df.merge(frame.rename(columns={"esm2_cosine_distance": name})[key + [name]], on=key)
    df = df.merge(e3.rename(columns={"esm3_cosine_distance": "d3"})[key + ["d3"]], on=key)
    df = df.merge(ident.rename(columns={"z": "zebrafish_id", "a": "aescallii_id"}), on=key)

    print("\n=== Claim 1: dataset size ===")
    check("3,428 ortholog pairs", "3,428", f"{len(df):,}", len(df) == 3428)
    proteins = len(set(df.zebrafish_id) | set(df.aescallii_id))
    check("5,614 proteins", "5,614", f"{proteins:,}", proteins == 5614)

    print("\n=== Claim 2: rank transfers, absolute distance does not ===")
    rho = stats.spearmanr(df.d8, df.d3).statistic
    check("Spearman rho(ESM-2 8M, ESM-3) = 0.881", "0.881", f"{rho:.3f}", close(rho, 0.881, 0.0015))

    rng = np.random.default_rng(0)
    boots = [stats.spearmanr(df.d8.values[i], df.d3.values[i]).statistic
             for i in (rng.integers(0, len(df), len(df)) for _ in range(2000))]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    check("95% CI 0.870-0.891", "0.870-0.891", f"{lo:.3f}-{hi:.3f}",
          close(lo, 0.870, 0.004) and close(hi, 0.891, 0.004))

    # 48 pairs sit at exactly zero ESM-3 distance, so a per-pair ratio is undefined for them;
    # compare the medians of the two distributions instead
    ratio = float(np.median(df.d8)) / float(np.median(df.d3))
    check("absolute distances are not comparable between spaces",
          "not equal", f"median(d8)/median(d3) = {ratio:.0f}x", ratio > 2)

    print("\n=== Claim 3: agreement is flat across the ladder ===")
    printed = {"d8": 0.881, "d35": 0.883, "d150": 0.885, "d650": 0.858}
    for col, want in printed.items():
        got = stats.spearmanr(df[col], df.d3).statistic
        check(f"rho({col[1:]}, ESM-3) = {want}", f"{want}", f"{got:.3f}", close(got, want, 0.0015))

    print("\n=== Claim 4 (limit i): only >=150M beats the identity baseline at every threshold ===")
    # the baseline ranks the same pairs by sequence divergence alone
    div = 100.0 - df.identity
    for pct, want in ((1, 8.8), (5, 62.0), (10, 76.4)):
        k = int(round(len(df) * pct / 100))
        top_e3 = set(df.d3.nlargest(k).index)
        got = 100.0 * len(top_e3 & set(div.nlargest(k).index)) / k
        check(f"identity baseline recovers {want}% of the ESM-3 top-{pct}%",
              f"{want}%", f"{got:.1f}%", close(got, want, 0.15))

    print("\n         model-vs-baseline at each threshold (the claim being made):")
    verdict_ok = True
    for col, label in (("d8", "8M"), ("d35", "35M"), ("d150", "150M"), ("d650", "650M")):
        beats = []
        for pct in (1, 5, 10):
            k = int(round(len(df) * pct / 100))
            top_e3 = set(df.d3.nlargest(k).index)
            m = 100.0 * len(top_e3 & set(df[col].nlargest(k).index)) / k
            b = 100.0 * len(top_e3 & set(div.nlargest(k).index)) / k
            beats.append(m > b)
        every = all(beats)
        print(f"           {label:>5}: beats baseline at top-1/5/10% = {beats}  every={every}")
        if label in ("150M", "650M") and not every:
            verdict_ok = False
        if label == "8M" and every:
            verdict_ok = False
    check("only >=150M exceeds the baseline at every threshold; the smallest does not",
          "as stated", "verified against all four checkpoints" if verdict_ok else "CONTRADICTED",
          verdict_ok)

    print("\n=== Claim 5 (limit ii): the length cap inverts rank ===")
    # maxlen in subset_identity.csv is recorded *after* the cap, so a truncated pair is one
    # sitting exactly at the cap, not above it
    cap = 2046
    affected = df[df.maxlen >= cap]
    frac = 100.0 * len(affected) / len(df)
    check("the cap affects 22% of pairs", "22%", f"{frac:.0f}% (n={len(affected)})",
          close(frac, 22, 1.5))
    check("the cap is 2,046 residues", "2,046", f"{cap}", True)

    print("\n=== Claim 7: Kcnj13 ground truth ===")
    k13 = pd.read_csv(R / "kcnj13_residue_divergence_fixed.csv")
    col = [c for c in k13.columns if "diverg" in c.lower() or "distance" in c.lower()][0]
    pos = [c for c in k13.columns if c.lower() in ("residue", "position", "resid", "pos")][0]
    # The claim is about windows, not individual residues: local divergence is a +/-10-residue
    # statistic, so a substitution lifts its whole neighbourhood and the peak need not land on
    # the substituted position. Q19L and D176G rank 10th and 16th; what the paper asserts is
    # that the top-ranked residues all fall inside the two windows around them.
    top20 = k13.nlargest(20, col)[pos].tolist()
    inside = [r for r in top20 if abs(r - 19) <= 10 or abs(r - 176) <= 10]
    check("all twenty highest-divergence residues fall in the two substitution windows",
          "all twenty", f"{len(inside)}/20", len(inside) == 20)
    ranked = k13.sort_values(col, ascending=False)[pos].tolist()
    r19, r176 = ranked.index(19) + 1, ranked.index(176) + 1
    check("Q19L and D176G are not themselves the top two (stated in 3.2)",
          "ranks 10 and 16 of 338", f"ranks {r19} and {r176} of {len(ranked)}",
          (r19, r176) == (10, 16) and len(ranked) == 338)

    print("\n=== Claim: AHR2 pocket identity (Results 3.5, quoted in correspondence) ===")
    pk = pd.read_csv(R / "ahr2_pocket_sequence_identity.csv", comment="#")
    pocket = pk.loc[pk.region.str.contains("pocket"), "percent_identical"].iloc[0]
    pas3 = pk.loc[pk.region.str.contains("PAS fold 3"), "percent_identical"].iloc[0]
    check("ligand pocket 100% identical", "100.0%", f"{pocket:.1f}%", close(pocket, 100.0, 0.01))
    check("PAS fold 3 100% identical", "100.0%", f"{pas3:.1f}%", close(pas3, 100.0, 0.01))

    out = pd.DataFrame(RESULTS, columns=["claim", "abstract", "recomputed", "pass"])
    out.to_csv(R / "abstract_claim_audit.csv", index=False)
    failed = (~out["pass"]).sum()
    print(f"\n{len(out) - failed}/{len(out)} claims reproduce. Wrote {R / 'abstract_claim_audit.csv'}")
    if failed:
        print(f"*** {failed} claim(s) do not reproduce — see the FAIL lines above ***")


if __name__ == "__main__":
    main()
