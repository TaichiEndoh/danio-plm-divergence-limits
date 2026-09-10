# data/reference/ — baseline data ported from the previous paper

Committed (small) reference data from
[`zebrafish-aescallii-go-plm`](https://github.com/TaichiEndoh/zebrafish-aescallii-go-plm).
These are the **ESM-2 8M baseline** and inputs that this paper's ESM-3 / AlphaFold work builds on.

| File here | Source (previous repo) | What it is |
|---|---|---|
| `ortholog_pairs_rbh.csv` | `zebrafish_aescallii_GO_pipeline/ortholog_table_relio_aesculapii.csv` | 68,971 RBH ortholog pairs + bidirectional bitscores + `A_to_B_qcov` (⚠ ratio-like, >1 possible — recompute/relabel before citing as coverage) |
| `esm2_8m_pair_distances.csv` | `plm_go_analysis/analysis/ortholog_plm_distance.csv` | ESM-2 8M per-pair cosine distance — **baseline for the 8M-vs-1.4B comparison** |
| `ahr_local_pocket_esm2_8m.csv` | `plm_go_analysis/analysis/ahr_local_pocket_plm.csv` | AHR subtype local-pocket distances (Takeda residues), ESM-2 8M |
| `local_identity_vs_plm_esm2_8m.csv` | `plm_go_analysis/analysis/local_identity_vs_plm_combined.csv` | local identity vs local PLM distance (targets), ESM-2 8M |
| `sequences/ahr1a_full.fasta` | `plm_go_analysis/analysis/ahr1a_full.fasta` | AHR1a full-length target sequence |
| `sequences/target_domain_sequences.fasta` | `plm_go_analysis/analysis/target_domain_sequences.fasta` | target domain sequences |

## Not ported (still in the previous repo)
- Full-proteome FASTAs (were external inputs, not committed there either).
- `plm_go_analysis/results/plm_go_table.tsv.gz` (4.4 MB integrated table) — fetch from the previous
  repo if the full GO-joined table is needed.

## Derived: `esm2_8m_pair_distances_dedup.csv` (provenance)

`esm2_8m_pair_distances.csv` is kept **byte-for-byte as ported** from the previous repo
(do not edit it). It has two artifacts of that repo's batched, append-mode PLM run — both
are *bookkeeping* issues, not distance errors:

1. **One corrupted value**: `4.547834396362305e-05git` (row `NP_001001403,XP_056310087`,
   present since the port commit). The stray `git` token forces the whole distance column
   to `dtype=str`.
2. **Duplicate rows**: 68,971 rows but only **50,921 unique pairs** — 18,050 pairs appear
   twice (append-mode runs over overlapping row ranges). Duplicates agree on distance to
   float noise (max spread ≈ 2.4e-7; only 8 pairs differ at all).

`analysis/build_baseline_dedup.py` produces `esm2_8m_pair_distances_dedup.csv`: numeric
distances, one row per unique pair (50,921). The corrupted pair's value is recovered from
its clean duplicate row, so 0 rows are NaN. **Use the dedup file for per-pair analyses.**

### ⚠ Percentiles: published values are row-wise (with duplicates)

The previous paper's stated percentiles reproduce **exactly** only over the FULL 68,971-row
file (duplicates included), because the duplicated pairs are disproportionately low-divergence
and pull the tail down. Per-unique-pair percentiles are higher:

| statistic | published (row-wise, full file) | per-pair (dedup, 50,921) |
|---|---|---|
| median | 0.000487 ✓ exact | 0.000496 |
| 95th   | 0.008432 ✓ exact | 0.009815 |
| 99th   | 0.039693 ✓ exact | 0.050938 |

Do not expect the dedup file to reproduce the published row-wise percentiles.

### Pair-count reconciliation (published 51,086 vs repo unique 50,921 → 165 gap)

The published "PLM-joined" pair count is **51,086**; this ported file dedups to **50,921**,
i.e. it is missing distances for **165** pairs the published analysis had (and, more broadly,
for 18,050 of the 68,971 RBH ortholog pairs). Those missing pairs are **not** explained by any
principled filter: across the ortholog table, `RBH` is `True` and `rescued` is `False` for all
pairs, and the `A_to_B_qcov` distribution of the missing pairs matches the present ones (the
present group even has *more* qcov < 0.8). So the gap is a **partial/append-assembled export**,
not a filter-stage difference. Pinpointing the exact 165 pairs requires the previous repo's
`plm_go_analysis/results/plm_go_table.tsv.gz` (the 51,086-pair published table) — not yet
imported here.

## Note on the ortholog table
The previous repo's `A_to_B_qcov` column ranged 0.60–1.808. A query coverage cannot exceed 1.0, so
this column is a ratio-like metric (mislabeled). The published coverage stats (mean 0.9901, median
1.0, min 0.8, max 1.0) came from a filter-consistent recomputation. Recompute/relabel before reuse.
