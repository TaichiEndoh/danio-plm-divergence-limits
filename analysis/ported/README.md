# analysis/ported/ — scripts brought from the previous paper (ESM-2 era)

Verbatim copies from
[`zebrafish-aescallii-go-plm`](https://github.com/TaichiEndoh/zebrafish-aescallii-go-plm),
kept here as a starting point to **adapt for ESM-3**. Do not assume they run unchanged (paths,
model name, and inputs differ in this repo).

| File | Original | Purpose | Adapt for ESM-3? |
|---|---|---|---|
| `plm_go_table.py` | `plm_go_analysis/plm_go_table.py` | ESM-2 embedding + cosine distance over the ortholog table | **Yes** — swap ESM-2 8M → ESM-3; keep pooling (final layer, exclude special tokens, mean) |
| `ahr_pocket_local_plm.py` | `plm_go_analysis/analysis/ahr_pocket_local_plm.py` | AHR pocket local-window distances (Takeda residues) | Yes — reuse residue defs; re-embed with ESM-3 |
| `kcnj13_local_plm.py` | `plm_go_analysis/analysis/kcnj13_local_plm.py` | Kcnj13 local-window distances | Yes |
| `domain_plm_analysis.py` | `plm_go_analysis/analysis/domain_plm_analysis.py` | domain-level distances (Pfam/InterPro) | Yes |
| `go_graph_plm_compare.py` | `plm_go_analysis/analysis/go_graph_plm_compare.py` | GO-graph distance vs embedding distance | Optional (only if GO analysis retained) |

New ESM-3 / AlphaFold scripts live one level up in `analysis/` (`esm3_embed.py`,
`alphafold_local_mapping.py`).
