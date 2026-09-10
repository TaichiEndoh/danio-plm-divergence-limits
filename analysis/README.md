# analysis/

New analysis scripts for this paper (ESM-3 scaling + AlphaFold structural validation).

- `build_baseline_dedup.py` — derive `esm2_8m_pair_distances_dedup.csv` (one row per unique
  pair) from the pristine ported baseline; see `data/reference/README.md` for provenance.
- `select_high_divergence_subset.py` — build the ESM-3 re-embedding subset as a labeled union
  of four selection modes: decile-**stratified** sample (unbiased backbone), **top**-N tail,
  **random** control, and curated **target** proteins (AHR subtypes, Kcnj13 + comparators).
  Each pair carries which mode(s) selected it, so the analysis can separate a real cross-model
  signal from regression-to-the-mean in the tail. CPU-only, seeded, deterministic.
- `esm2_variant_score.py` — ESM-2 variant-effect scoring (masked-marginal Δlog-likelihood):
  score specific substitutions, run a full single-substitution scan (L×19), and report a known
  variant's percentile rank among all substitutions (the "rediscover known functional variants"
  test). Core primitive for the expanded design — see `docs/PAPER_DESIGN_ja.md`. CPU-verified.
- `esm2_ladder_embed.py` — the paper's MAIN scaling axis: per-pair cosine distance across
  ESM-2 8M / 35M / 150M / 650M, one size per invocation, GPU/CPU/MPS, resumable, cached.
  Pooling is identical to the previous repo; verified end-to-end — the 8M distance for the
  AHR2 pair (NP_571339 ↔ XP_056303610) reproduces the ported baseline to float precision
  (0.0034116 vs 0.0034114). `--demo` runs the flow without weights.
- `esm3_embed.py` — embed sequences with ESM-3 (1.4B) and compute pairwise cosine distance.
  ESM-3 is the architecture-robustness confirmation (a different model family), not the
  primary evidence. Skeleton; `--demo` runs the flow with deterministic fake embeddings.
- `esmfold_predict.py` — predict a target's 3D structure with ESMFold (the chosen structure
  route: full-length AHR2 has no AlphaFold DB model, and ESMFold folds from a single sequence,
  staying in the ESM family). `--backend api` uses the ESM Atlas API (no GPU, seqs up to ~400 aa,
  good for Kcnj13); `--backend local` uses fair-esm ESMFold on a GPU for long proteins like AHR2.
- `render_structure.py` — render a B-factor-colored PDB to a PNG (matplotlib 3D Cα trace),
  so the structural figure needs no PyMOL/ChimeraX. Pairs with `alphafold_local_mapping.py`.
- `alphafold_local_mapping.py` — map local divergence onto a predicted structure (ESMFold `--pdb`
  or AlphaFold DB `--uniprot`); writes values into the B-factor column for PyMOL/ChimeraX. Use
  `--scale 1000` so small cosine distances give visible contrast. `--demo` uses synthetic values.
  Verified end-to-end: ESMFold-predicted AHR2 pocket fragment colored by real ESM-2 divergence.

Port the local-divergence / domain scripts from the previous repo into this folder as needed
(see `docs/REUSE_FROM_PREVIOUS.md`).
