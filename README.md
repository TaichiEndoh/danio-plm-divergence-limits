# Cross-model concordance and input-quality limits of protein language model divergence in two Danio species

Submission data and code snapshot `v0.1.0-review.20260928` for the manuscript by Taichi Endoh, Gerry Amor Camer, Kotetsu Kayama, Daiji Endoh and Hiroki Teraoka.

This version corresponds to the manuscript forwarded on 26 September 2026 and the submission figures supplied on 27 September. It replaces the outdated September 10 repository snapshot for purposes of checking this manuscript. It does not assert that journal submission or all-author approval is complete. The repository is currently **private**. A repository URL alone does not give an editor access.

## What can be checked

The release contains retained protein sequences, pair-level outputs from four ESM-2 checkpoints and ESM-3, pooled embedding vectors, residue-level outputs, predicted structures, six manuscript tables, ten final figures, and analysis code. It supports saved-output numerical verification and figure reconstruction without downloading model weights or running inference.

The primary set has 2,662 untruncated pairs from a selected union of 3,428 pairs. The ESM-2 8M–ESM-3 rank correlation is 0.909 (connected-component bootstrap 95% CI 0.897–0.921). Shared proteins and the selected sample limit generalization. Model concordance does not demonstrate functional-prediction accuracy.

Historical inference is not fully reproducible from the available provenance. Historical model-weight hashes, complete annotation snapshots and raw L×19 variant-score distributions are unavailable for some analyses. The retained vectors allow distances to be checked but cannot establish the exact historical model/input provenance. These limitations are not repaired by creating a release.

## Quick verification

Use Python 3.9 and the saved-output dependency versions in `requirements.txt` (the reference environment for the manuscript figures). In a virtual environment:

```sh
python -m pip install -r requirements.txt
python analysis/verify_release.py
```

The verifier checks the file manifest, input/output fingerprints, all reported whole-protein correlations and top-k overlaps, pooled-vector distances, and the current sequence/table/figure counts. It does not claim biological validation or recreate missing variant scores.

To rebuild descriptive statistics and figures, use a **working copy** so the frozen files remain intact:

```sh
python analysis/revision_statistics.py
python analysis/revision_figures.py
```

The bootstrap uses 2,000 connected-component replicates, seed 0. Rendering can differ across software versions even if numerical results agree. The release records the executed environment and verification in `metadata/verification.json` and `metadata/rebuild_check.json`. A verification failure after intentionally regenerating files can reflect changed file hashes; compare numerical results before interpreting it as a scientific discrepancy.

## Files

- `reports/paper_revision/pair_table.csv`, `ranking.json`: selected pairs and descriptive summaries.
- `reports/figures/final/`: the current ten figures; full legends in `FIGURE_LEGENDS.txt`.
- `tables/Table1.csv` through `Table6.csv`: exact cell text extracted from the supplied manuscript Word file.
- `data/reference/`: retained inputs, including subset FASTA sequences.
- `data/embeddings/`: five compressed NPZ archives of saved per-protein pooled vectors. Keys are normalized protein accessions, values are float32 arrays. These are not model weights.
- `reports/*_fixed.csv`: corrected ESM-2 residue profiles used for the current figures. Other historical profiles are retained for comparison and are **not interchangeable** with the corrected profiles.
- `reports/*_esmfold.pdb`: predicted structures, not experimentally determined structures.
- `analysis/`: saved-output regeneration and selected historical inference scripts. `--demo` options, where present in historical scripts, generate synthetic flow checks and must never be used as scientific data.
- `metadata/`: input fingerprints, manuscript hash, verification and source-vector hashes.
- `LIMITATIONS.md`: unresolved scientific and provenance issues.

## Inference scope

New inference is optional and was not run to prepare this release. The historical ESM-2 scripts use `fair-esm==2.0.0`; ESM-3 scripts use EvolutionaryScale `esm` in a separate environment because both import as `esm`. The historical dependency notes are in `requirements-inference-historical.txt`; they are not a tested new-inference lockfile. Model access and licenses must be obtained from the providers. Cache keys in historical scripts do not fully encode sequence/model provenance; use a new empty cache for any new inference and record its inputs separately. Do not present such a new run as a recovery of undocumented historical provenance.

## Data sources and permissions

Retained protein sequences originate from NCBI RefSeq and UniProt. Accession identifiers and the exact retained sequence bytes are provided. Some historical accessions were normalized by dropping version suffixes; do not infer a missing version or retrieval date. Annotation intervals derive from InterPro/UniProt and remain subject to the snapshot limitations described in `LIMITATIONS.md`.

The existing repository's MIT license is retained for code. Its existing CC BY 4.0 licensing statement is retained for author-generated computed outputs and figures (https://creativecommons.org/licenses/by/4.0/). Third-party sequences and resources retain their respective terms. Model weights and third-party journal PDFs are not redistributed. No new rights over third-party material are claimed.

## Access during review

The private GitHub repository requires explicit access. For confidential peer review, provide the matching ZIP through the journal submission system using the editorial office's designated confidential file channel; ordinary supplementary files may later be published, so the submitting author should label the review-stage status explicitly. Access is **not** automatically granted by this release, and no editor or reviewer has been invited. At public release, update the manuscript's availability wording to match the actual visibility. No DOI has been minted; the Git commit hash, release tag and SHA-256 checksums identify this snapshot.
