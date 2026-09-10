# Limits of protein language model ortholog divergence

Analysis code and computed outputs for:

> **Limits of protein language model ortholog divergence: model size, preprocessing and
> structural noise in two *Danio* species.**
> Taichi Endoh, Gerry Amor Camer, Kotetsu Kayama, Daiji Endoh, Hiroki Teraoka.
> Submitted to *Proteomes* (MDPI).

Follow-up to Endoh et al. 2026, *Proteomes* **14**(3):36.

This repository holds everything needed to check the paper's numbers: the scripts that produced
them, the per-analysis outputs they wrote, and the scripts that regenerate every figure. **Every
value reported in the paper is recomputed from the files in `reports/`.** No new experimental data
were generated.

## What the study asks

Protein language model (PLM) embeddings are increasingly used to quantify divergence between
orthologous proteins, but the conditions under which that quantity is trustworthy have not been
established. This is a methods-limits study: over 3,428 stratified *Danio rerio* / *Danio
aesculapii* ortholog pairs it compares an ESM-2 ladder (8M–650M parameters) against ESM-3 (1.4B)
under identical sequences, pooling and length handling, benchmarks the ranking against sequence
identity rather than against chance, and maps per-residue divergence onto ESMFold structures and
domain annotation.

Rank, not absolute distance, is the transferable quantity (Spearman ρ = 0.881, 95% CI
0.870–0.891). Three limits follow — one on model size, one on preprocessing, one on structural
grounding — and four preprocessing hazards are reported because each is easy to hit and hard to
notice.

## Layout

| Path | Contents |
|---|---|
| `analysis/` | Every analysis, figure and audit script. One concern per file. |
| `reports/` | The committed outputs. Every number in the paper comes from here. |
| `reports/figures/final/` | The composed figures at 300 dpi, with a README mapping them to source panels. |
| `reports/external/` | Third-party reference material, with provenance and licence in its README. |
| `data/reference/` | Input sequences (FASTA) and the ortholog pair table. |
| `pipeline/` | End-to-end driver. |
| `docker/`, `.devcontainer/` | GPU environment used for the embedding runs. |

## Reproducing

```bash
pip install -r requirements.txt

# recompute the sequence-level facts behind the AHR2 result
python analysis/pocket_sequence_identity.py

# recheck every quantitative claim the Abstract makes, from the committed CSVs
python analysis/audit_abstract_claims.py
```

`audit_abstract_claims.py` is the quickest way to confirm the paper against this repository: it
restates each claim as printed, recomputes it from `reports/`, and fails loudly on any that does
not reproduce.

The embedding runs themselves need a GPU. ESM-2 uses `fair-esm` 2.0.0 and ESM-3 uses the
`esm` package for `esm3-sm-open-v1`; the two import under the same name, so they are run from
separate environments. ESM-3 access requires accepting the model licence on Hugging Face.

## Data sources

All protein sequences are publicly available from NCBI RefSeq; accession numbers are given in the
paper and in `data/reference/`. Structures are ESMFold (`facebook/esmfold_v1`) predictions
computed here, not deposited experimental structures. Domain annotation is from InterPro and
UniProt.

## Licence

Code is released under the MIT Licence (`LICENSE`). The computed outputs and figures in
`reports/` are released under CC BY 4.0. `reports/external/` contains third-party material that
remains under its own licence, recorded in `reports/external/README.md`.

## Citation

Please cite the paper. If you use this repository directly, cite the archived release rather than
the branch.
