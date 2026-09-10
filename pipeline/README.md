# pipeline/

Ortholog inference and data-preparation pipeline, ported/adapted from
`zebrafish-aescallii-go-plm` (`zebrafish_aescallii_GO_pipeline/`).

Bring over only what this paper needs:
- ortholog pair table (MMseqs2 RBH; qcov ≥ 0.8, length ratio 0.8–1.25)
- protein FASTA extraction for D. rerio / D. aesculapii
- (optional) GO annotation transfer, if GO analysis is retained

See `docs/REUSE_FROM_PREVIOUS.md`.
