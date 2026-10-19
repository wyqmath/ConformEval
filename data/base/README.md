# ConformEval data snapshot v0.3.1

This is a portable export of the frozen 2026-10-14 development dataset plus the independently reviewed FD.03 representation-coverage overlay. Labels, masks, sample identities, and folds are copied without modification. Case-module rows point to `case_sources.tsv` through `source_ref`; that table records the upstream artifact name and available DOI, PMCID, or source notes without exposing project-local paths. `representation_coverage.tsv.gz` contains no local cache paths.

Distribution status: **draft, not cleared for public release**. Consult `sources.tsv` and the benchmark-level distribution notes before copying these files to a public host. No model weights or embedding caches are included.
