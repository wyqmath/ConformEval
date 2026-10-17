# Changelog

## Unreleased

- Add a runnable ESM-2 8M notebook with a Colab entry, full knot and synthetic modes, and independent rescoring checks.
- Recheck ESM-2 8M/650M and ProtT5 on all three synthetic task pipelines; retain untested full-task states explicitly.

- List compatible dependency ranges in the installation guide and keep exact historical requirements under `docs/reproduction/`.

- Remove snapshot-version suffixes from data, configuration, and result directory names.
- Update commands, scoring paths, tests, result indexes, and documentation to the new layout.
- Record original snapshot versions and source locations in [bundle metadata](../data/provenance/bundle_versions.json). Scientific versions, labels, folds, predictions, and saved execution records are preserved.
- Consolidate redundant ESM-2 8M arm prediction slices into the existing per-task `predictions.tsv.gz` tables.
- Add a single [data and table guide](data_tables.md) with public TSV/TSV.GZ purposes, fields, and reading examples.

## 0.5.1 — 2026-10-17

- Read UTF-8 inputs explicitly, accept an optional BOM, and write new JSON/FASTA text with LF line endings.
- Preserve frozen data, configuration and result bytes during Git checkout.
- Add text portability tests and a Windows timezone-data dependency; record actual operating-system installation/scoring coverage.

## 0.5.0 — 2026-10-16

- Add real ESM-2 and ProtT5 adapters with full coordinate projection and long-sequence windows.
- Configure layers, pooling, predictor candidates and nested selection; use one grouped trainer for sequence and residue tasks.
- Add per-task inputs/scoring, a complete three-task mini fixture, external prediction examples and optional inference dependencies.
- Add measured model-validation resources and two complete user tutorials.


## Repository organization — 2026-10-16

- Move the importable package to `src/confomeval` and use `confomeval-*` commands.
- Move usage documents to `docs`, protocols to `configs` and the task catalog to `data`.
- Add a minimal archive example and local scientific integrity tests.
- Remove workflows, delivery-only validators, root file inventories and internal publication records.
- Preserve the scientific labels, splits, predictions and fitted parameters.

## 0.4.0 — 2026-10-16 — initial scientific snapshot

- Preserve the v0.3.1 base scientific snapshot and two audited pooled ESM-2 reference tasks.
- Add the accepted disorder structure-reference data and 2,742,026 residue predictions across seven arms, with a portable scorer and fixed ESM-2 reader refit command.
- Add accepted paired-state fold-switch region controls as a separate annotation bundle.
- Add a reference prediction replay command, standalone requirements and source/export hash records.
- Refresh README, task catalog, data card, license scope, citation and contribution guidance.
- Keep the condition four-arm experiment deferred and unresolved source distribution entries explicit.

- Add the maintainer-authorized ConformEval public repository export, repository URL and citation metadata. At initial publication, scientific data and runtime source code matched the accepted local snapshot; source-specific data terms remain unresolved.

## Earlier local candidates

Version 0.3.1 replaced case-module source paths with portable source references. Version 0.3.0 fixed the original local data and pooled reference result snapshot. These candidate versions were not public releases.
