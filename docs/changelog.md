# Changelog

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
