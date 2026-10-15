# Repository snapshot

This repository hosts the maintainer-authorized v0.4.0 PhysInfoBench snapshot at https://github.com/wyqmath/ConformEval. The repository name is ConformEval; package commands and scientific task identifiers retain their PhysInfoBench names.

The snapshot replaces the previous main-branch file tree in a normal Git commit. Previous commits remain in Git history. The original accepted local directory and ZIP are preserved separately; SOURCE_LOCAL_DELIVERY_MANIFEST.json fingerprints that source delivery. PUBLICATION.json records the prepared export, original source manifest and changed paths.

## Scope

Three tasks have accepted ESM-2 reference evaluations. FS region stable controls are data-only; new condition four-arm evaluation is deferred. No encoder weights, cached representations, raw papers or PDB coordinates are included. Labels, groups, folds, predictions and reader parameter files retain accepted bytes.

## Licensing

MIT applies to original software. Source-specific terms for legacy KnotProt-derived annotations and some literature-derived records remain unresolved. The maintainer's publication instruction is recorded separately from source permissions; no blanket data license or source-clearance claim is made. Existing scientific data manifests retain their original review status. See SOURCES_AND_DISTRIBUTION.md.

## Verification

Use the install, tests and complete reference replay commands in README.md. The root manifest fingerprints this export. Earlier local audit records remain historical evidence; the repository export's verification is recorded in audits/PUBLICATION_VERIFY.json. GitHub Actions runs tests, file verification, reference rescoring and wheel build for this repository.
