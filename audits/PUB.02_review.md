Portable copy of the independent pre-publication content review. Source report SHA256: e25435823c0ba8eb59e6005d5f961d3fa08b4cf8cb28853018375eeb2bcc635c. This records readiness; publication is verified by the later remote commit/tree readback.

# PUB.02 ConformEval repository-export independent review

Review began: 2026-10-16T01:14:34+08:00. Content review completed: 2026-10-16T01:19:29+08:00. Reviewer: independent subagent stable_control_review.

**Verdict: PASS for the prepared repository-export content and source fidelity.** The final manifest, filled publication-verification sidecar, embedded PUB.02 review and actual Git checkout/tree must receive the planned final limited check before push. At this point the new version has not been published, and this review does not assert that it has. No new numerical model fits are needed because the scientific assets and runtime sources are unchanged.

## Authorization, repository access and rollback evidence

The delegated task records the user explicitly requesting replacement of https://github.com/wyqmath/ConformEval with the accepted benchmark, with deletion of the previous current-tree contents. The stated method is an ordinary main-branch commit retaining prior history, not repository deletion, visibility changes or force-pushing. This is consistent with the requested action; no additional user confirmation is inferred by this review.

Independently queried GitHub read-only: repository wyqmath/ConformEval is public, default branch main, authenticated permission push=true, main protected=false, current main SHA 0292bdbf7a2e1d7f3dfa5b3b121b966b462f3be6. This confirms the stated operational target and write capability; no target-permission discrepancy was found. Immediately before push the implementer must still use normal non-force semantics, so a concurrent main update cannot silently be overwritten.

Independently verified deliverables/ConformEval-before-v0.4.0-20261016.bundle with git bundle verify: complete history, 6 refs, including the exact current main SHA, the older codex branch and snapshot tag. The bundle is local backup material and is not included in the public export. The ordinary replacement commit must retain the old main as its parent/ancestor.

## Byte fidelity and declared publication changes

The accepted local source remains immutable: all 228 registered local files still match its source manifest, whose SHA256 is 6a7b657489f6c0c68a4292ef238d96625a7e364ee7ee8928b2ed68dd4aed19ca; the original local ZIP retains the reviewed 295d7e6fb3fee6fa41614ba9d2d9838650a9d35602548771a5193e6f166bdc4b identity. SOURCE_LOCAL_DELIVERY_MANIFEST.json is byte-identical to that original manifest.

Independently compared all shared paths: all **195 data/result/runtime-source/test files** are byte-identical to the accepted local delivery. This includes labels, masks, exact sequences, groups, inner/outer folds, predictions, metrics and fitted-reader parameter archives. No biological sample, protocol constant, model result or qualification was silently changed.

PUBLICATION.json correctly records the original manifest, target repository/branch, source research commit, maintainer instruction, unresolved source status, and all 13 changed publication-document/metadata paths. The root manifest is separately changed to identify a repository snapshot. The 13 declared changes are .gitignore, CHANGELOG, CITATION, CONTRIBUTING, DATA_CARD, DATA_REVIEWER_CHECKLIST, LICENSE_NOTICE, README, RELEASE_STATUS, SOURCES_AND_DISTRIBUTION, UPLOAD_GUIDE, VALIDATION and pyproject.toml. Their source/export hashes match. New files are PUBLICATION.json, SOURCE_LOCAL_DELIVERY_MANIFEST.json and the two publication audit/verification sidecars. No source-delivery path was removed from the export.

The package retains its PhysInfoBench command/task identifiers; the README explains that the repository is named ConformEval. Package metadata only adds the repository URL, and citation metadata uses the existing real repository URL without inventing a DOI. The .gitignore changes exclude generated runs, wheelhouse and virtual environments; they do not hide the required scientific files.

## Scientific and documentation consistency

Checked current README, catalog, dispositions, source/license notes and publication status against the accepted snapshot: exactly three scored references remain (KNOT-PRESENCE 734 samples; FS-L2 8,879 ranking units; disorder reference 1,279 records/1,029 groups/391,718 labeled rows/2,742,026 prediction rows). Strict DISORDER-RES remains distinct, FS region controls remain data-only, and the new condition four-arm experiment remains deferred. The corrected current task-dispositions page is preserved. KNOT-TYPE annotations are not described as a completed KNOT-TYPE score.

The source-review text still records missing legacy KnotProt lineage/terms and record-specific literature questions. MIT is limited to the original software; publication is not portrayed as a blanket dataset license. The original data manifests retain DRAFT_NOT_CLEARED_FOR_PUBLIC_RELEASE. The root repository-snapshot distribution status explicitly says source terms unresolved. This accurately separates the user's publication instruction from third-party source permissions; this review is not a source-rights clearance decision.

Local-only validation and prior scientific review records are retained as historical records with a current wrapper explaining their provenance. They are not used to claim that a remote push has already happened or that those old local acceptance records cleared source terms. DATA_CARD now distinguishes excluded encoder checkpoint weights from included fitted-reader parameters. No raw papers, original PDB coordinates, checkpoint weights, embedding caches, credentials or machine-local provenance paths were introduced.

The 17-item checklist preserves the original PDF fingerprint and all prior answers except 6a. The staged 6a=Yes is a code-publication statement for the intended published snapshot, not a source-permission statement. **For this version it must be substantiated after push by reading the target commit/tree back from GitHub; a prepared local export alone is not proof that this code version is publicly accessible.** No pre-push completion report should claim that the new snapshot has already been pushed.

## Functional verification

Independently ran all 37 focused tests: PASS. Independently ran the installed source-tree verifier against the export: PASS, checking registered hashes/file membership, data masks/coordinates, complete frozen grouped roles, source links, scientific counts, 18 Python files without comments/docstrings and absence of machine-local paths. Its initial snapshot had 231 registered files; the final review/verification sidecars are still to be locked. Checked 24 user-facing internal documentation links: all resolve.

The full three-task reference replay already passed for the accepted local artifact. Since the exported runtime code, data, predictions and metric references are identical, no new full replay or model refit was performed or claimed here. Publication-document changes cannot supply independent biological evidence or complete family isolation.

## Final-tree and post-push scope

The remaining pre-push check is mechanical: replace audit placeholders with the scoped review/actual verification result, regenerate the manifest, copy exactly into the target checkout, and compare the tracked Git blobs/path set to the reviewed export. The latter ensures no ignored files, Git-LFS filters, temporary files or old-tree remnants change the intended publication. After a normal push, read back GitHub main SHA and the published file/tree set; record the observed commit URL and result in the project task log before reporting publication complete. No extra user permission request is required for that already authorized next step.

## Preserved evidence

- logs/reviews/PUB.02_20261016_repository_access.json: independent target/permission readback.
- logs/reviews/PUB.02_20261016_bundle_verify.txt: complete-history backup verification.
- logs/reviews/PUB.02_20261016_content_checks.json: exact export changes, all prepared-file fingerprints, byte-identical scientific/source count and local-manifest fidelity.
- logs/reviews/PUB.02_20261016_tests.txt: 37-test PASS.
- logs/reviews/PUB.02_20261016_verifier_initial.json: independent initial export verification.

## Current content hashes

| File | SHA256 |
| --- | --- |
| `PUBLICATION.json` | `deaf282060b35be2be01507ca45be7383ca332966d24b5168d029d2172f059a8` |
| `SOURCE_LOCAL_DELIVERY_MANIFEST.json` | `6a7b657489f6c0c68a4292ef238d96625a7e364ee7ee8928b2ed68dd4aed19ca` |
| `README.md` | `44e761e244522463b5213f29e4771296c3ef7a0f121b29931931cbe6295bee0c` |
| `CHANGELOG.md` | `a0675336600807a63b3b117ec262eb6447d061df820856d3ee7d30a5d2b0ff08` |
| `CITATION.cff` | `87a8bfd844668ce151fbc73935b5b65311b73b0c01e8d502c6714291a2f6d2aa` |
| `pyproject.toml` | `a4f24d2beac072affabe7f8fb422f67142c651fb6735b8a6dc1f31d5b20268e3` |
| `DATA_CARD.md` | `ca3109f96900eee8e976a518fb069948e5876ef8fb0b846d5ac7e59930049eae` |
| `DATA_REVIEWER_CHECKLIST.json` | `991cfc8ac74b4f689fd59d5e7d181a25bbe6a03c65e59123f9d1384471a75006` |
| `SOURCES_AND_DISTRIBUTION.md` | `44811e9f682f3c4852b1f32411c16c17d7d52edcb5ed66483512dbfd64ce8377` |
| `LICENSE_NOTICE.md` | `1aab5833739619be94b04fe6081307552017aca7e779526660c4a7238b9b5708` |
| `RELEASE_STATUS.md` | `cd5fe556cd4dbf102b52798c2ee8fe6c0e869ef74acdb179b07a450798a0af60` |
| `UPLOAD_GUIDE.md` | `fabe25d368e020ad5a122ca619d2fff3732df5292a8c5a485c19e3f8c2ebfe56` |
| `VALIDATION.md` | `f6b18f65aef7bbc2be04fe3dd652932be597e2645ef6782d3de0f9ebffed39b4` |
| `scripts/verify_release.py` | `9f375bc369c55c23d858582f2fb90f69925d254942b856dc388395ccd67de6aa` |
| `scripts/evaluate_residue_predictions.py` | `e046e04b08152a607b06c3d8e2d77f65727b9efab2f45f8ce4e1cba06c804dfc` |
| `task_catalog.tsv` | `f45870b7d8aa47c911df2a45c16861e86a53adabadbb3ab7c813719c5529d170` |
