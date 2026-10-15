Portable copy of the independent research-workspace review. Original report SHA256: 08356df53426f2b0eac686b58c282de9c23dee412d15134b4f07740a8fb1ea74. Local paths are portable provenance references. Final ZIP acceptance is recorded separately after archive verification.

# PUB.01 independent local-delivery review

Content review began: 2026-10-16T00:30:16+08:00. Content review completed: 2026-10-16T00:39:42+08:00. Reviewer: independent subagent stable_control_review.

**Content and software review: PASS after corrections.** This is acceptance of the standalone local PhysInfoBench v0.4.0 content and functioning scoring/replay interfaces. The final ZIP, embedded PUB.01 audit copy and accompanying checksum will receive a separate final micro-review once generated. No GitHub repository was created or pushed by this review. No blanket third-party dataset redistribution approval is granted.

## Scope and resolved findings

Reviewed the source modules, package/entry-point metadata, README and model/prediction contracts, task catalog, manifests and scientific exports, source/licensing statements, validation document, 17-item reviewer checklist, release tests, installed base entry points and reference replay. Three tasks have scored references; FS region stable controls are data-only; the new condition four-arm experiment and new encoders/structure predictors are deferred.

The initial review retained its failing logs and identified: (1) a release-wrapper test monkeypatch targeted the wrapper rather than the imported legacy function globals; (2) replay directly indexed incompatible legacy/new FS metric names; (3) temporary build/egg-info artifacts entered the root file manifest; (4) MODEL_INPUT retained a source-project test path; (5) the root verifier checked listed hashes but did not reject unregistered distributed files; and (6) checklist 4a overclaimed support for new biological sample batches when the interfaces are scoped to frozen benchmark identities. The implementer corrected all six. Replay now maps legacy/new metrics explicitly, verifies KNOT and FS confidence intervals and all 70 FS case-rank rows; the file-set guard rejects omitted files and symlinks; checklist 4a is No with evidence and a next action. No scientific label or reference prediction was changed.

## Functional and failure verification

Independent final test suite: **37/37 PASS**. Coverage includes malformed biological masks, bounded/duplicate/missing residue coordinates, missing and non-finite predictions, sample/sequence/endpoint identities, inner-selection contamination by outer-test groups, training-only preprocessing, failed-run records, exact FS endpoint multiplicity, nearest-reference JSON serialization, source-reference closure, changed file hashes, and unregistered delivery files. Independently added read-only temporary checks also reject prediction mask=0, wrong labels/groups, coordinate=0, infinite scores, a file absent from the manifest, and a symlink. Temporary synthetic fixtures were outside the delivery.

The final source-tree verifier passed. It verifies exact distributed file-set coverage (runtime output directories excluded), registered hash/size checks, the 13-task catalog with exactly three scored tasks, frozen group roles, pooled prediction/sample/arm/fold joins, 37 portable case-source references, the separate disorder coordinates/labels/folds, and primary/strict FS masks and counts. The reviewed snapshot has 226 registered files; the final embedded audit will increase the final inventory. Root manifests are regenerated after final metadata is added, and runtime directories must be absent from the actual ZIP.

Checked all 13 installed source modules against the release copies in the clean external environment: byte-identical. All eight base CLI --help entry points succeeded. The optional archived residue-training entry point requires the documented PyTorch extra and was not executed or retrained during this review. The installed environment is recorded separately from the original cluster training environment. The base package supplies scoring/replay without PyTorch.

## Independent complete reference replay

Ran the corrected source replay independently into SOURCE_45a68528c135_PUB.01_independent_replay_20261016_v2, outside the delivery. **PASS** at tolerance 1e-10. For all seven pooled arms, replay reproduces KNOT-PRESENCE AUROC/AP and group-bootstrap intervals, FS rank/Recall statistics and intervals, and every field of all 70 FS positive case-rank rows. For all seven residue arms, replay reproduces the complete archived disorder metrics and 2,000 binding-group bootstrap intervals/differences. This consumes all 67,291 pooled and 2,742,026 residue prediction rows.

Replay is numerical rescoring of saved predictions. It does not rerun ESM inference, refit readers, or reconstruct every inner fit/early-stopping trajectory. The prior MB.02.03/DIS.REF.01 scientific reviews and this local replay address their distinct scopes; neither file presence nor the software verifier alone proves biological validity or the claimed training history of an arbitrary external model.

## Scientific source fidelity and interpretation

Independently verified all 66 SOURCE_ARTIFACTS.json export hashes and their original research-source hashes. Byte-identical records match both hashes; transformed metadata are separately identified. Source/reference labels, coordinates, groups, folds, prediction tables, metrics and fitted-reader archives retain accepted scientific values. Local path substitutions and JSON serialization do not upgrade source evidence. The original research directories were not changed.

- KNOT-PRESENCE: 734 samples, 283 groups; nested reference AUROC 0.8687413872301332 and AP 0.6872110825517245.
- FS-L2: 8,879 operational ranking items, including 10 confirmed positive cases; nested mean held-case rank percentile 0.47659152548553035. Background biological labels remain unknown and masked.
- Separate disorder reference: 1,279 protein records, 1,275 sequences, 1,029 groups, 391,718 labeled protein-residue rows versus 391,189 unique sequence-coordinate vectors. Positives are direct disorder; negatives comprise 35 direct-order rows and 317,917 PDB-observed reference rows. The strict DISORDER-RES view is preserved separately.
- FS region controls: 10 pairs / 20 endpoints / 5,718 positions; 124 inherited positive rows, 848 primary structure-reference negative rows and 580 strict rows. Three dual-class cases contain 356 primary negative endpoint-residue rows. No region model score is claimed.

The documentation discloses four same-sequence record pairs and 12 differing cross-record disorder-reference coordinates, distinguishes the two mismatching old local caches from intact cluster execution inputs, preserves run-time pending-review strings with later acceptance sidecars, and avoids prospective/family-independent claims. The registered binding map does not certify complete biological independence. Condition case records remain distinct from the deferred new four-arm evaluation.

## Distribution, documentation and checklist

MIT is explicitly limited to the original software. KnotProt source-lineage/terms and record-specific literature redistribution remain unresolved, and the base data manifest retains DRAFT_NOT_CLEARED_FOR_PUBLIC_RELEASE. This is suitable for the user-authorized local folder/ZIP delivery; it does not approve making the complete snapshot public. Weights, embedding caches, raw article text and original PDB coordinate files are absent. Trained reader parameter NPZ files are identified as fitted readers, not encoder checkpoints. No symlinks were present in the reviewed distributed tree; the largest file is the compressed residue prediction table at 38,459,643 bytes, below 100 MiB.

All 18 Python source/test files have zero COMMENT tokens and zero module/class/function docstrings. Command help remains explicit, and scientific attribution/licenses remain in Markdown. Reviewed the final validation statements against actual reviewer/main-agent logs; the full-retraining limitation is clear. The 17-item checklist preserves the exact original PDF fingerprint, includes evidence for every answer and next actions for every No, and does not label external-database or split-leakage questions N/A. It records 4a as No: fixed benchmark new-model representations are supported, general new biological batches are not claimed.

## Preserved review evidence

- logs/reviews/PUB.01_20261016_tests.txt and PUB.01_20261016_replay.txt: initial failures retained.
- logs/reviews/PUB.01_20261016_tests_v3.txt: independent 37-test PASS.
- logs/reviews/PUB.01_20261016_verifier_v4.json: latest source-tree verification PASS.
- logs/reviews/PUB.01_20261016_replay_v2.txt: complete independent replay output; REPLAY.json is in the external temporary replay directory.
- logs/reviews/PUB.01_20261016_installed_cli.json: installed source equality and eight base entry-point checks.
- logs/reviews/PUB.01_20261016_content_checks.json: complete source/export fingerprints, 17 checklist answers, independent failure tests and reviewed distributed-file snapshot.

## Reviewed core hashes

| File | SHA256 |
| --- | --- |
| `README.md` | `aece731cf8baefee85db3dec33f862f0287ecacc085201555d13645224042a9d` |
| `VALIDATION.md` | `ea3f524eaeeba4e3b75fc7127d417262aaa6c76f4c95ab2a2348d6d1d7c00b58` |
| `DATA_REVIEWER_CHECKLIST.json` | `02e5583b256d213878008283863857ab8d4c263919df7b3f707f09dc03656bf0` |
| `pyproject.toml` | `f18a7cf2b6ac58b8c86a96613c98ab5772c17000b7a6a8a632e69824e5e3dc72` |
| `MODEL_INPUT.md` | `4e7631faf618e70d96e7ac7cfbb60c5eef37ebbac112141ba6da41106aa5a186` |
| `SCORING_PROTOCOL.json` | `071a7a78b7f16850ac8610865548daff6f9258e5d0b0ca12c4c9d8222df84514` |
| `scripts/evaluate_predictions.py` | `cdc15a306395a7947dc1300563757f78ad26d97e57d2655a7b68aea7cd7ed576` |
| `scripts/run_nested_readers.py` | `36f38f058dceccc7d474b44c16c17b1a26d8de66a66dbb61f81a408aa0e2c5b5` |
| `scripts/evaluate_residue_predictions.py` | `e046e04b08152a607b06c3d8e2d77f65727b9efab2f45f8ce4e1cba06c804dfc` |
| `scripts/residue_scoring_core.py` | `66e2e649e7347b8d93b37ffcaa98d1c280617994aa5d36e2e6368e6837e841b5` |
| `scripts/replay_reference.py` | `7c7ac697ad0a9509ed3fa9bd290be0752cbd55ff7147c6fc089b3f768702be86` |
| `scripts/verify_release.py` | `9f375bc369c55c23d858582f2fb90f69925d254942b856dc388395ccd67de6aa` |
| `scripts/verify_legacy_snapshot.py` | `f5d5c13017ade9c17e20ede8e5ed3d7f52e38a3c73e48fad20f68a0e70c35d8c` |
| `SOURCE_ARTIFACTS.json` | `bcf97a6983ad34a8820a99dbd6bf6f97192738866f50cd0a3ac108d73feebbc8` |
| `COMMENT_CLEANUP.json` | `f5c0b397cdb496b887b8f944137ff25e767c99771401dcb2d413166311046a9e` |

## Later task-status document check

A subsequent cross-document check identified stale candidate-state statements in results/task_dispositions.md. The implementer synchronized that page with the accepted three-task catalog and FS region reference annotations, and the independent reviewer checked the corrected counts and scope. VALIDATION.md records the correction. SOURCES_AND_DISTRIBUTION.md also now distinguishes scored KNOT-PRESENCE labels from unscored KNOT-TYPE annotations. These documentation corrections change no scientific data, prediction, metric or runtime source code. They receive a new manifest and final archive check; the initial captured review body and its source fingerprint above are retained.
