Historical validation of the accepted local delivery. Repository export verification is in audits/PUBLICATION_VERIFY.json; repository publication does not change the scientific validation scope below.

# Local delivery validation — PUB.01

Validation date: 2026-10-16, Asia/Shanghai. This record covers the standalone local v0.4.0 delivery. No remote repository has been created or pushed.

## Installation and functional checks

A new virtual environment was created outside the delivery directory. Installation built and installed the package and all base numerical dependencies. The actual versions are recorded in ENVIRONMENT_VALIDATED.json. This environment is distinct from the archived Python 3.11.16 / PyTorch 2.5.1+cu121 training environment.

The final functional suite has 37 tests, including malformed masks, wrong residue coordinates, incomplete/duplicate predictions, non-finite scores, exact sample/group identities, outer-test contamination of inner selection, training-only preprocessing, failed-run records and unregistered delivery files. Both the main agent and the independent reviewer execute the checks. See audits/PUB.01_review.md for the independent result.

The release verifier checks all registered file hashes and sizes, exact registered file coverage, three scored task dispositions, frozen group roles, 37 portable case-source references, disorder sequence/coordinate identities, fold-switch primary/strict labels and comment removal. Generated runtime directories such as virtual environments, build output and runs are excluded from the distributed file set check. They are absent from the delivered ZIP.

## Reference prediction replay

The complete replay passed in the clean installed environment. All seven pooled arms reproduce KNOT-PRESENCE AUROC/AP and group-bootstrap intervals, FS-L2 rank/Recall metrics and intervals, and all 70 FS positive case-rank records. Seven residue arms reproduce all per-protein aggregate metrics and the 2,000 group-bootstrap intervals/differences in the archived disorder reference metrics. The comparison tolerance is 1e-10. A separate independent replay also passed.

The replay covers 67,291 pooled prediction rows and 2,742,026 residue prediction rows. Its scope and timestamps are in audits/REFERENCE_REPLAY.json. This is complete rescoring of archived predictions; encoder inference, all inner training trajectories and reader refitting were not repeated.

## Source fidelity and cleanup

Original accepted scientific artifacts are preserved in the research workspace. The delivery retains labels and predictions without changing scientific values. Portable provenance exports are distinguished from byte-identical source copies in SOURCE_ARTIFACTS.json. The independent reviewer checked all 66 listed export fingerprints.

All 18 Python source and test files contain zero COMMENT tokens and zero module/class/function docstrings. COMMENT_CLEANUP.json records removal of 11 comments and 17 docstrings, with syntax-tree equivalence after removal. Existing command help text is stored explicitly. Markdown documentation, source attribution and license text are retained.

The package contains no checkpoint weights, feature cache archives, article full text or raw PDB coordinates. Fitted reader parameters are included and must not be confused with encoder checkpoints. Final distribution files are inventoried by manifest.json; the ZIP has a separate SHA256 checksum.

## Scope of acceptance

The delivered references are retrospective grouped-development results. The accepted FS region annotation remains data-only. The new condition four-arm experiment is deferred. Source-specific full-data public redistribution review remains open; software MIT licensing and local delivery acceptance do not close it. DATA_REVIEWER_CHECKLIST.json retains all 17 answers, the original PDF fingerprint, and evidence/next actions for each No.

## Task-status document synchronization

A final cross-document check found that results/task_dispositions.md still described the old candidate states for the disorder reference task and FS region negatives. The page is synchronized with the accepted task catalog and README in this delivery. Labels, predictions, scoring code and reference metrics are unchanged. The final archive is rebuilt and independently checked after this documentation correction.
