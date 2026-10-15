# FS.CTRL.01 independent review

Reviewed at: 2026-10-15 22:38 Asia/Shanghai. Reviewer: independent subagent stable_control_review.

**Verdict: PASS, scoped to the v04 structure-derived annotation deliverable.** No blocking findings remain. This accepts the local stable-reference labels and reproducible annotation implementation. It does not approve new scoring folds, model evaluation, full benchmark completion, prospective independence, or public redistribution.

## Scope and correction

Reviewed the builder, registered configuration, targeted tests, v04 artifacts, task report, and the FS.CTRL.01 authorization/TODO entry. Earlier v03 omitted explicit derived evidence IDs for some control-only negative rows where the parent case evidence was NA. The implementer corrected this in v04; the earlier directory remains historical. v04 has 76 derived segment-evidence records covering both primary and strict variants. Every negative row resolves to its own structural evidence, with endpoint identities, experimental structure hashes, protocol hash, and coordinate ranges. Parent literature evidence may separately remain NA.

## Independent verification

- All 149 registered source/input hashes, current builder hash, and all registered output hashes match. Sequence hashes, lengths, complete 1-based model positions, original/extracted structure provenance, and development group identities were checked. Each paired endpoint remains in its same binding group. The map does not establish complete protein-family independence.
- Parsed the original experimental mmCIF sources independently, rather than reusing the production chain-extraction coordinates. Verified unique corresponding C-alpha identities, residue types, author/label/sequence coordinates, finite coordinates, occupancy, single-model policy, identical structured DSSP centers, 9-residue physical continuity, author/positive exclusion buffers, and expression-tag exclusions.
- Independently recomputed all 424 primary paired negative positions using a quaternion eigenvalue RMSD solution, separate from the production SVD Kabsch method; independently enumerated raw-coordinate neighborhoods and correspondence coverage. All pass the registered primary gates; all 290 strict paired positions pass the strict gates. Maximum difference from stored RMSD/distance/coverage metrics is 5.624806176385277e-13. A 1e-12 comparison tolerance handles floating-point representation of exact boundary fractions.
- Verified all 5,718 full-length rows in both variants: unique identities, valid mask/label combinations, all 124 original positive keys unchanged, no positive-negative overlap, unknown positions remaining masked, strict negatives a subset of primary negatives, and no gap bridging. Each negative segment has at least 5 consecutive centers in both states.
- Verified all 76 derived evidence records against original endpoint provenance, source/protocol hashes, model/label/auth boundaries, and all referencing primary/strict negative rows.
- Ran all 7 targeted unit tests successfully. Added a separate real-CLI contamination experiment: after changing one endpoint to a different synthetic binding group and supplying a matching synthetic test manifest, the builder rejects the pair before output creation. Original inputs were not modified.
- Rebuilt v04 into SOURCE_e37fff242d9f_fs_region_stable_controls_independent_review_20261015_v04. All 11 non-run artifacts match byte for byte; run.json input/code hashes, dependency versions, summaries, warnings, and output hashes match. Timestamps and output-directory arguments appropriately differ.

## Counts and scientific interpretation

The product covers 10 experimental structure pairs, 20 endpoints, and 5,718 model positions. Primary labels contain 124 inherited positives, 848 stable-reference negative endpoint rows (424 paired positions; 47 paired segments / 94 endpoint intervals), and 4,746 unknown rows. Strict sensitivity contains 580 negative endpoint rows. Three cases in three binding groups have both qualified inherited positives and stable references: 124 positive and 356 primary-negative endpoint rows in total. Endpoint rows and overlapping windows are not independent observations.

The negative meaning is local structural preservation between the two observed experimental states under fixed geometric rules. These are geometry-derived references, not author-asserted negative annotations or proof that a region never switches. Structured-center selection favors well-resolved helices and sheets; flexible regions and boundaries are underrepresented. The thresholds are operational project choices and were not selected using model scores. Frozen historical labels and scoring protocols remain unchanged. Current task/report wording accurately preserves these limits.

The decisions log still mentioned v03 at review time; the implementer was notified to append the v04 acceptance/version record when formally closing the task. This administrative update must not alter scientific artifacts.

## Reproduction evidence

- `logs/reviews/FS.CTRL.01_stable_controls_20261015_verify.py` — independent original-coordinate and evidence verifier. Run with `.venv_local/bin/python logs/reviews/FS.CTRL.01_stable_controls_20261015_verify.py v04`.
- `logs/reviews/FS.CTRL.01_stable_controls_20261015_geometry_v04.json` — complete geometric and source verification results.
- `logs/reviews/FS.CTRL.01_stable_controls_20261015_rebuild_v04.json` — deterministic rebuild and binding-contamination rejection results.

## Exact reviewed hashes

| File | SHA256 |
| --- | --- |
| `scripts/build_fs_region_stable_controls.py` | `8db6d846feccbb3efa85b440c6d979d287612120cb7b3ab32a5f6b368a831dfb` |
| `configs/fs_region_stable_controls_20261015_v01.json` | `2ce799bb419a6c03f9f0b4fac8ab2907e6fe6426c4f418d3cf70e9d210a1a43f` |
| `tests/test_fs_region_stable_controls.py` | `b9f22ce351164c1e4efa436d9f0101cdf04501c8c39ba07bc37a79334d338e90` |
| `logs/reviews/FS.CTRL.01_stable_controls_20261015_verify.py` | `6d7bcc75d7666f77582201badba7c0af9e56a80efc82317dbe3338d05e2d4f64` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/REPORT.md` | `ce0941338f54f6a8b1f9671238483499edd1937e053ef90d158af525b72be84e` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/case_coverage.tsv` | `399bc50596d5d5c447cc300123503a0c02d71bf5396969e5addf432e57689d33` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/endpoint_provenance.tsv` | `12f2549ff0728fe6a1ecbf875cc6274e9682cfd8900472464b660bc1973d8818` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/evaluation_eligibility.tsv` | `cfb7af57b298b27944935969eeec346b7db72a627e6885c64aba669bee2603fe` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/labels_primary.tsv.gz` | `612d2b3abd920af2ecf66cea88872f42025c5170cda5f9c443a754d4e83be376` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/labels_strict_sensitivity.tsv.gz` | `22ad3473f7d643ba902bbb54223dbac34257212362d63ceb150096cdb7a559be` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/pair_residue_metrics.tsv.gz` | `8fbd993e46d45c30c4b7ab62b324664578bf25a497f39c5766b92c6c242bdef1` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/protocol.json` | `2ce799bb419a6c03f9f0b4fac8ab2907e6fe6426c4f418d3cf70e9d210a1a43f` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/qc.json` | `7f36189878cae140801b0f246015472b56e1b3ddde3c53a7e3babbd981a9fcec` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/run.json` | `8be348d66dd5ca1a15d8ee9b5ce9cc2fa9843db4d3404c6e779f8f5f4b4a0715` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/stable_reference_intervals.tsv` | `25d9389339ea69ae119a8b0ab02d2d12c16e9900e5251ef4df5d360c4b51a863` |
| `results/benchmark_fixed/20261015/fs_region_stable_controls_v04/structural_reference_evidence.tsv` | `cc989bc2871984f5012a8c2380953b2b022fccdeea0405dc3536dddb90bbaf68` |

## Closure-record micro-review

Reviewed at 2026-10-15 22:40 Asia/Shanghai. PASS. The updated FS.CTRL.01 task report, only the FS.CTRL.01 section at the end of TODO.md, and the two FS.CTRL.01 decision records were checked. They correctly identify v04 as the sole accepted annotation product, retain v01/v02/v03 as intermediate records, preserve the 848/580 negative counts and 124 inherited positives, and limit acceptance to structure-derived labels. No model-scoring or public-release acceptance is asserted. The task report correctly explains that the immutable run directory retains its generation-time pending-review snapshot; acceptance is recorded separately.

All reviewed scientific output hashes, builder hash, configuration hash, and test-file hash remain unchanged from the PASS review. The final task report SHA256 is `6e9ea1e80ded16ee2cc14cb8095e83fd998b7b7e50a6caad5586e319c97c1442`. TODO and decision records are consistent with this scoped PASS; no remaining closure-record issue was found.

## Administrative timestamp correction micro-review

Reviewed at 2026-10-15 22:42 Asia/Shanghai. PASS, limited to the new metadata-correction JSON and final task-report paragraph. Independently checked the original configuration with os.stat(): local filesystem st_birthtime is 1792074413.7632504, equivalent to 2026-10-15T22:26:53.763250+08:00. The correction records this workspace filesystem creation time to whole-second precision and identifies the original nominal created_at=22:24:00 discrepancy. This filesystem observation is administrative provenance, not an external preregistration timestamp. Keeping the reviewed configuration bytes immutable and recording the correction separately is appropriate. Scientific output and builder hashes remain unchanged; no geometric recomputation or rebuild was required.

Metadata correction SHA256: `171e7310c947b969a6669579f93979a2ef26dcb6651ab6521273ab3873cc71db`. Final task report SHA256, superseding the preceding micro-review hash: `92a268d3513742806cc89fbe9c4ee0074fd45f169730c213acfaef630e9d085e`. The scoped scientific PASS remains in force.
