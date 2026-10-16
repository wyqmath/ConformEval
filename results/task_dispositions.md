# Task dispositions

This page explains why every registered task has a scored evaluation or is marked 未测试 in the current results index. No score entry means that the locked data or protocol did not define a valid evaluation for that task; it is not a model-performance result.

| Task | Data in the release | Disposition |
|---|---|---|
| KNOT-PRESENCE | 734 eligible samples, 283 binding groups | ESM-2 650M and 8M scored with 5 frozen outer folds. The FD.03 supplement completes the mean-representation coverage at 734/734. |
| KNOT-TYPE | 65 eligible samples: `3_1=57`, `4_1=6`, `5_2=2` | No scoring folds. The locked class set also requires `5_1`, which is absent. |
| FS-L2 grouped PU ranking | 8,879 evaluation units, 10 labeled positive cases, 10 outer folds | ESM-2 650M and 8M scored as a per-case rank within its own held-out fold pool. The rest of the background remains biologically unlabeled. |
| DISORDER-RES | Residue label tables for 2,279 proteins | No scoring folds. Only one eligible binding group contains both valid ordered and disordered labels; the locked grouped evaluation requires more. |
| DISORDER-STRUCTURE-REFERENCE-CANDIDATE | 1,279 protein records / 1,029 binding groups; 391,718 selected labeled residues | Accepted separate reference evaluation with 5 outer / 3 inner folds and seven ESM-2/traditional/capacity arms. 650M results are in `disorder_reference_esm2_v1/`, and 8M results in `esm2_8m_v1/residue/`; strict direct DISORDER-RES labels remain separate. |
| FS-L1-PAIRED | 16 candidate cases; 10 accepted cases (20 endpoint rows) | Case and input-evidence report only; no general performance score is defined. |
| FS-REGION | 10 experimental pairs / 20 endpoints; 124 inherited positives and 848 stable structure-reference negatives (580 strict); three dual-class cases contain 356 primary negatives | Accepted paired-state reference annotations in `data/fs_region_controls_v1/`. No region model evaluation or scoring folds have been run. These negatives describe locally preserved structure in the observed pair. |
| FS-L3-MATCHED | No qualified case-control pairs | No score. Matching constraints do not establish biological negative labels. |
| BPTI-REDOX | No case passes the accepted qualification rule | Prespecified qualification stop; not a negative result. |
| PYP-STATE | 2 state records | Case and input-evidence report only. |
| RNASE-ASSEMBLY | 4 records | Case and input-evidence report only. |
| DISORDER-CONDITION-pH | 7 condition records | Case report only; no independent scoring folds. |
| DISORDER-CONDITION-temperature | 2 condition records | Case report only; no independent scoring folds. |

The exact denominators, output pointers, and interpretation limits are in [`task_catalog.tsv`](../data/task_catalog.tsv). Base labels and masks are in [`data/v0.3.1/`](../data/v0.3.1/); case records are in `case_modules.tsv` and `case_samples.tsv.gz`. The accepted supplemental labels are in [`data/disorder_reference_v1/`](../data/disorder_reference_v1/) and [`data/fs_region_controls_v1/`](../data/fs_region_controls_v1/).

Condition case modules remain descriptive in this version. The new condition four-arm experiment is deferred. The separate supplemental label/result bundles update the current task state without rewriting the frozen base snapshot.

The [result summary](README.md) and [uniform result index](results_index.tsv) list actual scores, configurations, environments and prediction paths. All three complete ProtT5 task evaluations and the ten tasks without model evaluation are explicitly marked 未测试; completed synthetic examples do not fill those entries.
