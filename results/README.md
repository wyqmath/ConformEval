# Results and task status

This index covers complete biological task evaluations. **39 saved arm records** are reported with their original experiment settings; baseline repeats are retained by source run. **33 task/model combinations are 未测试 / Not tested**. The three-model synthetic smoke tests are separately recorded in [model validation](../docs/model_validation.md).

## Configuration and environment keys

| Key | Configuration source | Reader/extraction settings | Recorded environment | Evidence |
|---|---|---|---|---|
| 650M-pooled | [configuration](../results/esm2_existing_v0.3/configuration.json) | Mean sequence vectors; logistic C={0.1,1}; nested layers={17,33} or arm-fixed layer; nearest uses representation; train-only scaling | CPU reader training; {"numpy":"2.5.3","python":"3.14.6","scikit_learn":"1.9.1","scipy":"1.18.1"} | [actual run/environment](../results/esm2_existing_v0.3/provenance.json) |
| 650M-residue | [configuration](../configs/disorder_reference_v1/evaluation_protocol.json) | Residue logistic C={0.1,1}; fixed layer indicated by arm; 5x3 grouped folds; train-only scaling | Run device=cuda; MLP on GPU, sklearn linear readers on CPU; {"cuda":"12.1","numpy":"1.26.4","python":"3.11.16","scikit_learn":"1.6.1","scipy":"1.13.1","torch":"2.5.1+cu121"} | [actual run/environment](../results/disorder_reference_esm2_v1/run.json) |
| 8M-knot | [configuration](../results/esm2_8m_v1/knot/model_metadata.json) | Layers={3,6}; mean sequence pooling; logistic C={0.1,1}; train-only scaling; frozen grouped inner selection | CPU, macOS arm64; {"numpy":"2.5.4","python":"3.14.5","scikit-learn":"1.9.1","scipy":"1.18.1"} | [actual run/environment](../results/esm2_8m_v1/knot/run.json) |
| 8M-fs | [configuration](../results/esm2_8m_v1/fs/model_metadata.json) | Layers={3,6}; mean sequence pooling; logistic C={0.1,1}; train-only scaling; frozen grouped inner selection | CPU, macOS arm64; {"numpy":"2.5.4","python":"3.14.5","scikit-learn":"1.9.1","scipy":"1.18.1"} | [actual run/environment](../results/esm2_8m_v1/fs/run.json) |
| 8M-residue | [configuration](../results/esm2_8m_v1/residue/model_metadata.json) | Per-residue vectors, no sequence pooling; candidate layers={3,6} for nested or arm-fixed layer; logistic C={0.1,1}; train-only scaling; frozen grouped inner selection | CPU, macOS arm64; {"numpy":"2.5.4","python":"3.14.5","scikit-learn":"1.9.1","scipy":"1.18.1"} | [actual run/environment](../results/esm2_8m_v1/residue/run.json) |

The old 650M pooled `PLM-composition` arm includes log-length plus 20-AA frequencies (21 features), whereas the 8M arm uses 20-AA frequencies only. Their names are inherited from the original runs; the different feature sets are stated in each TSV row.

For 8M the exact extraction/reader configurations are embedded in each model metadata file. For 650M residue MLP, the configuration additionally specifies layer 33, width 64, ReLU, dropout 0.1, Adam 0.001, batch 8,192 and inner early stopping. Seeds 13, 42 and 2026 are all reported. Pooled 650M environment metadata describes the CPU downstream reader run; GPU details for its earlier pre-extraction are not recorded there.


## KNOT-PRESENCE

| Model / baseline | Arm | Config key | Scored data (available data if untested) | Metric | Value | Environment key | Prediction file | Metric source |
|---|---|---|---|---|---|---|---|---|
| ESM-2 650M | PLM-nested | 650M-pooled | 734 sample | AUROC | 0.868741 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| ESM-2 650M | PLM-hidden17 | 650M-pooled | 734 sample | AUROC | 0.801892 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| ESM-2 650M | PLM-hidden33 | 650M-pooled | 734 sample | AUROC | 0.868741 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| Sequence length + logistic | PLM-length | 650M-pooled | 734 sample | AUROC | 0.675945 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| Length + AA composition + logistic | PLM-composition | 650M-pooled | 734 sample | AUROC | 0.819491 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| ESM-2 650M | PLM-nearest | 650M-pooled | 734 sample | AUROC | 0.672987 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| Training-label prior | PLM-training-prior | 650M-pooled | 734 sample | AUROC | 0.494323 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_KNOT-PRESENCE.json) |
| AA composition + logistic | PLM-composition | 8M-knot | 734 sample | AUROC | 0.785140 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| ESM-2 8M | PLM-hidden3 | 8M-knot | 734 sample | AUROC | 0.821802 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| ESM-2 8M | PLM-hidden6 | 8M-knot | 734 sample | AUROC | 0.832955 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| Sequence length + logistic | PLM-length | 8M-knot | 734 sample | AUROC | 0.675945 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| ESM-2 8M | PLM-nearest | 8M-knot | 734 sample | AUROC | 0.724348 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| ESM-2 8M | PLM-nested | 8M-knot | 734 sample | AUROC | 0.814237 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| Training-label prior | PLM-training-prior | 8M-knot | 734 sample | AUROC | 0.494323 | 8M-knot | [predictions](../results/esm2_8m_v1/knot/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/knot/metrics.json) |
| ProtT5 XL encoder | 未测试 | 未测试 | 734 sample | AUROC; average precision secondary | 未测试 | 未测试 | 未测试 | 未测试 |

Shared prediction files are filtered by `task_id` and `arm`; these exact filters are in `results_index.tsv`.


## FS-L2-grouped-PU-ranking

| Model / baseline | Arm | Config key | Scored data (available data if untested) | Metric | Value | Environment key | Prediction file | Metric source |
|---|---|---|---|---|---|---|---|---|
| ESM-2 650M | PLM-nested | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.476592 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| ESM-2 650M | PLM-hidden17 | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.529372 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| ESM-2 650M | PLM-hidden33 | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.316165 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| Sequence length + logistic | PLM-length | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.711773 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| Length + AA composition + logistic | PLM-composition | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.498057 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| ESM-2 650M | PLM-nearest | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.500728 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| Training-label prior | PLM-training-prior | 650M-pooled | 8879 ranking_item | mean_held_case_rank_percentile | 0.500564 | 650M-pooled | [predictions](../results/esm2_existing_v0.3/predictions.tsv.gz) | [metrics](../results/esm2_existing_v0.3/metrics_FS-L2-grouped-PU-ranking.json) |
| AA composition + logistic | PLM-composition | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.480346 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| ESM-2 8M | PLM-hidden3 | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.220994 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| ESM-2 8M | PLM-hidden6 | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.355781 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| Sequence length + logistic | PLM-length | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.711773 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| ESM-2 8M | PLM-nearest | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.500737 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| ESM-2 8M | PLM-nested | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.304406 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| Training-label prior | PLM-training-prior | 8M-fs | 8879 ranking_item | mean_held_case_rank_percentile | 0.500564 | 8M-fs | [predictions](../results/esm2_8m_v1/fs/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/fs/metrics.json) |
| ProtT5 XL encoder | 未测试 | 未测试 | 8879 ranking_item | mean held-case rank percentile; Recall@k | 未测试 | 未测试 | 未测试 | 未测试 |

Shared prediction files are filtered by `task_id` and `arm`; these exact filters are in `results_index.tsv`.


## DISORDER-STRUCTURE-REFERENCE-CANDIDATE

| Model / baseline | Arm | Config key | Scored data (available data if untested) | Metric | Value | Environment key | Prediction file | Metric source |
|---|---|---|---|---|---|---|---|---|
| ESM-2 650M | ESM2-L11 | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.722277 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| ESM-2 650M | ESM2-L23 | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.732330 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| ESM-2 650M | ESM2-L33 | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.789856 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| ESM-2 650M | MLP-h33-seed13 | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.802670 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| ESM-2 650M | MLP-h33-seed2026 | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.803926 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| ESM-2 650M | MLP-h33-seed42 | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.802767 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| AA + position + local context + logistic | traditional-AA-context | 650M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.585906 | 650M-residue | [predictions](../results/disorder_reference_esm2_v1/residue_predictions.tsv.gz) | [metrics](../results/disorder_reference_esm2_v1/metrics.json) |
| ESM-2 8M | PLM-hidden3 | 8M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.740719 | 8M-residue | [predictions](../results/esm2_8m_v1/residue/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/residue/metrics.json) |
| ESM-2 8M | PLM-hidden6 | 8M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.744395 | 8M-residue | [predictions](../results/esm2_8m_v1/residue/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/residue/metrics.json) |
| ESM-2 8M | PLM-nested | 8M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.738289 | 8M-residue | [predictions](../results/esm2_8m_v1/residue/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/residue/metrics.json) |
| AA + position + local context + logistic | traditional-AA-context | 8M-residue | 1279 protein_record; 391,718 rows | mean_per_protein_AP | 0.585906 | 8M-residue | [predictions](../results/esm2_8m_v1/residue/predictions.tsv.gz) | [metrics](../results/esm2_8m_v1/residue/metrics.json) |
| ProtT5 XL encoder | 未测试 | 未测试 | 1279 protein | mean per-protein average precision | 未测试 | 未测试 | 未测试 | 未测试 |

Shared prediction files are filtered by `task_id` and `arm`; these exact filters are in `results_index.tsv`.


## Tasks without model evaluation

| Task | Input | Current data | Metric | Status / data needed |
|---|---|---|---|---|
| `KNOT-TYPE` | Full sequence | 65 eligible: 3_1=57, 4_1=6, 5_2=2 | Macro-F1 on locked classes ↑ | 未测试. Required 5_1 class is absent; 4_1 and 5_2 have only 6 and 2 examples. No scoring folds |
| `DISORDER-RES` | Sequence + residue coordinates | 2,279 registered proteins; 1 eligible dual-class binding group | Mean protein AP ↑ | 未测试. Additional proteins/groups with direct experimental order and disorder labels; locked protocol requires ≥4 dual-class groups. No scoring folds |
| `FS-L1-PAIRED` | Paired experimental sequences and structures | 16 candidates; 10 accepted cases / 20 endpoints | Case evidence and input-identity audit; no model metric | 未测试. Six pending candidates lack accepted endpoint inputs/evidence in the frozen sequence frame. Case report; no prediction target/scoring folds defined |
| `FS-REGION` | Sequence + switching-region coordinates | 10 cases / 20 endpoints; 124 positive + 848 reference-negative rows (580 strict); 3 dual-class cases | Label coverage; no model metric defined | 未测试. Switching-region positives cover 3/10 cases; 4,746 endpoint positions remain unknown. Annotations available; scoring protocol/folds not defined; model evaluation not run |
| `FS-L3-MATCHED` | Sequences of cases and matched controls | 0 qualified matched case-control pairs | AUROC/AP after qualification ↑ | 未测试. Experimentally qualified matched controls with usable biological labels. Control-label qualification pending |
| `BPTI-REDOX` | Same-variant sequence + redox conditions | 0 qualified experiment pairs | Not defined after qualification stop | 未测试. Qualified endpoint measurements and matched-state evidence. Prespecified qualification stop |
| `PYP-STATE` | Sequence + illumination condition | 2 state records | Case/input audit; no prediction metric defined | 未测试. Additional independent groups with matched illumination/state measurements. Case report; no independent scoring folds |
| `RNASE-ASSEMBLY` | Sequence + assembly condition | 4 state records | Case/input audit; no prediction metric defined | 未测试. Independent same-batch assembly pairs with measured conditions; current pairs carry cross-batch warnings. Case report; no independent scoring folds |
| `DISORDER-CONDITION-pH` | Same sequence/region + pH | 7 state records | Case/input audit; no prediction metric defined | 未测试. Additional qualified independent condition pairs; candidate records include cross-publication and region mismatches. Case report; no independent scoring folds |
| `DISORDER-CONDITION-temperature` | Same sequence/region + temperature | 2 state records | Case/input audit; no prediction metric defined | 未测试. Additional independent groups with matched two-temperature experimental labels. Case report; no independent scoring folds |

The [72-row TSV](results_index.tsv) supplies explicit rows for all three models on all 13 tasks, with every tested arm, full-precision primary and secondary values, metric direction, config and environment paths, sample/row units, prediction filters and data gaps. `model_id=none` identifies sequence-only baselines; `experiment_model` identifies the source run. An untested row records currently eligible data, **zero scored rows**, and `未测试` for configuration, environment, predictions and scores.

Condition modules and FS-L1 are case/evidence records; their prediction metrics and independent scoring folds are not defined. FS-REGION has paired structural annotations but no model evaluation. KNOT-TYPE lacks required 5_1 examples; strict DISORDER-RES has only one dual-class group against the locked minimum of four. These statuses are recorded separately from available data.

All scores are retrospective grouped-development results. Confidence intervals are conditional on saved out-of-fold predictions. See the [source and data definitions](../docs/data.md), [all-task dispositions](task_dispositions.md), and [resource measurements](../docs/resources.md).
