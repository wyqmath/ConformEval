# ConformEval data and table guide

This page is the single index for the public tabular files. Each row below is a
logical table: files with different labels, masks, folds, or provenance are
kept separate because they answer different questions. The ESM-2 8M arm
prediction slices were redundant copies of the combined prediction table and
have been removed; filter the combined table by `arm` instead.

## Reading files

The package readers accept UTF-8 TSV and UTF-8 gzip-compressed TSV, including an
optional BOM and either LF or CRLF line endings. From the repository root:

```bash
zcat data/base/task_samples.tsv.gz | head
head data/base/case_modules.tsv
```

For Python, use the package reader so that compression, encoding, and malformed
input checks are consistent with scoring:

```python
from pathlib import Path
from conformeval.evaluate_predictions import read_table

rows = read_table(Path("data/base/task_samples.tsv.gz"))
print(rows[0])
```

`read_table()` returns a list of dictionaries with string values. Convert
coordinates, labels, masks, folds, and numeric scores explicitly at the point
of use. FASTA is read with `conformeval.readers.parse_fasta()` and JSON manifests
with the standard `json` module.

## Frozen task data

All paths in this section are relative to `data/`.

| File | Purpose | Key fields |
|---|---|---|
| `task_catalog.tsv` | Registry of all 13 tasks, their units, eligibility, metrics, result links, and data gaps | `task_id`, `family`, `input_unit`, `target`, `evaluation_unit`, `eligible_units`, `scored_units`, `outer_folds`, `status`, `primary_metric`, `result_artifact` |
| `base/dataset_manifest.json` | Byte and schema manifest for the base snapshot | `data_version`, `files`, `sequences`, `coordinates`, `residue_label_masks`, `source_data_lock_sha256` |
| `base/sequences.fasta` | Unique sequence inputs addressed by stable sequence identifiers | FASTA identifier and amino-acid sequence |
| `base/task_samples.tsv.gz` | Main task-level sample registry consumed by pooled evaluation | `task_id`, `sample_id`, `unit_type`, `sequence_ids`, `sequence_sha256s`, `lengths`, `binding_group_id`, `learning_label`, `learning_mask`, `label`, `mask`, `data_eligible`, `feature_ready`, `score_ready` |
| `base/sample_sequences.tsv.gz` | Expands each task sample to its sequence or endpoint identities | `task_id`, `sample_id`, `endpoint_role`, `sequence_id`, `sequence_sha256`, `length` |
| `base/case_modules.tsv` | Case-level status for paired and condition modules that are not pooled scoring tasks | `module`, `case_id`, `status`, `qualification`, `sample_count`, `binding_groups`, `source_ref`, `detail` |
| `base/case_samples.tsv.gz` | Case-level endpoint and state records | `module`, `case_id`, `sample_id`, `parent_case_id`, `state_id`, `sequence_id`, `sequence_sha256`, `length`, `binding_group_id`, `use_policy`, `record_status`, `feature_ready` |
| `base/case_sources.tsv` | Resolves a case `source_ref` to its source artifact and citation evidence | `source_ref`, `module`, `case_id`, `source_artifact_id`, `source_artifact_sha256`, `evidence_ids`, `doi`, `pmcid`, `other_source_references` |
| `base/sources.tsv` | Source register and redistribution status | `source_id`, `name`, `doi_or_url`, `related_tasks`, `version`, `license`, `distribution_status` |
| `base/representation_coverage.tsv.gz` | Records whether a sample has a usable representation cache | `task_id`, `sample_id`, `sequence_sha256`, `length`, `layer_set`, `representation_shape`, `feature_ready`, `coverage_source` |
| `base/position_labels_strict.tsv.gz` | Direct experimental residue labels and masks | `sample_id`, `sequence_sha256`, `position`, `label`, `mask`, `reason`, `evidence_ids`, `condition_id`, `source_decision` |
| `base/position_labels_structure_reference_candidate.tsv.gz` | Wider residue view combining direct labels with PDB-observed structural reference evidence | `sample_id`, `sequence_sha256`, `position`, `label`, `mask`, `label_source`, `direct_label`, `direct_mask`, `pdb_observed`, `pdb_mask`, `structures_observed`, `decision` |
| `base/fs_region_full_length_labels.tsv.gz` | Full-length fold-switch endpoint labels, including unknown positions | `task_id`, `target_unit_type`, `target_id`, `endpoint_sample_id`, `start`, `end`, `sequence_position`, `label`, `mask`, `label_basis`, `condition_id`, `evidence_id`, `decision`, `observed`, `mapping_status` |
| `base/folds/outer_group_holdout.tsv` | Frozen outer test-group assignments | `task_id`, `outer_fold`, `group_id`, `role` |
| `base/folds/inner_validation_groups.tsv.gz` | Frozen inner validation-group assignments within each outer training partition | `task_id`, `outer_fold`, `inner_fold`, `group_id`, `role` |
| `base/folds/fold_qc.tsv` | Fold counts and class-coverage checks | `task_id`, `stage`, `outer_fold`, `inner_fold`, `valid`, `train_samples`, `held_samples`, `train_groups`, `held_groups`, `train_classes`, `held_classes`, `warning` |

`position_labels_strict.tsv.gz` and
`position_labels_structure_reference_candidate.tsv.gz` are intentionally not
merged: the latter adds a structural-reference label view and its extra
provenance columns. Likewise, the direct and reference masks must not be
treated as interchangeable targets.

## Accepted disorder-reference task

These files are relative to `data/disorder_reference/` and are consumed by
`conformeval-evaluate-residues`.

| File | Purpose | Key fields |
|---|---|---|
| `manifest.json` | Hashes and sizes for the accepted disorder bundle | `artifact_id`, `files` |
| `run_manifest.json` | Input and fold lineage for bundle construction | `configuration`, `input_sha256`, `output_sha256`, `status` |
| `selected_samples.tsv` | One row per selected protein record | `sample_id`, `sequence_id`, `sequence_sha256`, `length`, `binding_group_id`, `positive_positions`, `negative_positions`, `unknown_positions`, `label_view`, `use_policy` |
| `selected_labels.tsv.gz` | One row per labeled residue used by the residue scorer | `sample_id`, `sequence_id`, `sequence_sha256`, `position`, `label`, `label_source`, `direct_label`, `direct_mask`, `pdb_observed`, `decision` |
| `outer_fold_roles.tsv` | Outer grouped roles | `outer_fold`, `sample_id`, `binding_group_id`, `role` |
| `inner_fold_roles.tsv` | Inner grouped roles | `outer_fold`, `inner_fold`, `sample_id`, `binding_group_id`, `role` |
| `fold_qc.tsv` | Residue-task fold counts and overlap checks | `outer_fold`, `train_samples`, `test_samples`, `train_groups`, `test_groups`, `group_overlap`, `inner_validation_covered_train_samples` |
| `inputs/resid_index.tsv` | Residue coordinates available to the archived extraction run | `name`, `resid_indices` |

## Fold-switch region controls

These files are relative to `data/fs_region_controls/`. They are annotation and
coverage tables; the bundle does not claim a completed model score.

| File | Purpose | Key fields |
|---|---|---|
| `manifest.json` | Hashes and sizes for the structural-control bundle | `artifact_id`, `files` |
| `labels_primary.tsv.gz` | Primary endpoint residue labels with unknowns masked | `case_id`, `endpoint_sample_id`, `sequence_id`, `sequence_sha256`, `binding_group_id`, `position_1based`, `label`, `mask`, `label_basis`, `evidence_id`, `stable_segment_id`, `decision_reason` |
| `labels_strict_sensitivity.tsv.gz` | Strict negative sensitivity view of the same endpoint frame | Same identity and coordinate fields as primary labels, with strict mask decisions |
| `pair_residue_metrics.tsv.gz` | Pairwise structural correspondence and local geometry diagnostics | `case_id`, `binding_group_id`, `position_A`, `position_B`, `auth_A`, `auth_B`, `label_A`, `label_B`, `dssp_A`, `dssp_B`, `window_rmsd_A`, `symmetric_ca_distance_score`, `coverage_A`, `coverage_B`, candidate flags |
| `stable_reference_intervals.tsv` | Intervals selected as locally preserved structural reference segments | `variant`, `case_id`, `endpoint_sample_id`, `segment_id`, `evidence_id`, `start`, `end`, `n_positions`, sequence/auth bounds |
| `structural_reference_evidence.tsv` | Evidence records supporting each interval | `evidence_id`, `case_id`, `variant`, `evidence_kind`, `protocol_sha256`, paired positions, endpoint identities, PDB IDs, structure hashes |
| `endpoint_provenance.tsv` | Endpoint sequence and source-structure provenance | `case_id`, `endpoint_sample_id`, `pdb_id`, `auth_chain`, `sequence_id`, `sequence_sha256`, `length`, `binding_group_id`, source and extracted structure hashes, method |
| `case_coverage.tsv` | Per-case correspondence and label counts | `case_id`, `binding_group_id`, `correspondence_rows`, `usable_identity_pairs`, `inherited_positive_rows`, `primary_negative_rows`, `strict_negative_rows`, `mapping_exclusions` |
| `evaluation_eligibility.tsv` | Explicit status for future scoring eligibility | `case_id`, `binding_group_id`, `positive_rows`, `negative_rows`, `status` |
| `run.json` | Construction configuration, dependencies, input hashes, and output hashes | `configuration`, `dependencies`, `input_sha256`, `output_sha256`, `status` |

## Results tables

`results/results_index.tsv` is the single cross-task result index. It contains
one row per task/model/arm and points to the configuration, environment,
prediction, and metric files. Its fields are `task_id`, `status`, `model`,
`model_id`, `model_revision`, `arm`, `arm_role`, `configuration_id`,
`configuration_path`, sample and group counts, `primary_metric`, `direction`,
`primary_value`, secondary metrics, environment fields, and prediction/score
paths.

The pooled prediction tables use one schema per model family. The 650M
reference keeps its arm identity in every row and the 8M runs use the same
combined-table pattern for both pooled tasks:

| Directory | File | Purpose | Key fields |
|---|---|---|---|
| `results/esm2_reference/` | `predictions.tsv.gz` | Combined out-of-fold predictions for all 650M arms | `task_id`, `sample_id`, `outer_fold`, `binding_group_id`, `score`, `arm`, task labels and masks |
| `results/esm2_8m/{knot,fs}/` | `predictions.tsv.gz` | Combined out-of-fold predictions for all 8M arms; filter `arm` to select one | `task_id`, `sample_id`, `outer_fold`, `score`, `arm` |
| `results/esm2_reference/` | `inner_selection.tsv` | 650M inner-fold selection records | `outer_fold`, `inner_fold`, `arm`, `layer`, `C`, `selection_metric`, `metric_name` |
| `results/esm2_reference/` | `outer_selections.tsv` | 650M final choice for each outer fold | `task_id`, `outer_fold`, `arm`, `selected_layer`, `selected_C`, `inner_selection_metric` |
| `results/esm2_8m/{knot,fs}/` | `inner_selection.tsv` | 8M inner-fold selection records | `outer_fold`, `inner_fold`, `layer`, `predictor`, `C`, `validation_metric`, `validation_rows`, `aggregation_weight` |
| `results/esm2_8m/{knot,fs}/` | `outer_selections.tsv` | 8M final choice for each outer fold and arm | `outer_fold`, `arm`, `layer`, `predictor`, `C`, `inner_metric`, `train_rows`, `test_rows` |

The following audit tables are present only in the 650M pooled reference:
`fold_roles.tsv.gz` (`task_id`, `outer_fold`, `stage`, `inner_fold`,
`sample_id`, `binding_group_id`, `role`), `fold_metrics.tsv` (`task_id`,
`outer_fold`, `arm`, `selected_layer`, `selected_C`, train/test counts and
outer metrics), `case_ranks.tsv` (held-positive rank and recall fields), and
`warnings.tsv` (fold, arm, phase, warning type and message). They are kept
separate because the 8M runs record the equivalent selection information in
their own schema.

The residue result directories contain `protein_metrics.tsv` with equal-weight
per-protein aggregation (`task_id`, `sample_id`, `binding_group_id`,
`outer_fold`, `arm`, labeled/positive/negative counts, `average_precision`,
`auroc`) and, for the full disorder run, `residue_predictions.tsv.gz` with
`task_id`, `sample_id`, `binding_group_id`, `outer_fold`, `arm`,
`position_1based`, `label`, `mask`, and `score`.

Each run also has JSON metadata (`metrics.json`, `model_metadata.json`, and
`run.json`) and fitted parameter files. The arm-specific enriched scorer tables
under `results/esm2_8m/{knot,fs}/scores/<arm>/scored_predictions.tsv.gz` remain
separate because they include the task labels, masks, and binding groups needed
to audit a scorer output. Their columns are `task_id`, `sample_id`,
`outer_fold`, `binding_group_id`, `score`, `learning_label`, `learning_mask`,
`biological_target`, `biological_mask`, `label`, and `mask`. They are not
interchangeable with the compact prediction table.

The disorder-reference result directory additionally contains
`inner_selection.tsv` (inner-fold reader choices), `outer_fits.tsv` (outer-fold
fit counts and selected hyperparameters), `protein_metrics.tsv` (per-protein
residue scores), and `residue_predictions.tsv.gz` (the full residue prediction
table with `task_id`, `sample_id`, `binding_group_id`, `outer_fold`, `arm`,
`position_1based`, `label`, `mask`, and `score`).

To read one arm from a combined prediction file:

```python
from pathlib import Path
from conformeval.evaluate_predictions import read_table

rows = read_table(Path("results/esm2_8m/knot/predictions.tsv.gz"))
nested = [row for row in rows if row["arm"] == "PLM-nested"]
```

## Examples

The small files under `examples/mini/` have the same schemas as their full
counterparts and are intended for interface tests only. `examples/predictions/`
contains the two submission tables: `knot.tsv` has
`task_id`, `sample_id`, `outer_fold`, `score`; `residues.tsv` adds
`position_1based`. Both may be read with `read_table()` and passed directly to
the scoring commands described in `docs/prediction_submission.md`.

## What was consolidated

The 8M task directories previously contained seven files named after the arms,
each repeating a four-column slice of `predictions.tsv.gz`. Those 14 slices are
now removed. The combined table is the only compact prediction entry point;
the arm name is an explicit column, so no data or score is lost. Files with
different scientific targets, masks, fold roles, or provenance remain separate
and are documented above rather than silently concatenated.
