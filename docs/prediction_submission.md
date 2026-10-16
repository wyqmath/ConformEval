# Prediction submission and scoring

This page documents the pooled-task scorer. The third scored task has a separate residue interface in [residue_evaluation.md](residue_evaluation.md).

The scorer accepts complete out-of-fold predictions for the two pooled tasks with frozen scoring folds: `KNOT-PRESENCE` and `FS-L2-grouped-PU-ranking`. It does not fit a reader or inspect model representations. The model owner supplies one score per eligible sample, and the scorer joins those scores to the benchmark's exact sample IDs, labels, binding groups, and outer-fold assignments.

## Prediction table

Provide a tab-separated file, optionally gzip-compressed, with exactly these columns:

```text
task_id    sample_id    outer_fold    score
```

Each eligible sample must appear once. `outer_fold` must match the frozen assignment for its binding group. Scores must be finite numbers; they may be logits or probabilities as long as larger values mean stronger evidence for the positive class or higher ranking priority. Do not include labels in this input file.

For FS-L2, score each positive endpoint pair as one case using the benchmark protocol, and score each background protein as one unlabeled reference unit. The `learning_label=0` background is an operational PU ranking reference. Its biological target remains unknown and is never evaluated as a verified negative.

## Required model metadata

Supply a JSON file alongside the predictions:

```json
{
  "model_id": "organization/model-name",
  "model_revision": "immutable-commit-or-checkpoint-revision",
  "prediction_method": "nested balanced logistic regression",
  "prediction_code_sha256": "64-character-lowercase-sha256",
  "training_protocol": "frozen_nested_grouped_oof_v1",
  "candidate_layers": [17, 33],
  "pooling": "mean_over_residues_without_special_tokens"
}
```

The scorer records these fields as model-owner provenance. It checks the required fields and hash syntax, but it cannot independently prove which weights, training examples, or selection records produced the supplied scores. Preserve the inference and reader-fitting logs with the run for that audit.

The included nested-reader command follows the frozen outer and inner group assignments, fits preprocessing on each current training partition, selects KNOT-PRESENCE readers by mean inner-fold AUROC, and selects FS-L2 readers by the held-positive rank percentile within each inner validation pool (lower is better). It fits the selected reader on the outer training partition and scores only that fold's outer test pool. Do not compare positive cases across independently fitted outer-fold models.

The ESM-2 reference uses candidate layers 17 and 33, `C_grid=[0.1, 1.0]`, balanced logistic regression with the `liblinear` solver, `max_iter=5000`, seed 2026, and a `StandardScaler` fitted within each training partition. The current fixed protocol accepts candidate layers 17 and 33. A future protocol may declare other layers before scoring; never choose candidate layers from outer-fold results.

## Fit the frozen readers from representations

Use the representation archive contract in [model_input.md](model_input.md) to prepare an exact-sequence cache and manifest. Install the evaluation extra and pass explicit paths:

```bash
python -m pip install ".[evaluation]"
confomeval-fit-readers \
  --data-dir data/v0.3.1 \
  --protocol configs/scoring_protocol.json \
  --feature-manifest /path/to/representation_manifest.json \
  --output-dir /path/to/new_reader_run
```

The output directory must not exist. The runner writes a `RUNNING` record first and leaves a `FAILED` record with the error if a run stops. A successful run records hashes, versions, timestamps, warnings, inner and outer selections, fitted coefficients, nearest-neighbor references, seven complete OOF prediction tables, and scored metrics. Existing outputs are never overwritten. The runner can fit supplied model representations; it does not download weights or run model inference.

If you already have predictions, skip reader fitting and use `confomeval-evaluate` below.

## Run the scorer

Install the evaluation dependency from the benchmark directory, then provide explicit paths:

```bash
python -m pip install ".[evaluation]"
confomeval-evaluate \
  --data-dir data/v0.3.1 \
  --protocol configs/scoring_protocol.json \
  --predictions /path/to/predictions.tsv.gz \
  --model-metadata /path/to/model_metadata.json \
  --output-dir /path/to/new_scoring_run
```

The output directory must not exist. The scorer writes `metrics.json`, `scored_predictions.tsv.gz`, and `run.json`. The run record stores input, dataset, scorer, and output hashes; dependency versions; timestamps; row counts; and the task metrics. A complete prediction file is required for both tasks, so partial task runs cannot be presented as a complete ConformEval result.

## Reported metrics

KNOT-PRESENCE reports pooled outer-fold AUROC as the primary metric and average precision as a secondary metric. Its 95% interval resamples binding groups 2,000 times and is conditional on the fixed out-of-fold predictions.

FS-L2 reports each held positive case's rank percentile within its own outer test pool, then averages across the ten positive cases. It also reports Recall@1, @10, @100, @1,000, and @10,000. Its 95% interval resamples the ten positive binding groups 2,000 times and is conditional on the fixed predictions. Background proteins are ranking references, not verified biological negatives.

These are retrospective grouped-development results. The frozen group split audit and the independent historical-family audit answer different questions; this score file does not certify family-independent or prospective generalization.
