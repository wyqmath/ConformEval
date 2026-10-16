# Disorder reference evaluation

The task ID remains `DISORDER-STRUCTURE-REFERENCE-CANDIDATE` to preserve its registered identity. The separate v1 reference evaluation is now independently accepted. This does not rename or replace `DISORDER-RES`.

## Inputs and labels

`data/disorder_reference_v1/selected_samples.tsv` contains 1,279 protein records, 1,275 unique sequences and 1,029 binding groups. `selected_labels.tsv.gz` contains 391,718 scoring rows: 73,766 direct disorder positives, 35 direct order negatives and 317,917 PDB-observed reference negatives. A further 156,160 unknown positions are excluded. Positions are one-based coordinates in `inputs/sequences.fasta`.

Four pairs of records share identical sequences. Twelve sequence coordinates carry different labels across those records. Records remain distinct; all copies share a binding group and no exact sequence crosses outer or inner roles. The recorded interpretation follows the source evidence. There are 391,189 unique sequence-coordinate vectors behind the 391,718 scoring rows.

## Score a new prediction table

Provide one score for every selected coordinate, using these TSV columns:

```text
task_id    sample_id    position_1based    outer_fold    score
```

The optional `arm` column names multiple complete prediction arms. Without it, the arm name is the metadata `model_id`. Higher scores indicate disorder. Unknown or extra coordinates, incomplete coverage, duplicates, non-finite scores and incorrect outer folds are rejected. Labels are joined from the task data. Prediction labels are unnecessary; if present in a reference table, they must match.

The metadata JSON requires `model_id`, `model_revision`, `prediction_method`, `prediction_code_sha256` and `training_protocol`. Set the latter to `frozen_nested_grouped_oof_v1`. The code hash is a lowercase SHA256 fingerprint of the prediction implementation. Multiple arms must be declared in an `arms` list. This metadata records the model owner's account of inference and training; retain the actual logs for review.

```bash
confomeval-evaluate-residues \
  --task-dir data/disorder_reference_v1 \
  --predictions /path/to/residue_predictions.tsv.gz \
  --model-metadata /path/to/model_metadata.json \
  --output-dir runs/my-residue-scores
```

The scorer writes per-protein AP/AUROC, mean per-protein AP as the primary metric, secondary pooled AP, and 2,000 binding-group bootstrap intervals conditional on fixed predictions. Every sample has both classes. It saves hashes, dependency versions, timestamps and output records.

## Replay the archived predictions

```bash
confomeval-evaluate-residues \
  --task-dir data/disorder_reference_v1 \
  --predictions results/disorder_reference_esm2_v1/residue_predictions.tsv.gz \
  --reference \
  --output-dir runs/disorder-reference-replay
```

The archived arms are ESM-2 layers 11, 23 and 33, a 41-feature amino-acid/local-context baseline, and layer-33 MLP readers with seeds 13, 42 and 2026. All seeds are reported; outer scores do not select a seed.

## Refit the archived readers

Use the exact ESM-2 revision and residue archive format in [`configs/disorder_reference_v1/extraction_request.json`](../configs/disorder_reference_v1/extraction_request.json), [`evaluation_protocol.json`](../configs/disorder_reference_v1/evaluation_protocol.json) and [model_input.md](model_input.md). Model weights and 1,275 representation archives are external inputs. Prepare an `extract_manifest.tsv` and NPZ archives under one cache directory; exact sequence hashes, full coverage and requested coordinates are validated.

```bash
python -m pip install ".[reference-residue-training]"
confomeval-fit-reference-residues \
  --task-dir data/disorder_reference_v1 \
  --config-dir configs/disorder_reference_v1 \
  --cache-dir /path/to/esm2-residue-cache \
  --output-dir runs/disorder-reader-refit \
  --device cuda
```

This command preserves the archived ESM-2 model, revision, layers, feature dimension, C grid, folds and capacity diagnostic. It checks the fixed task configuration fingerprints before fitting; changing the scientific protocol requires a separately versioned configuration and runner support. It is not a configurable training runner for arbitrary encoders. It accepts CPU for reader training, although the archived cache metadata records CUDA extraction. The source research runner's hash is preserved separately from this importable package module. A fresh run needs its own review before it becomes an accepted result.

The original reference job ran for about 2 hours 21 minutes on its recorded cluster environment after representation extraction. Hardware and numerical environment affect refitting time and floating-point equality. The local reproduction command rescores all archived predictions; it does not rerun encoder inference or refit the full inner-selection trajectory.
