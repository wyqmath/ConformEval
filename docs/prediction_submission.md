# External predictions and metadata

Use the base package to score predictions from any model. The same scoring implementation is used by the supplied training command.

## Pooled tasks

Supply a TSV, optionally gzip-compressed, with exactly these columns:

```text
task_id    sample_id    outer_fold    score
```

Each selected eligible sample must appear exactly once in its assigned outer test fold. Higher scores indicate knot presence or stronger fold-switch ranking. Biological labels are joined from task data. The FS background is an unlabeled PU reference, not a verified negative class. Use `--task` once or repeatedly to select pooled tasks; omit it to score both.

A runnable external example is included:

```bash
confomeval-evaluate   --task KNOT-PRESENCE --data-dir examples/mini/base   --predictions examples/predictions/knot.tsv   --model-metadata examples/predictions/metadata.json   --protocol configs/scoring_protocol.json   --output-dir runs/external-knot
```

Expected: AUROC 1.0 over 20 synthetic rows, saved to `runs/external-knot/metrics.json`. For the full benchmark point `--data-dir` to `data/base` and supply complete predictions for its frozen sample set.

## Residue task

The coordinate-level columns are:

```text
task_id    sample_id    position_1based    outer_fold    score
```

Optional `arm` names multiple complete prediction sets. Every selected, labeled coordinate must occur once per arm. Unknown coordinates, missing rows, duplicates, non-finite values, incorrect groups or folds are rejected.

```bash
confomeval-evaluate-residues   --task-dir examples/mini/residue   --predictions examples/predictions/residues.tsv   --model-metadata examples/predictions/metadata.json   --output-dir runs/external-residue --example
```

Expected: 554 synthetic residue rows and mean per-protein AP 1.0, saved in `runs/external-residue/metrics.json`. For benchmark predictions use `--task-dir data/disorder_reference` and omit `--example`.

## Metadata template

[`examples/predictions/metadata.template.json`](../examples/predictions/metadata.template.json) lists the fields to fill with your actual model identity, immutable revision, inference/training method and prediction-code SHA-256. `candidate_layers` and `pooling` describe the model representation. Set `training_protocol` to `frozen_nested_grouped_oof_v1` only when those predictions follow the supplied nested group folds. Multiple residue arms must be declared in `arms`.

The executable synthetic example has its own metadata; its zero fingerprint is a placeholder for a schema example, and its `scope` is explicitly `software_example`. Real submissions require their actual implementation fingerprints. Metadata is recorded provenance; preserve the actual extraction and training logs with the prediction file.

Scoring writes joined predictions, metrics, input/output hashes, actual versions and timestamps to the new output directory. Confidence intervals use 2,000 binding-group bootstrap replicates conditional on the supplied outer predictions. Test scores do not select a layer, predictor or hyperparameter.
