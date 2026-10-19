# Evaluation

| Task | Unit | Primary metric | Direction | Outer / inner folds |
|---|---|---|---|---|
| Knot presence | Whole-chain sample | AUROC | Higher | 5 / 3 |
| Fold-switch ranking | Ranking item | Mean held-positive rank percentile | Lower | 10 / 3 |
| Residue disorder | Protein record | Mean protein average precision | Higher | 5 / 3 |

The released group assignments define every training/validation/test partition.
All models share samples, groups and scoring implementations. Preprocessing
is fit on training data only; inner validation selects layers and reader
hyperparameters. Saved predictions come from held-out outer groups.

## Scoring rules

Knot AUROC uses all eligible chain labels. Fold-switch ranking evaluates each
held-out confirmed case against the unlabeled background; its training labels
are operational PU labels, not verified biological negatives. Disorder uses
only known reference positions, scores each eligible protein record and then
averages protein AP. Unknown positions never become negative labels.

The scorers reject incomplete/duplicate predictions, sample or sequence identity
mismatches, incorrect folds, invalid coordinates and invalid label masks.
Confidence intervals use group resampling of fixed out-of-fold predictions.
They do not include uncertainty from retraining the encoder or selecting a
new dataset.

## Compare results

The homepage reports the archived nested pooled readers, the fixed-layer 650M
disorder reference and the 8M nested disorder reader. Detailed arms, baselines,
configurations and metric precision remain in the [result index](../results/results_index.tsv).
The models' candidate layer sets and original environments differ; these
numbers are not a controlled parameter-count comparison.

These are retrospective grouped-development results. Registered group isolation
is checked, but complete biological-family independence has not been established;
the source group audit retains 37 unresolved PDB-key gaps. Condition modules
and other unscored tasks remain separately listed in the task catalog.

## Local tests

```bash
python -m unittest discover -s tests -v
```

This runs interface, failure-case and scientific-input checks without downloading models. Three checks of archived predictions are skipped when those files are absent. To run them, extract the [reference snapshot](../results/README.md#full-reference-runs) and set `CONFORMEVAL_REFERENCE_ROOT` to its directory before running the same test command. Exact run environments, predictions and parameters remain in that snapshot.
