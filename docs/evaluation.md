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

The released splits support retrospective grouped-development evaluation. Registered group isolation
is checked, but complete biological-family independence has not been established;
the source group audit retains 37 unresolved PDB-key gaps. Condition modules
and other unscored tasks remain separately listed in the task catalog.

## Local tests

```bash
python -m unittest discover -s tests -v
```

This runs interface, failure-case and scientific-input checks without downloading models. Checks cover sample identities, label masks, residue coordinates, prediction completeness, and group isolation.
