# ConformEval

ConformEval tests how well protein sequence representations support prediction of topology, disorder and fold switching. It provides versioned sequence and label tables, frozen group splits, scoring tools, and audited ESM-2 reference predictions.

**Version 0.4.0 is the first benchmark snapshot in the [ConformEval repository](https://github.com/wyqmath/ConformEval).** Three tasks have reference model results. Paired-state region controls and condition case records are included as data modules. The condition four-arm experiment is deferred to a later version.

| Scored task | Evaluation data | Reference result |
|---|---|---|
| KNOT-PRESENCE | 734 samples; 283 binding groups; 5 outer × 3 inner folds | Nested ESM-2: AUROC 0.868741; AP 0.687211 |
| FS-L2 grouped PU ranking | 8,879 ranking units, including 10 positive cases; 10 outer × 3 inner folds | Nested ESM-2: mean held-case rank percentile 0.476592; lower is better |
| DISORDER-STRUCTURE-REFERENCE-CANDIDATE | 1,279 protein records; 1,029 binding groups; 391,718 labeled residues; 5 outer × 3 inner folds | ESM-2 layer 33: mean protein AP 0.789856; traditional baseline 0.585906 |

These are retrospective grouped development results. Unknown biological labels remain masked. In FS-L2, the ranking background is unlabeled. In the disorder reference task, PDB-observed positions supply a structural reference class; the strict direct experimental order/disorder task remains a separate data view.

## Install and run locally

Use Python 3.11 or newer. Clone the repository and install the package:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv
source .venv/bin/activate
python -m pip install .
python -m unittest discover -s tests -v
```

The Python import package and command prefix are `confomeval`. See [installation](docs/installation.md) for dependency ranges, validated versions, Windows activation and optional residue-reader training.

Run the minimal representation example or rescore the archived reference predictions manually:

```bash
python examples/representation_archive.py --output-dir runs/example-archive
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Each output directory must be new. Reference rescoring uses the included predictions on CPU and compares all seven arms per scored task, including group-bootstrap intervals, at a numerical tolerance of `1e-10`. Refitting readers additionally requires external representation archives.

## Evaluate a model

1. Prepare exact-sequence archives using the [model input contract](docs/model_input.md) and [minimal example](examples/README.md).
2. Fit the prescribed pooled readers or score complete outer-fold predictions with the [submission and scoring commands](docs/prediction_submission.md).
3. For residue predictions, use the [disorder reference interface](docs/residue_evaluation.md).

Sample identities, label masks, coordinate bounds, prediction completeness and registered group isolation are checked during input loading and scoring. [Local tests](docs/testing.md) cover valid inputs and rejection cases. Tests and reference rescoring are invoked manually.

## Repository layout

| Directory | Contents |
|---|---|
| [`src/confomeval/`](src/confomeval) | Importable model interfaces, reader fitting, scientific input checks and scorers |
| [`configs/`](configs) | [Scoring protocol](configs/scoring_protocol.json), task-specific protocols and model extraction settings |
| [`examples/`](examples) | A minimal representation archive example |
| [`docs/`](docs) | Installation, data sources, model input and evaluation instructions |
| [`data/`](data) | Versioned sequences, labels, groups, folds, source provenance and the [task catalog](data/task_catalog.tsv) |
| [`results/`](results) | Archived predictions, fitted reader parameters, metrics and run provenance; [task dispositions](results/task_dispositions.md) |
| [`tests/`](tests) | Local tests run with `python -m unittest discover -s tests -v` |

## Data, license and citation

[Data sources and licenses](docs/data.md) describes the actual source databases, literature evidence, data versions and source-specific terms. The [MIT License](LICENSE) covers the original software. Recorded source-lineage and redistribution gaps for KnotProt-derived and literature-derived annotations remain documented with the data.

Use [CITATION.cff](CITATION.cff) for the software citation and also cite the primary data sources used by your task. See [CONTRIBUTING.md](CONTRIBUTING.md) for changes and new model results, and the [changelog](docs/changelog.md) for version history.
