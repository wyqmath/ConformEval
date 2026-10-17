# ConformEval: Protein Conformation Evaluation Benchmark

ConformEval tests how well **protein language model representations** capture protein knots, fold switching, and residue-level disorder. It provides fixed sequences, labels, group folds, model adapters, training tools, and scoring commands. You can reproduce reference scores, evaluate a supported model, or supply your own representations and predictions.

## Tasks and reference results

| Task | Metric | ESM-2 650M | ESM-2 8M | ProtT5 |
| --- | --- | ---: | ---: | --- |
| Knot presence | AUROC ↑ | 0.868741 | 0.814237 | Not tested |
| Fold-switch ranking | Mean held-case rank percentile ↓ | 0.476592 | 0.304406 | Not tested |
| Disorder structure reference | Mean protein AP ↑ | 0.789856 | 0.738289 | Not tested |

The fold-switch background is unlabeled. The disorder task combines experimental disorder labels with ordered annotations and PDB-observed structural references. These are grouped-development results; model configurations differ. [Full results and baselines](results/README.md) document each configuration. The [task catalog](data/task_catalog.tsv) covers all 13 registered tasks, including the ten without model scores. The [data and table guide](docs/data_tables.md) indexes every public TSV/TSV.GZ file, its fields, and its reading method.

Adapters support **ESM-2 8M, ESM-2 650M, and ProtT5**. All three passed real-checkpoint extraction and small-example tests. The 650M scores above come from archived reference runs; full ProtT5 evaluation is not tested. See [model validation](docs/model_validation.md).

## Quick start

Use **Python 3.11+**. From a terminal on macOS or Linux:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

Windows is not tested. See [installation](docs/installation.md) for platform commands and optional dependencies.

### 1. Reproduce reference scores

No model download is required:

```bash
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Expected: `runs/reference-replay/REPLAY.json` reports `status: PASS`. Seven archived arms for each of the three tasks match their reference metrics and group-bootstrap intervals within `1e-10`. Use a new output directory for each run.

### 2. Run a model: ESM-2 8M

This small example downloads a pinned checkpoint, extracts representations from 30 short sequences, then trains and scores a knot-presence reader:

```bash
python -m pip install ".[esm2]"
confomeval-download --model-config configs/models/esm2_8m.json
confomeval-extract \
  --model-config configs/models/esm2_8m.json \
  --input-fasta examples/mini/sequences.fasta \
  --output-dir runs/esm8-mini-features --local-files-only
confomeval-fit-readers \
  --task KNOT-PRESENCE --data-dir examples/mini/base \
  --feature-manifest runs/esm8-mini-features/representation_manifest.json \
  --config configs/readers/linear.json --protocol configs/scoring_protocol.json \
  --output-dir runs/esm8-mini-knot --example
```

Expected: extraction produces **30 sequence archives with 320-dimensional representations**. Training reports `status: PASS` and writes `predictions.tsv.gz` (20 samples per arm), `metrics.json`, fitted parameters, and `run.json` under `runs/esm8-mini-knot/`. These synthetic samples test the interface; their scores do not measure biological performance.

[Detailed tutorials](docs/tutorials.md) cover all three small tasks and external predictions. [Full evaluation](docs/full_evaluation.md) provides commands for the biological datasets. To add another model, follow the [representation contract and adapter example](docs/model_input.md).

## Repository layout

```text
ConformEval/
├── src/confomeval/   # Model adapters, training, and scoring
├── configs/          # Model settings and evaluation protocols
├── data/             # Versioned task data and source metadata
├── results/          # Reference scores and predictions
├── examples/         # Small runnable examples and prediction templates
├── docs/             # Tutorials, data sources, and technical details
├── tests/            # Local tests
├── README.md
├── LICENSE
└── pyproject.toml
```

## Documentation and checks

- [Data definitions, sources, and reuse terms](docs/data.md)
- [Prediction format and scoring your own outputs](docs/prediction_submission.md)
- [Runtime, memory, and long-sequence handling](docs/resources.md)
- [Operating-system validation](docs/platforms.md)
- [Contributing](CONTRIBUTING.md)

Run local tests with `python -m unittest discover -s tests -v`. [LICENSE](LICENSE) covers original software; cite the primary data sources when using the benchmark.
