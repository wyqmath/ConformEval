# ConformEval: Protein Conformation Evaluation Benchmark

Evaluate protein language model representations on protein knots, fold switching, and residue-level disorder using fixed data splits and shared scorers.

[Download dataset](https://zenodo.org/records/23433394)

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/ESM-2.ipynb)

## Quick start

Open the [ESM-2](examples/notebooks/ESM-2.ipynb) notebook, or run a small scoring example locally with Python 3.9+:

```bash
git clone -c core.autocrlf=false https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m pip install .
conformeval-evaluate --example --output-dir runs/knot-example
```

This writes `runs/knot-example/metrics.json` with AUROC **1.0** on 20 synthetic samples. The [ESM-2](examples/notebooks/ESM-2.ipynb) notebook runs model inference and training; the [reference results](results/README.md) include the score-reproduction command.

## Tasks and results

| Task | Samples | Metric | ESM-2 650M | ESM-2 8M | ProtT5 |
|---|---:|---|---:|---:|---|
| Knot presence | 734 | AUROC ↑ | 0.868741 | 0.814237 | Not tested |
| Fold-switch ranking | 8,879 | Mean positive rank percentile ↓ | 0.476592 | 0.304406 | Not tested |
| Residue disorder | 1,279 records | Mean protein AP ↑ | 0.789856 | 0.738289 | Not tested |

The fold-switch background is unlabeled. Disorder uses experimental disorder labels and ordered/PDB structural references. Results use grouped development folds; model configurations differ. See [results and baselines](results/README.md) and [evaluation details](docs/evaluation.md).

[Supported models](docs/models.md): ESM-2 8M, ESM-2 650M, and ProtT5.

## Dataset

[Download the dataset](https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1) (3.85 MB) · [DOI: 10.5281/zenodo.23433394](https://doi.org/10.5281/zenodo.23433394). The archive includes sequences, labels, and fixed splits for the three tasks. [Data formats and sources](docs/data.md) describe the files and evaluator inputs.

## Repository layout

```text
src/conformeval/  Protein model inference, reader training, scoring
configs/         Model settings and evaluation protocols
examples/        Runnable notebook, small inputs, prediction templates
data/            Task inputs, labels, fixed cross-validation splits
results/         Model and baseline score summaries
docs/            Data formats, model setup, scoring rules
tests/           Automated checks of software and scoring
pyproject.toml   Package dependencies, installation, command names
```

## License

[MIT](LICENSE) covers the software. Dataset sources and their terms are listed in the [data guide](docs/data.md#sources).
