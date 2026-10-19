# ConformEval: Protein Conformation Evaluation Benchmark

Protein function depends on the structures a sequence can adopt and the behavior of regions that remain disordered. Understanding how protein language models encode these properties helps assess their usefulness for studying the physical behavior of proteins. ConformEval (Protein Conformation Evaluation Benchmark) evaluates model representations through three tasks: protein knot detection, fold-switch ranking, and residue-level disorder prediction. It combines curated sequences and reference annotations with tools for extracting representations, training downstream predictors, and computing task-specific scores. The datasets comprise 734 knot samples, 8,879 fold-switch ranking items, and 1,279 disorder records. Adapters for ESM-2 and ProtT5 and runnable examples let researchers evaluate their own models and examine how model layers and downstream predictors affect performance.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23433394.svg)](https://doi.org/10.5281/zenodo.23433394)

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/ESM-2.ipynb)

## Quick start

Open the [ESM-2](examples/notebooks/ESM-2.ipynb) notebook, or run a small scoring example locally with Python 3.9+:

```bash
git clone -c core.autocrlf=false https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m pip install .
conformeval-evaluate --example --output-dir runs/knot-example
```

This writes the example's scores to `runs/knot-example/metrics.json`. The [ESM-2](examples/notebooks/ESM-2.ipynb) notebook runs model inference and training.

## Tasks

| Task | Samples | Metric |
|---|---:|---|
| Knot presence | 734 | AUROC ↑ |
| Fold-switch ranking | 8,879 | Mean positive rank percentile ↓ |
| Residue disorder | 1,279 records | Mean protein AP ↑ |

The fold-switch background is unlabeled. Disorder uses experimental disorder labels and ordered/PDB structural references. See [evaluation details](docs/evaluation.md).

[Supported models](docs/models.md): ESM-2 8M, ESM-2 650M, and ProtT5.

## Dataset

[Download the dataset](https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1) (3.85 MB) · [DOI: 10.5281/zenodo.23433394](https://doi.org/10.5281/zenodo.23433394). The archive includes sequences, labels, and splits for the three tasks. [Data formats and sources](docs/data.md) describe the files and evaluator inputs.

## Repository layout

```text
src/conformeval/  Protein model inference, reader training, scoring
configs/         Model settings and evaluation protocols
examples/        Runnable notebook, small inputs, prediction templates
data/            Task inputs, labels, cross-validation splits
docs/            Data formats, model setup, scoring rules
pyproject.toml   Package dependencies, installation, command names
```

## License

[MIT](LICENSE) covers the software. Dataset sources and their terms are listed in the [data guide](docs/data.md#sources).
