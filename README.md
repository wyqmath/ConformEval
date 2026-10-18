# ConformEval: Protein Conformation Evaluation Benchmark

Evaluate protein language model representations on protein knots, fold switching,
and residue-level disorder using fixed data splits and shared scorers.

[Data](https://zenodo.org/records/23433394) · [Usage](docs/usage.md) · [Colab](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/esm2_8m.ipynb)

## Quick start

Open the [ESM-2 8M notebook](examples/notebooks/esm2_8m.ipynb), or run a small
scoring example locally with Python 3.11+:

```bash
git clone -c core.autocrlf=false https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m pip install .
conformeval-evaluate --task KNOT-PRESENCE \
  --data-dir examples/mini/base \
  --predictions examples/predictions/knot.tsv \
  --model-metadata examples/predictions/metadata.json \
  --protocol configs/scoring_protocol.json \
  --output-dir runs/knot-example
```

This writes `runs/knot-example/metrics.json` with AUROC **1.0** on 20 synthetic
samples. For real model inference, training, and reference-score reproduction,
follow the [usage guide](docs/usage.md).

## Tasks and results

| Task | Samples | Metric | ESM-2 650M | ESM-2 8M | ProtT5 |
|---|---:|---|---:|---:|---|
| Knot presence | 734 | AUROC ↑ | 0.868741 | 0.814237 | Not tested |
| Fold-switch ranking | 8,879 | Mean positive rank percentile ↓ | 0.476592 | 0.304406 | Not tested |
| Residue disorder | 1,279 records | Mean protein AP ↑ | 0.789856 | 0.738289 | Not tested |

The fold-switch background is unlabeled. Disorder uses experimental disorder
labels and ordered/PDB structural references. Results use grouped development
folds; model configurations differ. See [results and baselines](results/README.md)
and the [evaluation protocol](docs/evaluation.md).

[Supported models](docs/models.md): ESM-2 8M, ESM-2 650M, and ProtT5.

## Data

[Download the dataset](https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1)
(3.85 MB) · [DOI: 10.5281/zenodo.23433394](https://doi.org/10.5281/zenodo.23433394).
The archive includes sequences, labels, and fixed splits for the three tasks.
[Data formats and sources](docs/data.md) describe the files and reference bundle.

## Repository layout

```text
src/conformeval/  Model adapters, training, scoring
configs/         Model settings and evaluation protocols
examples/        Notebook, small datasets, prediction templates
data/            Task and download indexes
results/         Scores and result indexes
docs/            Usage, data, models, evaluation
tests/           Local tests
```

## License

[MIT](LICENSE) covers the software. Dataset sources and their terms are listed
in the [data guide](docs/data.md#sources).
