# ConformEval: Protein Conformation Evaluation Benchmark

Protein function depends on the structures a sequence can adopt and the behavior of regions that remain disordered. Understanding how protein language models encode these properties helps assess their usefulness for studying the physical behavior of proteins. ConformEval (Protein Conformation Evaluation Benchmark) evaluates model representations through three tasks: protein knot detection, fold-switch ranking, and residue-level disorder prediction. It combines curated sequences and reference annotations with tools for extracting representations, training downstream predictors, and computing task-specific scores. The datasets comprise 734 knot samples, 8,879 fold-switch ranking items, and 1,279 disorder records. Adapters for [ESM-2](https://doi.org/10.1126/science.ade2574) and [ProtT5](https://doi.org/10.1109/TPAMI.2021.3095381) and runnable examples let researchers evaluate their own models and examine how model layers and downstream predictors affect performance.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23433394.svg)](https://doi.org/10.5281/zenodo.23433394)

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/ConformEval_tutorial.ipynb)

## Quick start

Use Python 3.9 or later to run the scoring example:

```bash
git clone -c core.autocrlf=false https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m pip install .
conformeval-evaluate --example --output-dir runs/knot-example
```

The score is written to `runs/knot-example/metrics.json`. The [ConformEval tutorial](examples/notebooks/ConformEval_tutorial.ipynb) covers data download, representation extraction, predictor training, and scoring, using [ESM-2 8M](https://huggingface.co/facebook/esm2_t6_8M_UR50D) as an example model.

## Tasks

| Task | Input | Samples | Metric |
|---|---|---:|---|
| Knot presence | Protein sequence | 734 | AUROC ↑ |
| Fold-switch ranking | Sequences of confirmed cases and background items | 8,879 | Mean positive rank percentile ↓ |
| Residue disorder | Protein sequence and residue labels | 1,279 records | Mean protein AP ↑ |

Knot presence includes 112 knotted and 622 unknotted samples. Fold-switch ranking includes 10 confirmed cases and 8,869 unlabeled background items. Disorder includes 391,718 labeled residues: 73,766 experimental disorder labels, 35 direct ordered labels, and 317,917 positions observed in experimental PDB structures. Unknown positions remain masked.

## Models

| Model | Configuration | Dimensions | Small three-task example | Full benchmark |
|---|---|---:|---|---|
| [ESM-2 8M](https://huggingface.co/facebook/esm2_t6_8M_UR50D) | [esm2_8m.json](configs/models/esm2_8m.json) | 320 | Tested | All three tasks |
| [ESM-2 650M](https://huggingface.co/facebook/esm2_t33_650M_UR50D) | [esm2_650m.json](configs/models/esm2_650m.json) | 1,280 | Tested | Refitting with the new adapter: Not tested |
| [ProtT5 XL encoder](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc) | [prott5.json](configs/models/prott5.json) | 1,024 | Tested | Not tested |

Install the inference dependencies for your model:

```bash
python -m pip install '.[esm2]'
# For ProtT5:
python -m pip install '.[prott5]'
```

Model configurations select the checkpoint revision, layers, pooling, precision, and sequence window settings. Inference uses CPU float32 by default. The default window contains up to 1,022 residues with a 128-residue overlap. Overlapping residue vectors are averaged at their sequence coordinates before pooling; special and padding tokens are excluded. Set `long_sequence.policy` to `reject` to reject sequences longer than the configured window.

For another encoder, follow [adapter_usage.py](examples/adapter_usage.py) for adapter calls or [representation_archive.py](examples/representation_archive.py) for an external feature archive. Each archive contains one NPZ file per sequence and a JSON index recording sequence IDs, lengths, model revision, layers, dimensions, pooling, extraction settings, and filenames. Residue vectors use one-based sequence coordinates. Predictor candidates and regularization values are configured in [linear.json](configs/readers/linear.json).

## Data

[Download ConformEval-data.zip](https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1) (3.85 MB) from [Zenodo](https://zenodo.org/records/23433394).

```bash
curl --fail --location 'https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1' --output ConformEval-data.zip
unzip -q ConformEval-data.zip -d datasets
```

The archive extracts into `datasets/ConformEval-data/` with `knot_presence/`, `fold_switching/`, and `disorder/` directories. Each task includes `sequences.fasta`, `samples.tsv`, `test_splits.tsv`, and `validation_splits.tsv`; disorder also includes `labels.tsv`. The source register is `sources.tsv`.

Evaluator inputs in this repository use the following layout:

| Path | Contents |
|---|---|
| `data/base/` | Sequences, task labels, source records, and grouped folds |
| `data/disorder_reference/` | Residue labels and grouped folds |
| `data/fs_region_controls/` | Paired-structure annotations and supporting evidence |

Tables use tab-separated UTF-8 text. Read `.tsv.gz` files with `gzip.open(path, "rt", encoding="utf-8", newline="")` or `pandas.read_csv(path, sep="\t")`. Residue coordinates are one-based positions in the declared sequence. Task manifests list the required files. Case records use `source_ref` to identify entries in [case_sources.tsv](data/base/case_sources.tsv), which records source filenames and available DOI, PMCID, or source notes.

## Evaluation

Knot detection uses AUROC over chain-level predictions. Fold-switch ranking measures each held-out confirmed case's rank percentile within its test pool; the background is unlabeled and has no verified negative labels. Disorder scores known reference positions, calculates average precision for each protein record, and averages across records.

All models use the same group assignments and scoring implementation. Knot and disorder tasks have five outer folds; fold-switch ranking has ten. Each outer training set has three inner validation folds for choosing model layers and predictor settings. Preprocessing is fitted on training data. Predictions are scored on outer test groups. The scoring commands accept TSV prediction files; [predictions/](examples/predictions/) contains example files and a model metadata template.

The scorers check sample identities, prediction coverage, folds, coordinates, and label masks. Confidence intervals resample groups from the saved predictions and exclude uncertainty from model refitting or dataset selection. The current splits support retrospective grouped-development evaluation. Complete biological-family independence has not been established; the source group audit has 37 unresolved PDB-key gaps. The [task catalog](data/task_catalog.tsv) lists additional tasks and their data requirements.

## Repository layout

```text
README.md        Project description, installation, data, models, scoring
src/             Protein model inference, predictor training, scoring
configs/         Model settings and evaluation protocols
examples/        Notebook, example inputs, prediction templates
data/            Task inputs, labels, cross-validation splits
pyproject.toml   Package dependencies, installation, command names
```

The examples contain 30 synthetic sequences, 20 knot samples, 20 ranking items, and 554 labeled residue positions. Their scores demonstrate the software interface. The `representation_archive.py` example generates deterministic NumPy features for testing external model integration.

## Sources and license

Sequences and annotations come from UniProt, DisProt, RCSB PDB/SIFTS, KnotProt-derived legacy annotations, and primary publications. [sources.tsv](data/base/sources.tsv) records source versions, attribution requirements, and distribution terms; [case_sources.tsv](data/base/case_sources.tsv) records publication evidence. KnotProt redistribution terms and the original scored-label source snapshot remain unresolved. Literature-derived records retain their source-specific terms.

[MIT](LICENSE) covers the software. Dataset reuse follows the terms recorded for each source.
