# Installation

ConformEval requires Python 3.11 or newer. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

On Windows use `.venv\Scripts\activate` instead of `source`. For development use `python -m pip install -e .`. Commands accept explicit data, configuration and output paths, so they can be invoked outside the clone. The wheel contains the Python code; versioned data and archived results are read from the clone using those paths.

## Dependencies

The base package requires NumPy `>=1.26.4,<3`, SciPy `>=1.11,<2` and scikit-learn `>=1.5,<2`. These ranges are declared in [`pyproject.toml`](../pyproject.toml) and [`requirements.txt`](../requirements.txt). The standard library `unittest` runner is used for local tests.

The original reference computation used Python 3.11.16, NumPy 1.26.4, SciPy 1.13.1, scikit-learn 1.6.1 and PyTorch 2.5.1 with CUDA 12.1. [`requirements-reference.txt`](../requirements-reference.txt) records these numerical package versions. The actual training environment is also preserved in each archived run's provenance under `results`.

A clean local scoring install on macOS arm64 was validated with Python 3.14.5, NumPy 2.5.4, SciPy 1.18.1, scikit-learn 1.9.1, joblib 1.6.0 and threadpoolctl 3.7.0. This local environment is separate from the archived cluster training environment. A new run records its own actual versions.

## Optional residue-reader fitting

```bash
python -m pip install ".[reference-residue-training]"
```

This adds PyTorch `>=2.5,<3`. Install a CUDA or CPU build appropriate for your machine. The base package can score predictions and fit pooled readers without PyTorch. Residue-reader fitting needs the exact external representation cache described in [residue evaluation](residue_evaluation.md).

## Local checks

```bash
python -m unittest discover -s tests -v
python examples/representation_archive.py --output-dir runs/example-archive
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Run these manually. See [testing and reproduction](testing.md) for the checks and output scope.
