# ESM-2 8M notebook

[View on GitHub](https://github.com/wyqmath/ConformEval/blob/main/examples/notebooks/esm2_8m.ipynb)
or [open in Colab](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/esm2_8m.ipynb).
Run the cells in order, or choose **Runtime → Run all** in Colab. The notebook
uses CPU float32 and requires no Google Drive mount, access token or GPU.

For local JupyterLab use Python >=3.11:

```bash
python -m pip install "jupyterlab>=4,<5"
jupyter lab examples/notebooks/esm2_8m.ipynb
```

The notebook installs the base package and ESM-2 extra with compatible ranges
in its kernel environment. It fetches code/config/data from an immutable
ConformEval commit and verifies every downloaded Git blob. Model download uses
the pinned ESM-2 checkpoint revision. A complete run downloads about 20 MB of
benchmark/code inputs; the small-model checkpoint is about 30 MB. Existing
verified model caches are reused.

## Choose the data size

The default `FULL_BENCHMARK = True` evaluates the complete knot-presence task:
569 encoder sequences, 734 samples and 283 groups. Set it to `False` for the
20-sequence synthetic knot interface example. The mini bundle also includes
other task fixtures, but this notebook prepares only its 20 knot sequences.
Synthetic results are excluded from the biological model result table.

## Outputs and checks

Each execution creates a new `conformeval-notebook-*` workspace. Its `runs/`
directory contains downloaded-file hashes, prepared inputs, representation
archives, reader fits, out-of-fold predictions, a four-column submission TSV,
independent scorer output and logs. The final cell prints its location and
writes `runs/notebook_verification.json` with versions, sample counts, scores
and output hashes. Save required outputs before ending a temporary runtime.

Layer and C selection use inner validation; the scorer uses the fixed outer
predictions. The notebook asserts complete prediction coverage, five outer
folds, 320-dimensional features, and equality between independently rescored
AUROC/intervals and the reader output. It keeps full and mini evaluation
identities explicit.

## Actual validation on 2026-10-18

Both notebook modes were executed with nbclient in fresh local workspaces on
macOS arm64, Python 3.14.6, CPU float32. Data were downloaded from GitHub in
both runs. Full mode reused the existing model cache; mini mode downloaded the
pinned checkpoint into an empty cache. Exact dependencies and hashes are in
the linked records.

| Mode | Unique sequences | Scored samples | Result | Independent rescoring | Record |
|---|---:|---:|---|---|---|
| Complete knot benchmark | 569 | 734 | AUROC 0.814237 | PASS | [Full verification](../results/notebook_validation/full.json) |
| Synthetic interface example | 20 | 20 | Software example only | PASS | [Mini verification](../results/notebook_validation/mini.json) |

The full run took 49.66 seconds with cached weights, including data download,
installation, inference, training and rescoring. Mini mode took 43.58 seconds
with a fresh checkpoint cache. These are local measurements; hosted Colab
execution is **Not tested**. Runtime availability and package resolution may
differ in Colab. The source notebook is committed without machine-specific
outputs; verification records preserve the actual local execution evidence.

[Publication link checks](../results/notebook_validation/links.json) record the
GitHub file page, raw notebook and Colab URL separately from model execution.
Hosted Colab execution remains **Not tested**.

[Original validation records](../results/notebook_validation/archive/README.md)
retain the earlier source and notebook identities unchanged.
