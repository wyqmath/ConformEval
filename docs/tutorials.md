# Step-by-step tutorials

Run all commands from the cloned repository root. Use a new output directory for each run.

For an executable tutorial, open the [ESM-2 8M notebook](../examples/notebooks/esm2_8m.ipynb) in [Colab](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/esm2_8m.ipynb). Its default route downloads the complete knot data and runs inference, training and independent rescoring. [Notebook instructions and actual validation](notebook.md) include the small-example switch and output locations.

## Tutorial 1: reproduce reference scores from included predictions

This route needs the base package and the included data/results:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv
source .venv/bin/activate
python -m pip install .
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Expected output: `runs/reference-replay/REPLAY.json` has `status: PASS`. All seven arms for each of the three tasks match the archived metrics and group-bootstrap intervals at tolerance `1e-10`. Reference primary results are listed in [the results summary](../results/README.md). The output directory must be new. Windows activation commands and their validation status are in [installation](installation.md).

For a small, immediate external-prediction example, the actual TSV and metadata are included:

```bash
confomeval-evaluate   --task KNOT-PRESENCE   --data-dir examples/mini/base   --predictions examples/predictions/knot.tsv   --model-metadata examples/predictions/metadata.json   --protocol configs/scoring_protocol.json   --output-dir runs/external-knot-example
```

Expected: `runs/external-knot-example/metrics.json` reports AUROC **1.0** for 20 synthetic samples. This checks the scoring interface. [Residue predictions and metadata templates](prediction_submission.md) provide the equivalent coordinate-level example.

## Tutorial 2: extract representations, train readers and score a model

The following complete example uses the real ESM-2 8M checkpoint and 30 short synthetic sequences. It exercises all three task interfaces using the included small group folds. Start from a fresh clone and a separate environment:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv-model
source .venv-model/bin/activate
python -m pip install ".[esm2]"
confomeval-download --model-config configs/models/esm2_8m.json
confomeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta examples/mini/sequences.fasta   --output-dir runs/esm8-mini-features   --local-files-only

confomeval-fit-readers   --task KNOT-PRESENCE --data-dir examples/mini/base   --feature-manifest runs/esm8-mini-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-mini-knot --example

confomeval-fit-readers   --task FS-L2-grouped-PU-ranking --data-dir examples/mini/base   --feature-manifest runs/esm8-mini-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-mini-fs --example

confomeval-fit-readers   --task DISORDER-STRUCTURE-REFERENCE-CANDIDATE   --data-dir examples/mini/residue   --feature-manifest runs/esm8-mini-features/representation_manifest.json   --config configs/readers/linear.json   --output-dir runs/esm8-mini-residue --example
```

Expected: extraction reports **30 sequences, 320 dimensions**, with layers **3 and 6**. Each training command ends with `status: PASS`. Each pooled task writes 20 out-of-fold rows per arm; the residue task writes **554 labeled rows per arm**. Outputs are `predictions.tsv.gz`, `model_metadata.json`, `metrics.json`, `inner_selection.tsv`, `outer_selections.tsv`, `parameters/` and `run.json` in the requested directories. Synthetic scores have no biological performance interpretation.

The model config chooses layers, pooling and long-sequence policy. The reader config chooses candidate predictors and regularization. Inner validation selects settings; outer groups are used only for the final predictions. For ESM-2 650M use `configs/models/esm2_650m.json`. For ProtT5 install `.[prott5]` and use `configs/models/prott5.json`, with its residue-spaced tokenizer and 1,024-dimensional outputs. Run each model in its own output directories.

For the complete benchmark, use the full task data rather than the synthetic fixture. This is a complete single-task example:

```bash
confomeval-prepare-inputs   --task KNOT-PRESENCE --data-dir data/base   --output-dir runs/knot-inputs
confomeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta runs/knot-inputs/sequences.fasta   --output-dir runs/esm8-knot-features --sequence-only --local-files-only
confomeval-fit-readers   --task KNOT-PRESENCE --data-dir data/base   --feature-manifest runs/esm8-knot-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-knot
```

Expected: **569 unique sequences**, **734 samples**, **283 groups**, and **5 × 3** folds. [Full-task commands](full_evaluation.md) cover PU ranking and generic residue training. The [small NumPy example](../examples/README.md) also exercises training without model downloads.
