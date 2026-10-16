# ConformEval

ConformEval evaluates **protein language model representations** on three scoring tasks. It provides fixed sequences and labels, biological binding-group folds, model adapters, nested reader training, scoring tools and reference predictions. Models use the same task samples, group assignments and scoring implementation.

| Scoring task | Input and target | Fixed evaluation | ESM-2 650M reference |
|---|---|---|---|
| `KNOT-PRESENCE` | Full sequence; whole-chain knot presence | 734 samples, 283 groups; 5 outer × 3 inner folds | Nested AUROC **0.868741**, AP 0.687211 |
| `FS-L2-grouped-PU-ranking` | Each positive case's two endpoint sequences, or one background sequence; held-case ranking | 8,879 ranking items, including 10 positives; 10 outer × 3 inner folds | Mean held-case rank percentile **0.476592**, lower is better |
| `DISORDER-STRUCTURE-REFERENCE-CANDIDATE` | Full sequence with residue coordinates; direct disorder versus an ordered structural reference | 1,279 protein records, 1,029 groups, 391,718 labeled rows; 5 outer × 3 inner folds | Layer 33 mean protein AP **0.789856** |

Unknown biological labels remain masked. FS background items have unknown biological targets. The disorder reference task includes experimentally observed PDB positions as reference negatives. The evaluation is retrospective grouped development; task definitions and data sources are in [data documentation](docs/data.md) and the [task catalog](data/task_catalog.tsv).

## Supported models and inputs

| Model | Adapter | Transformer layers / dimension | Tested use |
|---|---|---|---|
| ESM-2 650M | `esm2`, [`esm2_650m.json`](configs/models/esm2_650m.json) | 33 / 1,280 | Archived reference scores for all three tasks; real-weight extraction and all three small task examples |
| ESM-2 8M | `esm2`, [`esm2_8m.json`](configs/models/esm2_8m.json) | 6 / 320 | Real-weight extraction with layers 3/6; all three small examples and all three complete tasks passed; [results](docs/model_validation.md) |
| ProtT5 encoder | `prott5`, [`prott5.json`](configs/models/prott5.json) | 24 / 1,024 | Real-weight extraction and all three small tasks passed; complete biological tasks not run; [details](docs/model_validation.md) |

Inputs are an exact-sequence FASTA and a pinned model configuration. Extraction writes sequence/residue NPZ archives and `representation_manifest.json`. Training reads those archives, task labels and fixed folds, and writes out-of-fold predictions, fitted parameters, selected settings, metrics and `run.json`. Complete external prediction TSVs can also be scored directly. See [the archive and adapter contract](docs/model_input.md).

Base scoring and linear/ridge reader training require Python **3.11+**, NumPy, SciPy and scikit-learn. Model dependencies are optional: install `.[esm2]` for ESM-2 or `.[prott5]` for ProtT5. See [installation](docs/installation.md) for tested versions and [resources](docs/resources.md) for measured time, peak memory, cache sizes and long-sequence handling.

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

Expected output: `runs/reference-replay/REPLAY.json` has `status: PASS`. All seven arms for each of the three tasks match the archived metrics and group-bootstrap intervals at tolerance `1e-10`. The table above shows the reference primary results. The output directory must be new. On Windows activate with `.venv\Scripts\activate`.

For a small, immediate external-prediction example, the actual TSV and metadata are included:

```bash
confomeval-evaluate   --task KNOT-PRESENCE   --data-dir examples/mini/base   --predictions examples/predictions/knot.tsv   --model-metadata examples/predictions/metadata.json   --protocol configs/scoring_protocol.json   --output-dir runs/external-knot-example
```

Expected: `runs/external-knot-example/metrics.json` reports AUROC **1.0** for 20 synthetic samples. This checks the scoring interface. [Residue predictions and metadata templates](docs/prediction_submission.md) provide the equivalent coordinate-level example.

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
confomeval-prepare-inputs   --task KNOT-PRESENCE --data-dir data/v0.3.1   --output-dir runs/knot-inputs
confomeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta runs/knot-inputs/sequences.fasta   --output-dir runs/esm8-knot-features --sequence-only --local-files-only
confomeval-fit-readers   --task KNOT-PRESENCE --data-dir data/v0.3.1   --feature-manifest runs/esm8-knot-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-knot
```

Expected: **569 unique sequences**, **734 samples**, **283 groups**, and **5 × 3** folds. [Full-task commands](docs/full_evaluation.md) cover PU ranking and generic residue training. The [small NumPy example](examples/README.md) also exercises training without model downloads.

## Repository and local checks

`src/confomeval` contains importable code; `configs` holds protocols/model settings; `data` holds versioned task data; `results` contains actual reference results; `examples` contains runnable examples; `docs` contains usage/source information. All checks run locally:

```bash
python -m unittest discover -s tests -v
```

Use [CITATION.cff](CITATION.cff) for the software citation, and cite the primary data sources for scientific use. [LICENSE](LICENSE) covers original software; [data documentation](docs/data.md) records actual source terms. See [CONTRIBUTING.md](CONTRIBUTING.md) for adding models, data or results.
