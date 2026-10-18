# Runnable small examples

The [ESM-2 8M notebook](notebooks/esm2_8m.ipynb) runs a complete knot benchmark by default, with a switch for the small synthetic knot example. It includes data/checkpoint download, inference, training and independent scoring. [Open in Colab](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/esm2_8m.ipynb) or see [validation and local instructions](../docs/notebook.md).

The `mini` directory contains 30 synthetic protein sequences and complete small group-fold tables for all three task interfaces. The pooled tasks have 20 samples each, and the residue task has 554 known positions. The data manifests mark them `software_example`.

## No model download

```bash
python -m pip install .
python examples/representation_archive.py --output-dir runs/numpy-example-features
conformeval-fit-readers --task KNOT-PRESENCE --data-dir examples/mini/base --feature-manifest runs/numpy-example-features/representation_manifest.json --config configs/readers/linear.json --protocol configs/scoring_protocol.json --output-dir runs/numpy-example-knot --example
```

Expected: `status: PASS`, 20 scored samples per arm, and output predictions/metrics/parameters. The vectors are deterministic synthetic values, not a model-performance result.

## Real models

The [tutorials](../docs/tutorials.md) supply the complete ESM-2 8M extraction/training/scoring walkthrough for all three tasks. To directly exercise either adapter:

```bash
python examples/adapter_usage.py --model-config configs/models/esm2_8m.json
```

Install the model's optional dependency extra first. For ESM-2 650M or ProtT5 choose its supplied model config. The script downloads/loads the pinned checkpoint, represents a short sequence and prints layer/residue/dimension shapes. [Validation](../docs/model_validation.md) records the checkpoints and commands actually tested.

## External predictions

`predictions/knot.tsv` and `predictions/residues.tsv` are complete synthetic prediction files. `predictions/metadata.json` makes them runnable; `metadata.template.json` shows the values to replace for an actual submission. [Scoring instructions](../docs/prediction_submission.md) include copyable commands and expected scores.
