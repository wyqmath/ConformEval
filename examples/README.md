# Minimal representation example

After `python -m pip install .`, run from the repository root:

```bash
python examples/representation_archive.py --output-dir runs/example-archive
```

The script creates a synthetic six-residue sequence archive with layers 17 and 33, then reads it through `confomeval.model_interface`. The output is a validated `2 × 4` matrix. The output directory must be new. The vectors illustrate the schema; use actual encoder representations to evaluate a model.

See the [model input contract](../docs/model_input.md) for residue archives, sequence-coordinate alignment and complete cache manifests, and [prediction submission](../docs/prediction_submission.md) for evaluation commands.
