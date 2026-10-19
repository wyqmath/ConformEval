# Models

| Model | Config | Dimensions | Small three-task example | Complete tasks |
|---|---|---:|---|---|
| ESM-2 8M | `configs/models/esm2_8m.json` | 320 | Tested | All three |
| ESM-2 650M | `configs/models/esm2_650m.json` | 1,280 | Tested | Complete refitting with the new adapter: Not tested |
| ProtT5 XL encoder | `configs/models/prott5.json` | 1,024 | Tested | Not tested |

Install `.[esm2]` for ESM-2 or `.[prott5]` for ProtT5. Base training and scoring
require only numerical dependencies. Model configs pin checkpoint revisions
and select layers, pooling, precision and long-sequence policy.

## Inference settings

CPU float32 is the default. `conformeval-extract` accepts `--device`, `--dtype`,
`--cache-dir`, `--snapshot`, and `--local-files-only`. Weights are downloaded
outside the source repository. Keep device, precision, windows and pooling
consistent when comparing models.

The default window is 1,022 residues with 128-residue overlap. Overlapping
vectors are averaged at their original coordinates before sequence pooling.
Special/padding tokens are excluded. `long_sequence.policy="reject"` produces
an explicit error for overlength inputs. The implementation preserves all
residue positions; this windowing policy limits their available context.

## Connect another model

Follow [`examples/adapter_usage.py`](../examples/adapter_usage.py) for the built-in
API or [`examples/representation_archive.py`](../examples/representation_archive.py)
for an external representation archive. Keep model loading and tokenization
inside the adapter; use the shared task loaders, reader trainer and scorers.

An archive contains one NPZ file per sequence and a JSON manifest with the model ID, revision, sequence IDs, layers, dimensions, pooling, extraction configuration, and array filenames. Each sequence provides pooled layer vectors and, when requested, a residue matrix whose rows map to one-based sequence coordinates. Use the array keys and manifest schema shown in the example.

Reader config selects predictor candidates and regularization. Inner validation
selects settings; outer held-out scores do not select a model. Layers are
model-specific configuration values, with no fixed global layer numbers.
