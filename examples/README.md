# Examples

- [ESM-2 8M notebook](notebooks/esm2_8m.ipynb): full knot task or a synthetic example.
- `mini/`: 30 synthetic sequences, with 20 knot samples, 20 ranking items and 554 labeled residue positions.
- `predictions/`: external prediction files and metadata templates.
- `adapter_usage.py`: direct model-adapter calls.
- `representation_archive.py`: deterministic NumPy representations without model downloads.

Run the [ESM-2 8M notebook](notebooks/esm2_8m.ipynb) for model inference and training, or follow the [model guide](../docs/models.md) to connect another encoder. Synthetic scores are interface examples and are excluded from biological results.
