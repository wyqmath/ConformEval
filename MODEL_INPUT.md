# Model representation input contract

PhysInfoBench stores biological labels and splits independently of any model. A model adapter must create one archive per exact input sequence and call `scripts/physinfobench_model_interface.py` to validate the representation before alignment or scoring.

## Mean representation archive

An `.npz` archive contains:

- `meta`: a JSON object stored as a NumPy scalar string. Required fields are `seq_sha256`, `len_full`, `len_used`, `truncated`, `model`, `revision`, and `layer_set`.
- `mean_layers`: finite numeric matrix with shape `[n_layers, dimension]`.

The `layer_set` field must identify the 1-based layer order, for example `mean_5_11_17_23_29_33`. `load_mean_representation` checks the archive SHA-256, exact sequence hash and length, truncation state, model revision, requested layers, numeric values, and array shape. It returns vectors in the requested layer order as `float32`.

For KNOT-PRESENCE, supply the full chain sequence. For FS-L2, supply each experimental endpoint separately, validate both sequence identities, then average the two endpoint vectors to form the paired case vector. Do not collapse endpoints before validating their identities.

## Residue representation archive

An `.npz` archive contains the same `meta` fields plus:

- `resid_layers`: finite numeric array with shape `[n_layers, n_represented_residues, dimension]`.
- `resid_indices`: integer array listing the exact, strictly increasing 1-based full-sequence coordinate for every represented residue vector.

The metadata identifies the 1-based residue-layer order, for example `mean_5_11_17_23_29_33+resid_11_23_33|domidx`. `load_residue_representation` rejects duplicate, out-of-range, or unordered coordinates. `align_residue_labels` joins by sequence SHA-256 and exact 1-based position. Masked unknown labels are omitted; if any known label lacks a representation coordinate, the alignment fails.

## Build the representation manifest

The nested reader requires one archive for each unique eligible input sequence across the two pooled scored tasks. Its JSON manifest records the model ID and immutable revision, pooling description, extractor source hash, frozen candidate layers `[17, 33]`, and an `archives` list. Each list item contains `sequence_sha256`, full `length`, a manifest-relative `file`, and the archive's `archive_sha256`. Paths may not escape the manifest directory.

After writing the one-sequence `.npz` archives, generate and validate the manifest rather than writing thousands of rows by hand:

```bash
physinfobench-build-feature-manifest \
  --data-dir /path/to/physinfobench/data/v0.3.1 \
  --archive-root /path/to/model_run/features \
  --output-manifest /path/to/model_run/features/representation_manifest.json \
  --model-id organization/model-name \
  --model-revision immutable-checkpoint-revision \
  --pooling mean_over_residues_without_special_tokens \
  --extractor-code-sha256 64-character-extractor-source-hash
```

The command verifies the frozen dataset and folds, exact sequence and endpoint identities, unique eligible sequence coverage, selected layers, lengths, truncation state, finite values, model revision, and feature dimensions. It refuses to overwrite the manifest and rejects conflicting duplicate archives.

The archive metadata must use the same model and revision as the top-level manifest. Supply the exact full input sequence without truncation and save mean vectors for layers 17 and 33 in that declared 1-based order. For a fold-switch pair, provide both endpoint sequences separately; the reader averages their validated vectors only after loading. Do not put labels in the representation manifest or archive.

## Scope of the supplied code

The copied interface module validates cached model representations and label-coordinate joins. It does not download model weights or run model inference. `physinfobench-fit-readers` fits the prescribed nested readers from the manifest and representation archives, writes predictions and fitted parameters for all seven arms, and invokes the official scorer for KNOT-PRESENCE and FS-L2. `physinfobench-evaluate` remains available for complete externally generated outer-fold prediction files. Scores cannot prove which weights or training examples were used, so preserve model and extractor provenance with each run. See [PREDICTION_SUBMISSION.md](PREDICTION_SUBMISSION.md) for commands, metrics, and required metadata.

The interface was independently tested against actual ESM-2 mean and residue caches in the source project. Its tests are included in this delivery and can be run with:

```bash
python -m unittest discover -s tests -v
```

The scored disorder structure-reference task has its own frozen residue folds, label bundle and scorer. See [RESIDUE_EVALUATION.md](RESIDUE_EVALUATION.md). The pooled runner described above is for topology and PU ranking.
