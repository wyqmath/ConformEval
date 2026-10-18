# Model adapter and representation contract

The task data and folds are independent of the encoder. The built-in `esm2` adapter loads ESM-2; `prott5` loads the ProtT5 encoder. Both produce the same archive structure and use the same training/scoring entrypoints.

## Configuration and extraction

A model configuration specifies `adapter`, `model_id`, immutable `revision`, positive 1-based transformer `layers`, `pooling` (`mean` or `max`), `dtype`, `device`, `store_residues` and `long_sequence`. Examples are in [`configs/models`](../configs/models).

```bash
conformeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta examples/mini/sequences.fasta   --output-dir runs/esm8-features
```

The input FASTA may contain repeated sequence strings under distinct sample identifiers. Extraction deduplicates by the SHA-256 of the original, full input sequence. Labels and folds are never passed to the encoder.

## Uniform NPZ archives

Each archive contains:

| Array | Shape | Meaning |
|---|---|---|
| `meta` | scalar JSON string | Original sequence hash/length, model ID/revision, layers, dimension, pooling, normalization, extraction config and windows |
| `sequence_layers` | `[n_layers, dimension]` | Configured pooling over all full-sequence residue vectors |
| `resid_layers` | `[n_layers, n_positions, dimension]` | Selected transformer-layer vectors at exact residue coordinates |
| `resid_indices` | `[n_positions]` | Strictly increasing, unique, one-based coordinates in the original full sequence |

`--sequence-only` omits the last two arrays. `representation_manifest.json` records the model revision, output dimension, candidate layers, pooling, complete extraction configuration and each archive's relative path, sequence identity, length and SHA-256. Relative archive filenames are serialized with `/` on every operating system. Absolute paths, escapes outside the manifest directory, mixed dimensions, model/config mismatches, missing labeled coordinates and non-finite values are rejected by the training loader.

Every archive uses `len_full == len_used` and `truncated: false`. A sequence that cannot be represented is an explicit error. For a fold-switch pair, keep both endpoint identities, even when their sequences are identical; training validates them before averaging the two vectors.

## Tokens and coordinates

The ESM-2 adapter uses its protein tokenizer. The ProtT5 adapter inserts one space between residues and applies the published rare/ambiguous-residue normalization to `X`. Original sequence hashes and coordinate numbers remain unchanged; the normalized sequence hash and substitution rule are recorded.

Special and padding tokens are excluded through tokenizer masks. The token pieces must reconstruct the normalized sequence exactly. A missing, merged or unexpected token cannot silently shift coordinates. The adapter supports mapping multiple residues to a subword token, while the reconstructed residue string remains the identity check.

ProtT5 normalization follows the [encoder model card](https://huggingface.co/Rostlab/prot_t5_xl_half_uniref50-enc). ESM-2 model definitions are documented by the [model authors](https://github.com/facebookresearch/esm).

## Long sequences

`long_sequence.policy` is `window` or `reject`. The supplied configs use windows of at most 1,022 residues with 128-residue overlap. Each window is tokenized independently. Overlapping residue vectors are averaged at the exact original coordinates, then full-sequence pooling is calculated once per residue. No residue is dropped and overlap positions are not counted twice. This uses local window context rather than a single whole-chain forward pass; apply the same declared policy when comparing models.

Set `policy: reject` to fail on sequences longer than the configured maximum. The ESM-2 positional limit is checked before inference. Window parameters and actual half-open window coordinates are stored in each archive.

## Supplying your own representations

[`examples/representation_archive.py`](../examples/representation_archive.py) is a runnable NumPy writer for this contract. It creates complete synthetic archives and a manifest for the small data fixture. Replace its synthetic arrays with your encoder's vectors and preserve identity/coordinate metadata.

[`examples/adapter_usage.py`](../examples/adapter_usage.py) shows the real adapter API: load a config, instantiate the adapter, call `encode(sequence)` and validate output shapes. Both supplied architectures implement this contract.

For a directory of external sequence archives, build a validated manifest with explicit layer IDs:

```bash
conformeval-build-feature-manifest   --task KNOT-PRESENCE --data-dir data/base   --archive-root /path/to/model-archives   --output-manifest /path/to/model-archives/representation_manifest.json   --model-id organization/model-name --model-revision immutable-revision   --pooling mean --layers 3 6   --extractor-code-sha256 0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef
```

Replace the model identity and extractor fingerprint with the actual values. The command verifies sequence coverage and archive identity against the selected task before writing its manifest.
