# Measured resources

These measurements were made on 2026-10-16 on an Apple M5, 10 CPU cores and 16 GiB RAM, macOS 26.4 arm64. Inference used CPU float32 and four PyTorch threads. Times include model loading and output writing, but exclude checkpoint downloads. Peak RSS is process high-water memory, not an incremental allocation. Processes shared this workstation; these are observed runs, not hardware-independent guarantees.

## Environments

| Use | Actual versions |
|---|---|
| Encoder extraction | Python 3.14.6; NumPy 2.5.3; SciPy 1.18.1; scikit-learn 1.9.1; torch 2.14.0; transformers 5.17.0; huggingface-hub 1.32.0; SentencePiece 0.2.2; protobuf 6.33.6 |
| Base training and scoring | Python 3.14.5; NumPy 2.5.4; SciPy 1.18.1; scikit-learn 1.9.1; no inference packages |
| Second base installation | Python 3.11.15; NumPy 2.4.6; SciPy 1.17.1; scikit-learn 1.9.1; no inference packages |

Both base installations passed 50 local tests. Python 3.11 encoder inference, CUDA, MPS, float16 and bfloat16 were not tested in this validation. Each saved run contains its actual environment. The older archived 650M reference environment is recorded in its original results and [historical requirements](reproduction/requirements-reference.txt); it is separate from the new adapter validation. [Reproduction records](reproduction/README.md) index the exact environments; [installation](installation.md) lists the compatible dependency ranges.

## Real-checkpoint small examples

All models processed the same 30 sequences, 867 total residues, lengths 24–35. The three reader runs use 20 knot samples, 20 ranking items and 554 known residue rows per arm. Their labels are synthetic software fixtures.

| Model | Extraction seconds | Peak RSS MiB | Representation archives MiB | Knot / ranking / residue reader seconds |
|---|---:|---:|---:|---:|
| ESM-2 8M | 2.24 | 482.6 | 2.08 | 6.15 / 0.54 / 0.96 |
| ESM-2 650M | 19.06 | 2,624.4 | 8.18 | 6.65 / 0.43 / 3.61 |
| ProtT5 XL encoder | 27.06 | 4,760.5 | 6.57 | 5.98 / 0.36 / 1.85 |

The ProtT5 checkpoint is 2,416,373,051 bytes before loading. Downloads and extraction caches live outside the public repository. Download time depends on the network and is excluded from this table. [Validation records](../results/adapter_validation/validation.json) contain input/checkpoint hashes, immutable revisions, dimensions, archive hashes, actual times and loading diagnostics.

## Complete ESM-2 8M benchmark

Sequence-only archives were used for knot and ranking. Residue archives retain both selected layers for every coordinate. Extraction caches include manifests and provenance but exclude model weights.

| Task | Unique sequences | Extraction seconds / peak RSS MiB | Cache MiB | Training + scoring seconds / peak RSS MiB |
|---|---:|---:|---:|---:|
| Knot presence | 569 | 11.27 / 570.7 | 2.13 | 12.97 / 210.1 |
| Fold-switch PU ranking | 8,888 | 215.48 / 667.3 | 33.29 | 87.36 / 376.8 |
| Disorder reference | 1,275 | 85.36 / 640.2 | 1,244.27 | 828.75 / 4,708.7 |

The full residue run fits the configured inner/outer readers and scores 391,718 rows per arm, including 2,000 group bootstrap replicates. All three complete 8M runs together took about 20 minutes of measured process time, excluding downloads and preparation. Full runs of the new 650M and ProtT5 adapters were not measured; do not extrapolate these small-example times linearly to long proteins or full datasets.

## Long sequences and memory

The default configuration uses windows of at most 1,022 residues with 128-residue overlap. Each window is independently encoded. Overlapping residue vectors are averaged at their original coordinates, then mean or max pooling is applied once across the reconstructed sequence. Special and padding tokens are excluded. No residues are silently truncated. This policy changes the available sequence context and is recorded in each archive.

An actual ESM-2 8M run on a 1,025-residue synthetic sequence retained all 1,025 coordinates: 2.76 seconds and 557.3 MiB peak RSS. A separate real-weight run used layers 1/6, max pooling, 16-residue windows and four-residue overlap on all 30 example sequences. These complement the token projection and overlap unit tests. ProtT5/650M long-sequence inference was not measured.

Set `long_sequence.policy` to `reject` for an explicit overlength error. Adjust `long_sequence.max_residues` and `long_sequence.overlap` in the model JSON before extraction to reduce memory. Long-sequence settings and dtype must be reported with model comparisons. Residue archives and flat training arrays can dominate disk and host memory even if inference windows fit.
