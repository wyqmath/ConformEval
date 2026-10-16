# Model validation

ConformEval 0.5.0 supplies real encoder adapters and one generic nested reader/scoring path. The following tests were actually run on 2026-10-16. The small fixtures test software behavior; the complete-task results use the versioned biological datasets and their saved groups.

| Model | Pinned checkpoint revision | Real extraction | Three small tasks | Complete biological tasks |
|---|---|---|---|---|
| `facebook/esm2_t6_8M_UR50D` | `c731040fcd8d73dceaa04b0a8e6329b345b0f5df` | PASS, layers 3/6, dimension 320 | PASS, knot/ranking/residue | All three newly run with this adapter and generic trainer |
| `facebook/esm2_t33_650M_UR50D` | `08e4846e537177426273712802403f7ba8261b6c` | PASS, layers 17/33, dimension 1,280 | PASS, knot/ranking/residue | Existing archived reference predictions replayed; full inference/refitting with the new adapter not run |
| `Rostlab/prot_t5_xl_half_uniref50-enc` | `94a6abc029ae13029317b140b7424e012bf8dfbf` | PASS, layers 12/24, dimension 1,024 | PASS, knot/ranking/residue | Not run |

All real small extractions cover 30 sequences and 867 residue coordinates. Sequence hashes, archive hashes, full coordinate coverage, model identity and dimensions were checked. ESM encoder loading ignores the checkpoint's unused language-model head; no missing/mismatched encoder parameters or loading errors were reported. ProtT5 reported no missing, mismatched or unexpected parameters. The actual SentencePiece tokenizer was checked with residue-spaced amino acids and EOS exclusion. Model-specific normalization and original-sequence identities are described in [the input contract](model_input.md).

[Machine-readable records](../results/adapter_validation_v1/validation.json) retain extraction manifests and completed small-task run records. [Resource measurements](resources.md) distinguish tested hardware/precision from untested devices. Base-only installed-package tests passed on Python 3.11 and 3.14. GPU inference and full ProtT5/650M runs through the new adapter remain untested.

## Complete ESM-2 8M results

These are frozen grouped-development scores with the same samples, folds and scoring code used by the archived reference. Candidate layers 3/6 and C values 0.1/1.0 were selected only within inner validation. Each outer fold fits its own training-only scaler and reader. The configured baseline arms use the same held groups.

| Task | Nested-reader primary result | Saved output |
|---|---:|---|
| Knot presence | AUROC **0.814237** | [metrics](../results/esm2_8m_v1/knot/metrics.json), [run](../results/esm2_8m_v1/knot/run.json) |
| Fold-switch PU ranking | Mean held-case rank percentile **0.304406**, lower is better | [metrics](../results/esm2_8m_v1/fs/metrics.json), [run](../results/esm2_8m_v1/fs/run.json) |
| Disorder reference | Mean protein AP **0.738289** | [metrics](../results/esm2_8m_v1/residue/metrics.json), [run](../results/esm2_8m_v1/residue/run.json) |

Each directory includes complete out-of-fold predictions, fitted coefficients and preprocessing parameters, inner candidate scores, outer selections, model metadata, extraction provenance and actual environment/timestamps. The residue directory also contains per-protein metrics. These results do not establish independent prospective generalization. Layers and search grids differ from the older archived 650M experiment; this is a functioning smaller-model reference, not a controlled model-size ablation.

## Provenance and export

Scientific prediction/parameter/table files are copied byte-for-byte from completed runs. JSON display spelling was normalized to ConformEval and absolute machine paths were replaced by portable source labels. Input hashes and executed-code hashes remain unchanged. `source_run_sha256` and `source_output_sha256` preserve the original run/output digests; `output_sha256` covers exported outputs. `exported_package_code_sha256` identifies final source after the display-string correction. Model weights and representation caches are external inputs; their complete hashes remain recorded.

The local unit suite covers dimensions/depths below layer 17, special tokens, multi-residue tokenizer pieces, window stitching, mean/max pooling, coordinate/mask/identity failures, incompatible extraction metadata, split contamination, train-only scaling, validation-only selection and anchored numerical ties. The historical three-task seven-arm reference replay still matches archived metrics and confidence intervals at absolute tolerance `1e-10`.
