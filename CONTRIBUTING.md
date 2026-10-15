# Contributing to PhysInfoBench

## Current scope

Version 0.4.0 retains the base v0.3.1 scientific lock and adds the separately accepted disorder reference and fold-switch region reference-label versions. Frozen labels, masks, sequences, binding groups, folds and predictions must not change in place. Scientific changes require a new version and independent review.

Source-specific redistribution review remains open, as recorded in
[`SOURCES_AND_DISTRIBUTION.md`](SOURCES_AND_DISTRIBUTION.md). Publication of this
maintainer-authorized snapshot does not resolve those terms or grant a blanket
data license. A future source-cleared export requires documented terms or
source exclusions, a new manifest and review.

## Before opening a pull request

Run the package checks from this directory:

```bash
python -m pip install ".[test]"
python -m unittest discover -s tests -v
physinfobench-verify --benchmark-root .
physinfobench-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Changes to data or evaluation code must include focused tests for both accepted
inputs and rejection cases. Preserve unknown labels as masked unknowns, validate
sequence and coordinate identities, keep preprocessing inside the training
fold, and verify that test groups never enter selection. Update the data card,
task catalog, release notes, and manifests when the declared release contents
change.

## Fitting a reader from model representations

Use `MODEL_INPUT.md` to provide exact-sequence archives and a manifest, then
run `physinfobench-fit-readers` as documented in `PREDICTION_SUBMISSION.md`.
The runner uses the frozen layers, C grid, inner groups, outer groups, and
metrics. It does not run the encoder, so model inference remains a separate,
traceable step. Preserve the inference revision, extractor source, archive
hashes, and run outputs together.

## Adding a model result

Use an explicit, versioned adapter and configuration. Record the model and
revision, feature archive hashes, software versions, frozen data and split
hashes, timestamps, warnings, predictions, selections, parameters, and output
hashes. Results must cover the same qualified samples, groups, folds, and task
labels as the reference. Keep PU learning labels separate from biological
targets and masks. New-model runs require the established review path before
their results are described as benchmark results or distributed.
