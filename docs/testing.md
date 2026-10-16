# Local tests and reference reproduction

Install the package from the clone before testing:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

The suite checks sequence hashes, archive identities and layers, label masks, coordinate bounds, complete predictions, registered group isolation, train-only preprocessing and failed-run records. Scientific validation functions live in `confomeval.data_validation`, `confomeval.model_interface` and the scoring modules. They are called through ordinary Python imports, input loading and local tests.

## Rescore archived predictions

```bash
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

The output directory must be new. The command rescores all seven arms for each of the three scored tasks: 67,291 pooled prediction rows and 2,742,026 residue prediction rows. It compares pooled point metrics, binding-group intervals and all FS case ranks; the residue comparison covers the full metrics and 2,000 group-bootstrap intervals. The tolerance is `1e-10`. Results and timestamps are saved to the requested output directory, including `REPLAY.json` and per-task scoring records.

Rescoring uses archived predictions. Reader refitting and encoder inference require external representations or model weights; see [model input](model_input.md) and [residue evaluation](residue_evaluation.md).

## Historical provenance

Archived `run.json`, `provenance.json` and related input/output hashes describe the original scientific runs. Their source-code fingerprints and source-project paths are historical identifiers, rather than paths to the reorganized package. [`data/provenance/upstream_artifacts.json`](../data/provenance/upstream_artifacts.json) preserves the upstream data/result fingerprints. Data and result manifests verify only their scientific artifacts. Documentation, source comments and new examples are editable without a repository inventory.
