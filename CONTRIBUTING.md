# Contributing to ConformEval

Install from the repository root with `python -m pip install -e .`, then run:

```bash
python -m unittest discover -s tests -v
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

The replay output directory must be new. These are local, manually invoked commands.

## Code and examples

Place importable code in `src/confomeval`, protocols in `configs`, minimal examples in `examples` and usage documentation in `docs`. Python comments and docstrings are welcome. Add focused tests for functional changes, including rejection cases for malformed masks, coordinates, prediction identities and split contamination. Run commands from a clone after installing the package.

## Data and results

Scientific data changes require a new version with explicit source identifiers, versions, exclusions and provenance. Preserve historical labels, masks, coordinates, groups and results. Keep unknown biological targets masked, and distinguish PU learning labels from biological targets. Check complete biological binding groups and historical use before claiming independent evaluation.

Data integrity manifests apply to the versioned scientific artifacts they describe. They do not register source files or documentation. Update the affected data or result hashes when intentionally changing a scientific artifact; adding a source comment or example needs no manifest entry.

For a new model result, record the model revision, extractor source, input and code hashes, actual dependency versions, timestamps, warnings, predictions, selected settings and fitted parameters. Fit preprocessing on training data and select settings only through the declared inner validation. Store the actual run configuration and outputs under `results`; keep reusable model settings under `configs`. Review the result before describing it as an accepted benchmark result.

Use the [model input contract](docs/model_input.md) and [scoring instructions](docs/prediction_submission.md). Resolve source-specific permissions before asserting a data license; see [data sources and licenses](docs/data.md).
