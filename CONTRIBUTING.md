# Contributing to ConformEval

Install with `python -m pip install -e .` and run local tests:

```bash
python -m unittest discover -s tests -v
conformeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Use new output directories. Code is in `src/conformeval`, model/protocol settings in `configs`, examples in `examples`, documentation in `docs`, versioned scientific inputs in `data`, and actual runs in `results`. Comments and docstrings are welcome.

## Adding a model adapter

Keep model loading/tokenization inside the adapter. Implement full-sequence/residue projection with explicit special-token exclusion, normalized sequence checks, exact one-based coordinates and a declared long-sequence policy. Emit the [common archive contract](docs/model_input.md), pinned model revision, dimensions, pooling and extraction config. Add a runnable config/example and list exactly which checkpoints/tasks were tested. The [`adapter_usage.py`](examples/adapter_usage.py) example calls the built-in API directly.

Use the existing task loaders, group folds, nested reader selection and scorers. Configure layers and predictors; do not select them on outer-test metrics. Include valid and rejection tests for identities, masks, coordinates and group contamination. Base scoring must remain independent of model-specific dependencies.

## Data and results

Version scientific data changes and record actual source versions, exclusions and provenance. Preserve historical labels, masks, coordinates, groups and model results. Unknown biological targets stay masked; PU learning labels are distinct from biological targets. Audit full biological binding groups and historical use before claiming independent testing.

New runs save actual model revisions, code/input hashes, dependency versions, timestamps, warnings, predictions, selections and fitted parameters. Review results before describing them as accepted benchmark results. Scientific artifact manifests do not register source files or documentation. Data sources and their source-specific terms are documented in [docs/data.md](docs/data.md).
