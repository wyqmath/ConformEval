# Local tests and reproduction

Install the base package, then invoke tests manually:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

The suite runs without Torch, Transformers or downloaded models. It checks sequence/archive identities, special-token projection, coordinate coverage, long-sequence stitching, label masks, predictions, group isolation, train-only preprocessing and nested candidate selection. The small three-task fixture is marked `software_example`; tests run complete training with layers below 17 and more than one representation dimension.

The replay rescores 67,291 pooled prediction rows and 2,742,026 residue prediction rows. It compares pooled metrics, group intervals, all FS case ranks, and full disorder metrics with 2,000 bootstrap intervals at tolerance `1e-10`. It writes `REPLAY.json` and task scoring records under the requested output directory.

[Model validation](model_validation.md) separately records real-checkpoint inference, small task runs and full benchmark runs. [Resources](resources.md) lists the actual environments, memory, time and cache size for those tests. Archived original code fingerprints identify their source runs; each new extraction/training run records its own actual code/input hashes and configuration. Data/result integrity manifests are scoped to scientific assets; source comments and new files need no repository inventory.
