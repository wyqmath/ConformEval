# Reference results

These are recorded ESM-2 evaluation scores on the benchmark's fixed grouped folds. Use them to compare a model's scores with the reference models and simple baselines.

| Task | Metric | ESM-2 650M | ESM-2 8M | ProtT5 |
|---|---|---:|---:|---|
| Knot presence | AUROC ↑ | 0.868741 | 0.814237 | Not tested |
| Fold-switch ranking | Mean positive rank percentile ↓ | 0.476592 | 0.304406 | Not tested |
| Residue disorder | Mean protein AP ↑ | 0.789856 | 0.738289 | Not tested |

[summary.tsv](summary.tsv) preserves full-precision headline values. [results_index.tsv](results_index.tsv) lists model configurations, baselines, sample counts, environments, and untested tasks. Only these small score tables are kept in the current repository.

The 650M knot and ranking values use nested layer selection; its disorder value is the fixed layer-33 linear reference. The 8M headline values use inner layer selection. For simple baselines, filter `arm_role=baseline` in the index.

## Full reference runs

Full predictions, fitted parameters, uncertainty estimates, and original run records are preserved in the [reference snapshot](https://github.com/wyqmath/ConformEval/tree/5496e710ff202392330e50cb1bf75318d1b3c275/results). The index columns prefixed with `archive_` refer to paths within that snapshot. These files are not required to evaluate a new model.

To reproduce the archived scores, [download the snapshot](https://github.com/wyqmath/ConformEval/archive/5496e710ff202392330e50cb1bf75318d1b3c275.zip), extract it, and supply the extracted directory:

```bash
conformeval-replay-reference --benchmark-root /path/to/extracted-snapshot --output-dir runs/reference-replay
```

The command rescores the saved predictions and writes `runs/reference-replay/REPLAY.json` with `status: PASS` when the archived metrics match. Use a new output directory for each run.

See [evaluation](../docs/evaluation.md) for label semantics and comparison scope.
