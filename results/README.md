# Reference results

| Task | Metric | ESM-2 650M | ESM-2 8M | ProtT5 |
|---|---|---:|---:|---|
| Knot presence | AUROC ↑ | 0.868741 | 0.814237 | Not tested |
| Fold-switch ranking | Mean positive rank percentile ↓ | 0.476592 | 0.304406 | Not tested |
| Residue disorder | Mean protein AP ↑ | 0.789856 | 0.738289 | Not tested |

[summary.tsv](summary.tsv) preserves full-precision headline values.
[results_index.tsv](results_index.tsv) includes every saved arm and baseline,
model/config/environment references, sample counts and untested tasks. Paths
in that table are relative to the repository root. Complete predictions,
parameters, uncertainty estimates and run metadata accompany the indexed runs.

The 650M knot and ranking values use nested layer selection; its disorder
value is the fixed layer-33 linear reference. The 8M headline values use inner
layer selection. Configuration differences are recorded in each indexed run.
For baseline comparisons, filter `arm_role=baseline` in the full index.

Run this command from the repository root to reproduce the archived metrics:

```bash
conformeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

The command rescores the saved predictions and writes `runs/reference-replay/REPLAY.json` with `status: PASS` when the archived metrics match. Use a new output directory for each run.

See [evaluation](../docs/evaluation.md) for label semantics and comparison scope.
