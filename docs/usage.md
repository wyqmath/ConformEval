# Usage

Use Python 3.11+. Run commands from the repository root, with a new output
directory for each run. Install the base package with `python -m pip install .`.
Clone with `git clone -c core.autocrlf=false https://github.com/wyqmath/ConformEval.git`
to preserve the small examples' recorded file hashes.
For an interactive example, open the [notebook](../examples/notebooks/esm2_8m.ipynb)
in [Colab](https://colab.research.google.com/github/wyqmath/ConformEval/blob/main/examples/notebooks/esm2_8m.ipynb).
Its default is the full 734-sample knot task; set `FULL_BENCHMARK=False` for
its 20-sequence synthetic example.

## Run ESM-2 8M

```bash
python -m pip install ".[esm2]"
conformeval-download --model-config configs/models/esm2_8m.json
conformeval-extract --model-config configs/models/esm2_8m.json \
  --input-fasta examples/mini/sequences.fasta \
  --output-dir runs/esm8-features --local-files-only
conformeval-fit-readers --task KNOT-PRESENCE \
  --data-dir examples/mini/base \
  --feature-manifest runs/esm8-features/representation_manifest.json \
  --config configs/readers/linear.json --protocol configs/scoring_protocol.json \
  --output-dir runs/esm8-knot --example
```

Extraction writes 30 archives with 320-dimensional vectors. Training writes
20 predictions per arm, metrics, fitted parameters, and run metadata under
`runs/esm8-knot/`. These synthetic scores test the interface.
For ranking use `--task FS-L2-grouped-PU-ranking`. For residue scoring use
`--task DISORDER-STRUCTURE-REFERENCE-CANDIDATE` and
`--data-dir examples/mini/residue`, omitting `--protocol`.

## Reproduce the reference scores

The checkout includes original evaluator inputs under `data/` and complete
reference predictions and fitted parameters under `results/`.

```bash
conformeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

`runs/reference-replay/REPLAY.json` reports `status: PASS`. The command rescores
saved predictions for all three tasks and compares them with reference metrics.
It does not run a protein encoder or refit a reader.

## Train on a full task

From the repository root:

```bash
conformeval-prepare-inputs --task KNOT-PRESENCE \
  --data-dir data/base \
  --output-dir runs/knot-inputs
conformeval-extract --model-config configs/models/esm2_8m.json \
  --input-fasta runs/knot-inputs/sequences.fasta --sequence-only \
  --output-dir runs/knot-features --local-files-only
conformeval-fit-readers --task KNOT-PRESENCE \
  --data-dir data/base \
  --feature-manifest runs/knot-features/representation_manifest.json \
  --config configs/readers/linear.json --protocol configs/scoring_protocol.json \
  --output-dir runs/knot
```

This task has 569 unique encoder sequences, 734 samples and 283 groups.
Ranking uses the same base directory and its task ID. Disorder uses
`data/disorder_reference`, its residue task ID,
and extraction without `--sequence-only`.

## Score external predictions

The README scoring command uses the complete example files in
[`examples/predictions/`](../examples/predictions/). Sequence-level TSV files
require `task_id`, `sample_id`, `outer_fold`, and `score`. Residue TSV files
require `task_id`, `sample_id`, `position_1based`, `outer_fold`, and `score`.
Use the provided metadata template to identify the actual model and run.
For residues, run `conformeval-evaluate-residues` with `--task-dir`,
`--predictions`, `--model-metadata`, and `--output-dir`; add `--example`
only for software fixtures. See the [evaluation guide](evaluation.md).

## Resources

CPU float32 measurements on 30 short sequences: ESM-2 8M used about 0.5 GiB
peak memory, ESM-2 650M 2.6 GiB, and ProtT5 4.7 GiB. The measured complete
ESM-2 8M three-task run took about 20 minutes, excluding downloads and input
preparation. Long proteins and residue archives can require substantially more
memory and disk. Exact environments remain with each saved reference run.
