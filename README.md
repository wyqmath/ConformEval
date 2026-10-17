# ConformEval

ConformEval evaluates **protein language model representations** on three scoring tasks. It provides fixed sequences and labels, biological binding-group folds, model adapters, nested reader training, scoring tools and reference predictions. Models use the same task samples, group assignments and scoring implementation.

## Task and model results

The catalog registers **13 tasks**. **Three have complete biological model scores**; the other ten are marked **Not tested** below. Small synthetic examples are interface tests and are excluded from these performance tables. The full [result index](results/results_index.tsv) and [readable results](results/README.md) include every saved arm, configuration, sample unit, environment and prediction path.

| Task / primary metric | Data | ESM-2 650M | ESM-2 8M | ProtT5 |
|---|---|---|---|---|
| KNOT-PRESENCE / AUROC ↑ | 734 samples | 0.868741 (nested L17/33) | 0.814237 (nested L3/6) | Not tested |
| FS-L2 / mean held-case rank percentile ↓ | 8,879 items; 10 positives | 0.476592 (nested L17/33) | 0.304406 (nested L3/6) | Not tested |
| Disorder structure reference / mean protein AP ↑ | 1,279 records; 391,718 rows | 0.789856 (fixed L33) | 0.738289 (nested L3/6) | Not tested |

Each cell names its reader/layer rule. The 650M disorder value is the fixed L33 linear reference, not a nested-layer result. Its three MLP seeds are all listed in the [full results](results/README.md). This matrix reports existing experiments; their different layer grids do not constitute a controlled model-size comparison. Higher is better for AUROC/AP; lower is better for the FS rank percentile. The FS background is biologically unlabeled.

Knot presence uses a whole-chain label. FS ranks confirmed switching cases against unlabeled background items. The disorder reference combines direct experimental disorder labels with direct ordered labels and positions observed in experimental PDB structures as the ordered reference.

### All 13 registered tasks

| Task | Input | Current data | Metric | Status / data needed |
|---|---|---|---|---|
| `KNOT-PRESENCE` | Full sequence | 734 samples; 283 groups | AUROC ↑ | Tested: 650M and 8M; ProtT5 Not tested. 5 outer × 3 inner folds; evaluated |
| `KNOT-TYPE` | Full sequence | 65 eligible: 3_1=57, 4_1=6, 5_2=2 | Macro-F1 on locked classes ↑ | Not tested. Required 5_1 class is absent; 4_1 and 5_2 have only 6 and 2 examples. No scoring folds |
| `FS-L2-grouped-PU-ranking` | Two endpoint sequences per positive case; one per background item | 8,879 items; 10 positive cases; 8,767 groups | Mean held-case rank percentile ↓ | Tested: 650M and 8M; ProtT5 Not tested. 10 outer × 3 inner folds; evaluated |
| `DISORDER-RES` | Sequence + residue coordinates | 2,279 registered proteins; 1 eligible dual-class binding group | Mean protein AP ↑ | Not tested. Additional proteins/groups with direct experimental order and disorder labels; locked protocol requires ≥4 dual-class groups. No scoring folds |
| `DISORDER-STRUCTURE-REFERENCE-CANDIDATE` | Sequence + residue coordinates | 1,279 records; 1,029 groups; 391,718 labeled rows | Mean protein AP ↑ | Tested: 650M and 8M; ProtT5 Not tested. 5 outer × 3 inner folds; evaluated |
| `FS-L1-PAIRED` | Paired experimental sequences and structures | 16 candidates; 10 accepted cases / 20 endpoints | Case evidence and input-identity audit; no model metric | Not tested. Six pending candidates lack accepted endpoint inputs/evidence in the frozen sequence frame. Case report; no prediction target/scoring folds defined |
| `FS-REGION` | Sequence + switching-region coordinates | 10 cases / 20 endpoints; 124 positive + 848 reference-negative rows (580 strict); 3 dual-class cases | Label coverage; no model metric defined | Not tested. Switching-region positives cover 3/10 cases; 4,746 endpoint positions remain unknown. Annotations available; scoring protocol/folds not defined; model evaluation not run |
| `FS-L3-MATCHED` | Sequences of cases and matched controls | 0 qualified matched case-control pairs | AUROC/AP after qualification ↑ | Not tested. Experimentally qualified matched controls with usable biological labels. Control-label qualification pending |
| `BPTI-REDOX` | Same-variant sequence + redox conditions | 0 qualified experiment pairs | Not defined after qualification stop | Not tested. Qualified endpoint measurements and matched-state evidence. Prespecified qualification stop |
| `PYP-STATE` | Sequence + illumination condition | 2 state records | Case/input audit; no prediction metric defined | Not tested. Additional independent groups with matched illumination/state measurements. Case report; no independent scoring folds |
| `RNASE-ASSEMBLY` | Sequence + assembly condition | 4 state records | Case/input audit; no prediction metric defined | Not tested. Independent same-batch assembly pairs with measured conditions; current pairs carry cross-batch warnings. Case report; no independent scoring folds |
| `DISORDER-CONDITION-pH` | Same sequence/region + pH | 7 state records | Case/input audit; no prediction metric defined | Not tested. Additional qualified independent condition pairs; candidate records include cross-publication and region mismatches. Case report; no independent scoring folds |
| `DISORDER-CONDITION-temperature` | Same sequence/region + temperature | 2 state records | Case/input audit; no prediction metric defined | Not tested. Additional independent groups with matched two-temperature experimental labels. Case report; no independent scoring folds |

FS-REGION counts are endpoint-residue rows: 848 reference negatives across all ten cases, including 356 negatives in the three dual-class cases; 580 is the strict subset. Condition counts are state records, not independent proteins or cases. The strict direct-label disorder task and the structural-reference task retain separate labels and results.

### Simple baselines on the same held groups

| Task | Simple baseline | 650M reference run | 8M run |
|---|---|---|---|
| KNOT-PRESENCE | Length logistic | 0.675945 | 0.675945 |
| KNOT-PRESENCE | Composition logistic (650M: length + AA20; 8M: AA20) | 0.819491 | 0.785140 |
| KNOT-PRESENCE | Representation 1-nearest neighbor | 0.672987 | 0.724348 |
| KNOT-PRESENCE | Training prior | 0.494323 | 0.494323 |
| FS-L2-grouped-PU-ranking | Length logistic | 0.711773 | 0.711773 |
| FS-L2-grouped-PU-ranking | Composition logistic (650M: length + AA20; 8M: AA20) | 0.498057 | 0.480346 |
| FS-L2-grouped-PU-ranking | Representation 1-nearest neighbor | 0.500728 | 0.500737 |
| FS-L2-grouped-PU-ranking | Training prior | 0.500564 | 0.500564 |
| Disorder structure reference | AA + position + local context logistic | 0.585906 | 0.585906 |

The columns identify source runs. Length, composition, training prior and traditional residue features do not use a language model. The 650M composition baseline uses log-length plus 20 amino-acid frequencies (21 features); the 8M composition baseline uses only the 20 frequencies. The nearest-neighbor baseline uses that run's highest candidate layer (33 or 6). All arms and full-precision values, including fixed-layer diagnostics, appear in [the result index](results/results_index.tsv). Prediction intervals remain in the linked metric files.

Data definitions and sources are documented in [data documentation](docs/data.md) and the [task catalog](data/task_catalog.tsv).

## Supported models and inputs

| Model | Adapter | Transformer layers / dimension | Tested use |
|---|---|---|---|
| ESM-2 650M | `esm2`, [`esm2_650m.json`](configs/models/esm2_650m.json) | 33 / 1,280 | Archived reference scores for all three tasks; real-weight extraction and all three small task examples |
| ESM-2 8M | `esm2`, [`esm2_8m.json`](configs/models/esm2_8m.json) | 6 / 320 | Real-weight extraction with layers 3/6; all three small examples and all three complete tasks passed; [results](docs/model_validation.md) |
| ProtT5 encoder | `prott5`, [`prott5.json`](configs/models/prott5.json) | 24 / 1,024 | Real-weight extraction and all three small tasks passed; complete biological tasks not run; [details](docs/model_validation.md) |

Inputs are an exact-sequence FASTA and a pinned model configuration. Extraction writes sequence/residue NPZ archives and `representation_manifest.json`. Training reads those archives, task labels and fixed folds, and writes out-of-fold predictions, fitted parameters, selected settings, metrics and `run.json`. Complete external prediction TSVs can also be scored directly. See [the archive and adapter contract](docs/model_input.md).

Base scoring and linear/ridge reader training require Python **3.11+**, NumPy, SciPy and scikit-learn. Actual operating-system installation and scoring coverage is listed in the [platform matrix](docs/platforms.md); Windows is **Not tested**. Model dependencies are optional: install `.[esm2]` for ESM-2 or `.[prott5]` for ProtT5. See [installation](docs/installation.md) for tested versions and [resources](docs/resources.md) for measured time, peak memory, cache sizes and long-sequence handling.

## Tutorial 1: reproduce reference scores from included predictions

This route needs the base package and the included data/results:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv
source .venv/bin/activate
python -m pip install .
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

Expected output: `runs/reference-replay/REPLAY.json` has `status: PASS`. All seven arms for each of the three tasks match the archived metrics and group-bootstrap intervals at tolerance `1e-10`. The table above shows the reference primary results. The output directory must be new. Windows activation commands and their validation status are in [installation](docs/installation.md).

For a small, immediate external-prediction example, the actual TSV and metadata are included:

```bash
confomeval-evaluate   --task KNOT-PRESENCE   --data-dir examples/mini/base   --predictions examples/predictions/knot.tsv   --model-metadata examples/predictions/metadata.json   --protocol configs/scoring_protocol.json   --output-dir runs/external-knot-example
```

Expected: `runs/external-knot-example/metrics.json` reports AUROC **1.0** for 20 synthetic samples. This checks the scoring interface. [Residue predictions and metadata templates](docs/prediction_submission.md) provide the equivalent coordinate-level example.

## Tutorial 2: extract representations, train readers and score a model

The following complete example uses the real ESM-2 8M checkpoint and 30 short synthetic sequences. It exercises all three task interfaces using the included small group folds. Start from a fresh clone and a separate environment:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv-model
source .venv-model/bin/activate
python -m pip install ".[esm2]"
confomeval-download --model-config configs/models/esm2_8m.json
confomeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta examples/mini/sequences.fasta   --output-dir runs/esm8-mini-features   --local-files-only

confomeval-fit-readers   --task KNOT-PRESENCE --data-dir examples/mini/base   --feature-manifest runs/esm8-mini-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-mini-knot --example

confomeval-fit-readers   --task FS-L2-grouped-PU-ranking --data-dir examples/mini/base   --feature-manifest runs/esm8-mini-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-mini-fs --example

confomeval-fit-readers   --task DISORDER-STRUCTURE-REFERENCE-CANDIDATE   --data-dir examples/mini/residue   --feature-manifest runs/esm8-mini-features/representation_manifest.json   --config configs/readers/linear.json   --output-dir runs/esm8-mini-residue --example
```

Expected: extraction reports **30 sequences, 320 dimensions**, with layers **3 and 6**. Each training command ends with `status: PASS`. Each pooled task writes 20 out-of-fold rows per arm; the residue task writes **554 labeled rows per arm**. Outputs are `predictions.tsv.gz`, `model_metadata.json`, `metrics.json`, `inner_selection.tsv`, `outer_selections.tsv`, `parameters/` and `run.json` in the requested directories. Synthetic scores have no biological performance interpretation.

The model config chooses layers, pooling and long-sequence policy. The reader config chooses candidate predictors and regularization. Inner validation selects settings; outer groups are used only for the final predictions. For ESM-2 650M use `configs/models/esm2_650m.json`. For ProtT5 install `.[prott5]` and use `configs/models/prott5.json`, with its residue-spaced tokenizer and 1,024-dimensional outputs. Run each model in its own output directories.

For the complete benchmark, use the full task data rather than the synthetic fixture. This is a complete single-task example:

```bash
confomeval-prepare-inputs   --task KNOT-PRESENCE --data-dir data/v0.3.1   --output-dir runs/knot-inputs
confomeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta runs/knot-inputs/sequences.fasta   --output-dir runs/esm8-knot-features --sequence-only --local-files-only
confomeval-fit-readers   --task KNOT-PRESENCE --data-dir data/v0.3.1   --feature-manifest runs/esm8-knot-features/representation_manifest.json   --config configs/readers/linear.json --protocol configs/scoring_protocol.json   --output-dir runs/esm8-knot
```

Expected: **569 unique sequences**, **734 samples**, **283 groups**, and **5 × 3** folds. [Full-task commands](docs/full_evaluation.md) cover PU ranking and generic residue training. The [small NumPy example](examples/README.md) also exercises training without model downloads.

## Repository and local checks

`src/confomeval` contains importable code; `configs` holds protocols/model settings; `data` holds versioned task data; `results` contains actual reference results; `examples` contains runnable examples; `docs` contains usage/source information. All checks run locally:

```bash
python -m unittest discover -s tests -v
```

Cite the primary data sources for scientific use. [LICENSE](LICENSE) covers original software; [data documentation](docs/data.md) records actual source terms. See [CONTRIBUTING.md](CONTRIBUTING.md) for adding models, data or results.
