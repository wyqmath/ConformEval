# Full task evaluation

Install the appropriate model extra, then download the pinned model selected in your config. Run each task independently. All models consume the same saved group folds and scorer.

## KNOT-PRESENCE

```bash
confomeval-download --model-config configs/models/esm2_8m.json
confomeval-prepare-inputs --task KNOT-PRESENCE --data-dir data/base --output-dir runs/knot-inputs
confomeval-extract --model-config configs/models/esm2_8m.json --input-fasta runs/knot-inputs/sequences.fasta --output-dir runs/esm8-knot-features --sequence-only --local-files-only
confomeval-fit-readers --task KNOT-PRESENCE --data-dir data/base --feature-manifest runs/esm8-knot-features/representation_manifest.json --config configs/readers/linear.json --protocol configs/scoring_protocol.json --output-dir runs/esm8-knot
```

Expected coverage: 569 unique sequences, 734 samples, 283 groups, 5 outer and 3 inner folds. The primary metric is pooled out-of-fold AUROC; AP is secondary.

## FS-L2 grouped PU ranking

```bash
confomeval-prepare-inputs --task FS-L2-grouped-PU-ranking --data-dir data/base --output-dir runs/fs-inputs
confomeval-extract --model-config configs/models/esm2_8m.json --input-fasta runs/fs-inputs/sequences.fasta --output-dir runs/esm8-fs-features --sequence-only --local-files-only
confomeval-fit-readers --task FS-L2-grouped-PU-ranking --data-dir data/base --feature-manifest runs/esm8-fs-features/representation_manifest.json --config configs/readers/linear.json --protocol configs/scoring_protocol.json --output-dir runs/esm8-fs
```

Expected coverage: 8,888 unique endpoint/background sequences, 8,879 ranking items including 10 positive cases, and 10 outer folds. Each positive case is ranked only within its assigned outer test pool. A case vector averages its two validated endpoint representations. Unknown background biological targets stay masked.

## Disorder structural reference

```bash
confomeval-prepare-inputs --task DISORDER-STRUCTURE-REFERENCE-CANDIDATE --data-dir data/disorder_reference --output-dir runs/residue-inputs
confomeval-extract --model-config configs/models/esm2_8m.json --input-fasta runs/residue-inputs/sequences.fasta --output-dir runs/esm8-residue-features --local-files-only
confomeval-fit-readers --task DISORDER-STRUCTURE-REFERENCE-CANDIDATE --data-dir data/disorder_reference --feature-manifest runs/esm8-residue-features/representation_manifest.json --config configs/readers/linear.json --output-dir runs/esm8-residue
```

Expected coverage: 1,275 sequences, 1,279 protein records, 391,718 labeled rows and 5 × 3 group folds. Residue archives must include every labeled position. Unknown residues are excluded. The primary selection and scoring metric is mean per-protein AP.

## Configuration and outputs

Model configs declare layers and pooling before extraction. Reader configs declare candidate predictors (`logistic`, `ridge`) and `C_grid`; ridge uses `alpha=1/C`. Optional `candidate_layers` in the reader config selects a subset of the extracted layers. Residue validation aggregates AP equally over all inner held proteins. Numerical ties are measured within `1e-12` of the exact best validation metric, with deterministic ties: smaller C, shallower layer and predictor name. The final outer predictions are never used to select these settings.

Each task writes `predictions.tsv.gz`, `metrics.json`, `model_metadata.json`, `inner_selection.tsv`, `outer_selections.tsv`, fitted `parameters/` and a provenance `run.json`. Pooled tasks store all compact arm predictions in one `predictions.tsv.gz` with an `arm` column; select an arm by filtering that column. Pooled runs also retain enriched scorer tables under `scores/<arm>/` for label, mask and binding-group audits. For residue tasks `protein_metrics.tsv` records equal-weight per-protein aggregation. New output directories are required.

See [resources](resources.md) for actual timings and memory. The original ESM-2 650M reference configurations remain with its run outputs. Model performance comparisons must identify their predeclared layer/predictor grids and representation extraction settings.
