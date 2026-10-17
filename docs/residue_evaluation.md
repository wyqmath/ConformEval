# Generic residue evaluation

The scored task `DISORDER-STRUCTURE-REFERENCE-CANDIDATE` contains 1,279 protein records, 1,275 unique sequences, 1,029 groups and 391,718 known protein-residue rows. Labels include 73,766 direct disorder positives, 35 direct order negatives and 317,917 PDB-observed reference negatives; 156,160 positions remain unknown. The strict direct order/disorder view is a separate task.

## Train a residue reader for any supported adapter

```bash
confomeval-prepare-inputs   --task DISORDER-STRUCTURE-REFERENCE-CANDIDATE   --data-dir data/disorder_reference --output-dir runs/residue-inputs
confomeval-extract   --model-config configs/models/esm2_8m.json   --input-fasta runs/residue-inputs/sequences.fasta   --output-dir runs/esm8-residue-features
confomeval-fit-readers   --task DISORDER-STRUCTURE-REFERENCE-CANDIDATE   --data-dir data/disorder_reference   --feature-manifest runs/esm8-residue-features/representation_manifest.json   --config configs/readers/linear.json --output-dir runs/esm8-residue
```

Use any supported model config, including ProtT5. The reader consumes the supplied layers and dimension; there is no global layer requirement. The same command supports logistic or ridge candidates from the reader config. Scaling is fit on the current training residues. Candidate layers, predictors and C values are selected by mean per-protein AP over the saved inner held proteins; the outer reader is refit without outer test groups.

Outputs include complete residue predictions for every known coordinate, per-protein AP/AUROC, mean protein AP, pooled secondary AP, binding-group intervals, all inner scores, selected settings and fitted parameters. The traditional reference baseline uses amino-acid one-hot, relative position and local 15-residue composition.

The `--example` flag selects the explicitly synthetic bundle in `examples/mini/residue`. Without it the loader requires the versioned benchmark data. Four pairs in the full data share exact sequences, with 12 coordinates having source-label differences across records; all copies share a binding group. Each record is scored separately. There are 391,189 unique sequence-coordinate feature vectors behind the scoring rows.

## Score external predictions or reproduce archived scores

See [prediction submission](prediction_submission.md) for a complete residue TSV and metadata example. To score the archived seven ESM-2 reference arms:

```bash
confomeval-evaluate-residues   --task-dir data/disorder_reference   --predictions results/disorder_reference_esm2/residue_predictions.tsv.gz   --reference --output-dir runs/disorder-reference-replay
```

The archived reference includes ESM-2 layers 11/23/33, a traditional baseline and three MLP seeds. Their original configurations, parameters and numerical environments remain in `results/disorder_reference_esm2`. New adapter runs use their own configured candidates and separate output directories. Reference rescoring checks the archived predictions; resource measurements and current-model test coverage are in [resources](resources.md) and [validation](model_validation.md).
