# PhysInfoBench

PhysInfoBench tests how well protein sequence representations support prediction of topology, disorder and fold switching. It provides versioned sequence and label tables, frozen group splits, scoring tools, and audited ESM-2 reference predictions.

**Version 0.4.0 is the first benchmark snapshot in the [ConformEval repository](https://github.com/wyqmath/ConformEval).** Three tasks have reference model results. Paired-state region controls and condition case records are included as data modules. The condition four-arm experiment is deferred to a later version.

| Scored task | Evaluation data | Reference result |
|---|---|---|
| KNOT-PRESENCE | 734 samples; 283 binding groups; 5 outer × 3 inner folds | Nested ESM-2: AUROC 0.868741; AP 0.687211 |
| FS-L2 grouped PU ranking | 8,879 ranking units, including 10 positive cases; 10 outer × 3 inner folds | Nested ESM-2: mean held-case rank percentile 0.476592; lower is better |
| DISORDER-STRUCTURE-REFERENCE-CANDIDATE | 1,279 protein records; 1,029 binding groups; 391,718 labeled residues; 5 outer × 3 inner folds | ESM-2 layer 33: mean protein AP 0.789856; traditional baseline 0.585906 |

These are retrospective grouped development results. Unknown biological labels remain masked. In FS-L2, the ranking background is unlabeled. In the disorder reference task, PDB-observed positions supply a structural reference class; the strict direct experimental order/disorder task remains a separate data view.

## Install

Use Python 3.11 or newer. Clone this repository, then install from its root:

```bash
git clone https://github.com/wyqmath/ConformEval.git
cd ConformEval
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install .
```

Windows users can activate `.venv` with `.venv\Scripts\activate`. Base installation supports validation, prediction scoring and the pooled nested readers. The archived residue-reader training additionally requires `python -m pip install ".[reference-residue-training]"`. `requirements-reference.txt` records the exact numerical and PyTorch versions used by that reference computation (Python 3.11.16, PyTorch 2.5.1 with CUDA 12.1); install a PyTorch build appropriate for your compute environment when retraining.

## Verify and reproduce the scores

The following commands run on CPU and use the included predictions. They do not download model weights or fit new readers.

```bash
physinfobench-verify --benchmark-root .
python -m unittest discover -s tests -v
physinfobench-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

The replay checks all seven reference arms for each scored task. It compares KNOT-PRESENCE and FS-L2 point metrics, group-bootstrap intervals and FS case ranks to the archived values, then compares the disorder reference metrics and its 2,000 group-bootstrap intervals. Numerical comparison uses a tolerance of `1e-10`. The output directory must be new. Rescoring and full model retraining are separate reproduction steps; retraining needs model representations that are excluded from this bundle.

## Evaluate a model

1. Read [MODEL_INPUT.md](MODEL_INPUT.md) for sequence identity and representation archive formats.
2. Use the frozen data and groups. Fit preprocessing within each training partition and select readers with inner validation only.
3. For the pooled topology and ranking tasks, run `physinfobench-fit-readers` from exact-sequence representations or submit complete outer-fold predictions to `physinfobench-evaluate`. See [PREDICTION_SUBMISSION.md](PREDICTION_SUBMISSION.md).
4. For disorder residue predictions, use `physinfobench-evaluate-residues`. See [RESIDUE_EVALUATION.md](RESIDUE_EVALUATION.md). The included residue training command reproduces the fixed ESM-2 reference configuration; another encoder supplies its own frozen OOF predictions and inference provenance.

Scores identify which samples were evaluated and how their metrics were computed. Training records, model revisions and preprocessing provenance accompany the predictions. A score table alone does not prove how the model was trained.

## Contents

| Path | Contents |
|---|---|
| `data/v0.3.1/` | Frozen base dataset, sequence mappings, source records, strict labels and pooled task folds |
| `data/disorder_reference_v1/` | Accepted separate disorder reference task, exact sequences, selected residue labels and frozen 5 × 3 roles |
| `data/fs_region_controls_v1/` | Ten experimental pairs; full-length labels, provenance and stable reference intervals |
| `results/esm2_existing_v0.3/` | Pooled ESM-2 reference predictions, metrics, selections and fitted readers |
| `results/disorder_reference_esm2_v1/` | Residue predictions, per-protein metrics, 35 fitted readers and run metadata |
| `scripts/`, `tests/` | CLI implementations and checks for identities, masks, coordinates and split contamination |
| `task_catalog.tsv` | All 13 registered tasks with counts, status and result links |
| `manifest.json`, `SOURCE_ARTIFACTS.json` | Delivery hashes and the relation to reviewed research artifacts |
| `audits/`, `VALIDATION.md` | Review records and local delivery validation |

The fold-switch region module preserves 124 positive endpoint-residue labels and adds 848 stable structural-reference negatives, with 580 negatives in the stricter sensitivity view. Three dual-class cases contain 356 primary negative endpoint-residues. This is an annotation release; region model scoring has not run. KNOT-TYPE, strict DISORDER-RES and condition modules retain the status recorded in the catalog. No new encoder or structure model result is claimed.

## License and data access

Original benchmark software is licensed under the [MIT License](LICENSE), copyright 2026 Jiayuan Ma. Third-party data retain their source terms; [SOURCES_AND_DISTRIBUTION.md](SOURCES_AND_DISTRIBUTION.md) records the database and literature sources. Model weights, embedding caches, raw papers and PDB coordinate files are excluded.

Source-specific redistribution terms remain **unresolved** for legacy KnotProt-derived annotation lineage/terms and some literature-derived tables. Publishing this snapshot does not establish those permissions or grant a blanket dataset license. The original data manifests retain their recorded source-review status. See [UPLOAD_GUIDE.md](UPLOAD_GUIDE.md) and `PUBLICATION.json` for this repository export and its relation to the accepted local delivery.

## Citation and contributions

[CITATION.cff](CITATION.cff) contains the software citation. Also cite the primary data sources used by each task. The repository is https://github.com/wyqmath/ConformEval. No DOI has been assigned. See [CONTRIBUTING.md](CONTRIBUTING.md) for validation and versioning rules, [DATA_CARD.md](DATA_CARD.md) for data details, and [CHANGELOG.md](CHANGELOG.md) for this delivery's changes.
