# DIS.REF.01 independent review

Review began: 2026-10-15T23:59:35+08:00. Review completed: 2026-10-16T00:08:47+08:00. Reviewer: independent subagent stable_control_review.

**Verdict: PASS for the separately versioned, retrospective disorder structure-reference evaluation.** No blocking execution, mask, split, prediction, metric, or parameter-consistency error remains. This is not acceptance of strict direct ordered/disordered ground truth, a prospective/family-independent test, public redistribution, or the whole benchmark.

## Actual computation and provenance

SLURM read-only accounting reports job 820827 COMPLETED, exit 0:0, scheduled 2026-10-15 20:29:28–22:50:19 (02:20:51). The scientific runner records 20:29:34–22:50:18 and 8443.627 seconds. The differing seconds reflect job versus Python-run boundaries. The run used CUDA, Python 3.11.16, NumPy 1.26.4, SciPy 1.13.1, scikit-learn 1.6.1, PyTorch 2.5.1+cu121, CUDA 12.1, with no recorded warnings. Candidate-builder dependency versions are a distinct preparation environment; saved exact fold assignments were consumed during model evaluation.

Local runner and model-interface SHA256 match both the beginning and ending execution hashes recorded on the cluster. All 12 run-input hashes match the local frozen sources or the read-only copied extraction manifest. All 41 registered output hashes match. The frozen candidate preparation manifest matches its registered outputs, and original data_v06 source input hashes match. Independently reread all 1,275 cluster cache archives to SHA256-check them against the exact execution run record; all match, as do all 12 cluster inputs and both code files. No cluster files were changed, no GPU job launched, and no models retrained during review.

## Domain, masks and splits

The selected 1,279 source records match the eligible frozen data_v06 task rows and correspond to 1,275 exact sequences and 1,029 binding groups. All sequence hashes/lengths, complete sequence eligibility (length <=1022), selected residue coordinates, reference-label sources, and original source masks were independently checked. There are 391,718 scored protein-residue rows: 73,766 direct experimental disorder positives, 35 direct experimental order negatives, and 317,917 PDB-observed reference negatives. The 156,160 unknown positions remain excluded from fitting and scoring. Shared sequences permit reuse of vectors, giving 391,189 unique sequence-coordinate vectors; these are different denominators.

Independently reconstructed all outer and inner role tables: every sample is tested exactly once across 5 outer folds; 3 inner validation folds partition only each outer training set. All roles cover their required sample universes, with no binding group split across roles. Identical sequences remain in one binding group. This confirms the registered grouped development split, not complete biological-family isolation or absence of exposure during ESM pretraining.

Streamed every one of the 2,742,026 outer prediction rows. Each row has the correct task, arm, sample, binding group, held-out fold, original label, mask=1, valid 1-based coordinate, finite probability, and unique sample-coordinate-arm identity. All seven arms cover all 391,718 eligible rows exactly once; no unknown coordinate was scored.

## Independent metric recomputation

Computed all 8,953 protein-arm AP/AUROC pairs with a separate threshold-grouped precision integral and rank-pair ROC calculation, without importing the production metric functions. Maximum absolute difference is 6.661338147750939e-16. Independently recomputed every arm primary macro AP, macro AUROC, AP standard deviation, and pooled AP. Reimplemented the 2,000 paired binding-group bootstrap through group sums and sample counts; all reported intervals and deltas match within numerical roundoff (maximum macro-CI difference 5.551115123125783e-16). Intervals condition on fixed outer predictions; they do not include refit or dataset-construction uncertainty.

| Arm | Macro protein AP | 95% group-bootstrap interval |
| --- | ---: | --- |
| ESM2-L11 | 0.722277223 | [0.704348194, 0.740443617] |
| ESM2-L23 | 0.732329698 | [0.713915094, 0.750152319] |
| ESM2-L33 | 0.789855986 | [0.773245197, 0.806096025] |
| traditional-AA-context | 0.585905812 | [0.567035718, 0.605202764] |
| MLP-h33-seed13 | 0.802669763 | [0.785611797, 0.818674296] |
| MLP-h33-seed42 | 0.802767482 | [0.785898215, 0.818654487] |
| MLP-h33-seed2026 | 0.803925855 | [0.786990061, 0.819858049] |

## Fitting, selection and saved parameters

Code dataflow was reviewed for all arms: StandardScaler and class balancing are fit using current training residues only; outer test labels do not enter C selection or MLP early stopping. The 20 linear outer fits use the frozen 3 inner folds and C={0.1,1.0}, averaging held-protein AP with protein-count weighting across folds and selecting the smaller C on ties. The 15 MLP outer fits use the registered 64-unit ReLU/dropout architecture, all three registered seeds, training-only class weights, inner validation AP, and the median of the three best-epoch counts for an independent outer-training refit. All seeds remain reported. Fixed numeric settings in this versioned implementation match the evaluation protocol; the script is not a general executor of arbitrary changed protocol fields.

Checked all 205 selection rows against actual inner domains, validation counts, training-residue counts, weighted aggregate selection values, chosen C values, and median epoch counts. Checked all 35 NPZ parameter archives for model/scaler shapes, finite values, class order, positive scales, revision, arm/fold/seed/epoch identity and agreement with outer-fit records. Independently computed the five traditional-baseline training-domain scaler means/scales; maximum difference is 9.144351942325102e-13.

Reconstructed all 391,718 held-out traditional-baseline probabilities from the saved coefficients and independently generated full-sequence AA/context features: exact agreement. Using 23 exact execution cache archives, reconstructed 7,129 held-out scores for each of the three ESM linear arms and each MLP seed. Linear probabilities agree exactly; MLP differences are at most 4.172325134277344e-07 across CPU NumPy versus original CUDA float32 evaluation. No re-fitting was performed. Inner validation predictions and complete epoch trajectories were not saved; the review checks their selection records and audited dataflow, not an independent numerical reconstruction of every inner fit or best-epoch trajectory. Full ESM prediction reconstruction over all 1,275 sequences and full model-refit reproduction were not performed or claimed.

All four targeted existing tests passed, including positive-evidence/mask rejection, binding-group contamination rejection, both-class metric qualification, and traditional feature definitions. See the preserved test output.

## Source-record differences and local cache finding

Four pairs of distinct DisProt records share an exact sequence. At 12 sequence coordinates, two frozen records carry different reference labels. Their same binding group prevents train/test leakage. These are inherited source-record differences; this review does not adjudicate condition-dependent biology or relabel them as newly discovered errors. A sequence-only deterministic representation cannot distinguish the same sequence-coordinate solely by its record identity. This reference-task ambiguity should be disclosed alongside the distinct protein-record versus unique-vector denominators. The complete records and label sources are in DIS.REF.01_20261015_cross_record_label_differences.json. Frozen source data and evaluated labels were not modified.

The local candidate cache directory contains 23 NPZ-named files; 21 match execution hashes and 2 do not (a809ecfe0a4f8c4d16f0aa28.npz, f46de316667f529b97d5b58a.npz; one is not readable as a ZIP archive). They were not the cluster inputs: all 1,275 actual cluster execution archives match run.json. The two local originals were preserved as found. Only the corresponding two correct cluster archives were copied into logs/reviews/DIS.REF.01_20261015_cache_recovery/ for parameter reconstruction and verified against execution hashes. These two preserved local mismatches must not be described as valid execution caches or used for local resumption; use the matching cluster archives/review copies. This local transfer/cache issue does not invalidate the completed cluster result.

## Review artifacts and reproduction limits

- DIS.REF.01_20261015_verify.py: independent domain, prediction, metric, bootstrap, selection-record and saved-parameter inference verification; run with OPENBLAS_NUM_THREADS=2 .venv_local/bin/python logs/reviews/DIS.REF.01_20261015_verify.py. Requires the preserved local sources, matching review recovery copies, and saved model outputs.
- DIS.REF.01_20261015_verification.json: complete numerical review results, all registered output hashes, local-source hashes and execution code hashes.
- DIS.REF.01_20261015_remote_hashes.json: independent read-only cluster hash results for 1,275 caches, 12 inputs and 2 code files.
- DIS.REF.01_20261015_extract_manifest.tsv: byte-identical read-only cluster manifest copy.
- DIS.REF.01_20261015_local_cache_discrepancies.json: preserved local mismatch names, sizes and old/execution hashes.
- DIS.REF.01_20261015_cross_record_label_differences.json: complete source-record reference differences.
- DIS.REF.01_20261015_tests.txt: targeted test log.

The numerical findings support an accessible residue-level signal in this retrospective reference task. They do not settle strict direct experimental disorder-versus-order discrimination, biological state changes under conditions, or family-independent generalization. MLP gain intervals and all seeds are conditional on the fixed prediction set; there was no best-seed selection.

## Core exact hashes

| File | SHA256 |
| --- | --- |
| `scripts/run_disorder_structure_reference_evaluation_v01.py` | `ade8dbae1a34de020e47595a69ffa5f636f36cd6c071f65fa5ef85d37d811d30` |
| `scripts/physinfobench_model_interface.py` | `c512de6261ebf51cf92fc453a2ab74149c1ca329136cdf5de13d15c9c8e5b7aa` |
| `tests/test_disorder_structure_reference_evaluation.py` | `c930e509ac484b62a9e8237e7516dd6b214e25524552a3797ea9a7516fb56c36` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_protocol_v01.json` | `77df789cf2c9d5767b4b8f38574db461e066194899f1d65055d4664cc4f105aa` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_v01/run.json` | `44ce9b34a37c1e99a46a186b64bc78e526b01fe4f4786932966b6cb77552e34c` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_v01/residue_predictions.tsv.gz` | `66b0fe39d684b592a061a4e0b6e418e2bf5a6da832312fd3c0a929521e72e0ab` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_v01/protein_metrics.tsv` | `4498571c4867280494b0997b31497c6cb06defa27d8c78efab1980833abfe63b` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_v01/metrics.json` | `3bea11cb477bed063f72963518ae13c96c48d80d748f749bc91d6ed28cda2f96` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_v01/inner_selection.tsv` | `64e7485b2d319881395818d00e66c9fd251294565944407af746923f842d8c28` |
| `results/benchmark_fixed/20261015/disorder_structure_reference_candidate_v02/evaluation_v01/outer_fits.tsv` | `44da7b9cddf9d125c7e4da05289b312e8373148600fca42189cf3cf66cdd9849` |
| `logs/reviews/DIS.REF.01_20261015_verify.py` | `7c9cb8b99787cb18e2e7bf091a6569cd38f94ee65fa9ae3792a3875e67a7a8d1` |
| `logs/reviews/DIS.REF.01_20261015_verification.json` | `b35fe171dd98b1defb6e158882adc30a517ff2e15cc75775a052bf5edffb7750` |
| `logs/reviews/DIS.REF.01_20261015_remote_hashes.json` | `950c049f6ee750a485c2d4b9abb4af7365b10b77eb2ad8c34f265f3c5cc9c97b` |

## Closure-record micro-review

Reviewed at 2026-10-16T00:19:29+08:00. PASS. Checked the new DIS.REF.01 task report and only the DIS.REF.01 section at the end of TODO.md. Counts, timing, evaluation scope, independent recomputation, parameter reconstruction and explicit absence of a full refit match this review. Both required disclosures are present: 12 inherited cross-record differing sequence coordinates in four same-sequence record pairs, and two preserved local cache mismatches distinguished from the 1,275 valid cluster execution archives. The report accurately separates structure-reference evaluation from strict DISORDER-RES, family independence, and public-package review. It correctly retains the scientific directory generation-time pending-review state with acceptance recorded separately.

All 41 scientific output hashes and both execution-code hashes remain unchanged. Final task-report SHA256: `f5e68bb5aab70a0d33762873eb401a0a77106d0f75c78ca1846fe7b07355ea27`. The scoped scientific PASS remains in force; no closure-record issue remains.
