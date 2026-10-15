# PhysInfoBench data card

## Purpose and version

PhysInfoBench examines protein sequence representations for topology, disorder, fold switching and condition-dependent states. This repository snapshot uses software version 0.4.0 and derives from the accepted first local delivery. It retains the base `data/v0.3.1/` scientific lock and adds separately versioned accepted disorder reference and fold-switch region annotation bundles.

## Data inventory

| Bundle | Contents | Evaluation state |
|---|---|---|
| `data/v0.3.1/` | 11,684 unique sequences; 14,223 task-to-sequence rows; 14,937 task sample rows; strict labels, source references and group folds | KNOT-PRESENCE: 734 samples/283 groups. FS-L2: 8,879 ranking units/8,767 groups. Other task rows retain their individual eligibility status. |
| `data/disorder_reference_v1/` | 1,279 records/1,275 sequences/1,029 groups; 391,718 labeled protein-residue rows | Accepted separate 5 × 3 grouped development evaluation; 2,742,026 predictions across seven arms |
| `data/fs_region_controls_v1/` | 10 experimental pairs/20 endpoints/5,718 positions; 124 inherited positives; 848 primary negatives; 4,746 unknown | Accepted structural-reference annotation. Three dual-class cases; no scoring folds or model result. Strict sensitivity has 580 negatives. |

Counts refer to different units. The task catalog declares those units explicitly. A repeated sequence record or two endpoints of one case are not independent proteins.

## Labels and coordinates

Unknown labels are masked. Strict direct experimental disorder/order annotations remain in the base snapshot. The separate disorder reference task combines 73,766 direct disorder positives with 35 direct ordered and 317,917 experimentally observed PDB reference negatives. It excludes 156,160 unknown positions. A PDB-observed reference position does not establish the state of that residue under every experimental condition.

The disorder reference task has 391,189 unique sequence-coordinate feature vectors. Four identical-sequence record pairs retain 12 coordinates with different source labels across records. All copies stay in the same binding group. See RESIDUE_EVALUATION.md for the coordinate and scoring contract.

Fold-switch region negatives describe locally preserved structured segments in the two supplied experimental states. All 124 earlier accepted positive endpoint-residue labels are retained. Of 848 primary negative endpoint-residues, 356 occur in the three positive-supported cases. The rest occur in cases without accepted positive region labels. Strict negatives are a subset of primary negatives. Unknown positions are not negative labels.

FS-L2 uses operational PU learning labels to rank known fold-switch cases against an unlabeled background. The background biological target remains unknown with biological mask zero. It is not a verified negative class.

## Splits and historical use

KNOT-PRESENCE and the disorder reference task use five outer folds; FS-L2 uses ten. All use three inner group folds. Preprocessing is fit within the current training partition. The released folds fix each sample's role and prevent its registered binding group from entering both training and test or inner validation.

These are retrospective grouped development evaluations. The declared group checks do not establish complete biological-family independence. The broader group graph still has 37 unresolved PDB-key gaps. Historical source and model usage remain part of the interpretation; these results do not create a new prospective test set.

## Reference model and uncertainty

The encoder is `facebook/esm2_t33_650M_UR50D`, revision `08e4846e537177426273712802403f7ba8261b6c`. The pooled reference compares layers 17/33 and conventional sequence/nearest/prior baselines. The residue reference compares layers 11/23/33, amino-acid/local context and MLP capacity diagnostics with three seeds.

Encoder checkpoint weights and representation archives are excluded. Predictions, selections, fitted reader parameters, environments and hashes are included. Confidence intervals condition on the fixed outer predictions. Reader refitting uncertainty and new-model inference are not part of the reference replay test.

Independent disorder review checked all prediction identities, metrics, parameter archives and actual extraction hashes, reconstructed every traditional-baseline prediction, and numerically sampled 23 actual ESM-2 caches. Two mismatched early local cache copies were excluded; the cluster input hashes were verified. No caches are distributed. This review did not refit every inner reader.

## Sources and traceability

Sources include DisProt, UniProt, RCSB PDB/SIFTS, KnotProt-derived legacy annotations and primary literature. The base source table is `data/v0.3.1/sources.tsv`; 37 portable case source references resolve through `case_sources.tsv`. Supplementary bundle manifests and SOURCE_ARTIFACTS.json preserve the relation between exported files and accepted research artifacts. Local provenance paths are replaced by portable references without altering labels or predictions.

The scored knot annotations inherit 2026-09-14/15 legacy topology tables. The separately captured 2026-10-03 KnotProt lists were candidate-audit inputs and do not reconstruct the scored label snapshot. The original legacy source snapshot remains unresolved. This distinction is retained rather than assigning a later download to earlier labels.

The data favor well-studied proteins with published structural or experimental evidence. Label availability differs across tasks. The strict DISORDER-RES task has only one eligible dual-class binding group; KNOT-TYPE lacks a required class. The condition modules and paired input audits retain case-level status; the new condition four-arm experiment is deferred.

## Licensing and distribution

The repository snapshot contains protein sequences, scientific annotations, identifiers and references. It includes no encoder checkpoint weights, cached embeddings, article full text or PDB coordinates. MIT applies to original software. Third-party data retain source-specific terms. Source-specific data terms remain unresolved as recorded in SOURCES_AND_DISTRIBUTION.md. Repository publication does not supply a blanket data license or upgrade those source-review answers.

## Updating the benchmark

Create a new version for scientific changes. Preserve the frozen labels, coordinates, groups, masks and results. Record new source versions, exclusions and hashes; review the complete biological binding and historical usage before claiming independent evaluation. Update the data card, task catalog and manifests and obtain independent review.
