# ConformEval data card

## Purpose and version

ConformEval examines protein sequence representations for topology, disorder, fold switching and condition-dependent states. This repository snapshot uses software version 0.4.0 and retains the original scientific data versions. It retains the base `data/v0.3.1/` scientific lock and adds separately versioned accepted disorder reference and fold-switch region annotation bundles.

## Data inventory

| Bundle | Contents | Evaluation state |
|---|---|---|
| `data/v0.3.1/` | 11,684 unique sequences; 14,223 task-to-sequence rows; 14,937 task sample rows; strict labels, source references and group folds | KNOT-PRESENCE: 734 samples/283 groups. FS-L2: 8,879 ranking units/8,767 groups. Other task rows retain their individual eligibility status. |
| `data/disorder_reference_v1/` | 1,279 records/1,275 sequences/1,029 groups; 391,718 labeled protein-residue rows | Accepted separate 5 × 3 grouped development evaluation; 2,742,026 predictions across seven arms |
| `data/fs_region_controls_v1/` | 10 experimental pairs/20 endpoints/5,718 positions; 124 inherited positives; 848 primary negatives; 4,746 unknown | Accepted structural-reference annotation. Three dual-class cases; no scoring folds or model result. Strict sensitivity has 580 negatives. |

Counts refer to different units. The task catalog declares those units explicitly. A repeated sequence record or two endpoints of one case are not independent proteins.

## Labels and coordinates

Unknown labels are masked. Strict direct experimental disorder/order annotations remain in the base snapshot. The separate disorder reference task combines 73,766 direct disorder positives with 35 direct ordered and 317,917 experimentally observed PDB reference negatives. It excludes 156,160 unknown positions. A PDB-observed reference position does not establish the state of that residue under every experimental condition.

The disorder reference task has 391,189 unique sequence-coordinate feature vectors. Four identical-sequence record pairs retain 12 coordinates with different source labels across records. All copies stay in the same binding group. See residue_evaluation.md for the coordinate and scoring contract.

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

Sources include DisProt, UniProt, RCSB PDB/SIFTS, KnotProt-derived legacy annotations and primary literature. The base source table is `data/v0.3.1/sources.tsv`; 37 portable case source references resolve through `case_sources.tsv`. Supplementary bundle manifests preserve scientific artifact integrity. [Upstream provenance](../data/provenance/upstream_artifacts.json) records the original source fingerprints. Local provenance paths are replaced by portable references without altering labels or predictions.

The scored knot annotations inherit 2026-09-14/15 legacy topology tables. The separately captured 2026-10-03 KnotProt lists were candidate-audit inputs and do not reconstruct the scored label snapshot. The original legacy source snapshot remains unresolved. This distinction is retained rather than assigning a later download to earlier labels.

The data favor well-studied proteins with published structural or experimental evidence. Label availability differs across tasks. The strict DISORDER-RES task has only one eligible dual-class binding group; KNOT-TYPE lacks a required class. The condition modules and paired input audits retain case-level status; the new condition four-arm experiment is deferred.

## Licensing and distribution

The repository snapshot contains protein sequences, scientific annotations, identifiers and references. It includes no encoder checkpoint weights, cached embeddings, article full text or PDB coordinates. MIT applies to original software. Third-party data retain source-specific terms. Source-specific data terms are recorded below. Repository publication does not supply a blanket data license or upgrade those source-review answers.

## Updating the benchmark

Create a new version for scientific changes. Preserve the frozen labels, coordinates, groups, masks and results. Record new source versions, exclusions and hashes; review the complete biological binding and historical usage before claiming independent evaluation. Update the data card, task catalog and manifests and obtain independent review.


## Sources and licenses

This file records the source review inherited by the v0.4.0 repository snapshot. It does not grant a license. Source attribution and permission to redistribute are tracked separately for database records, derived labels, and literature-derived case annotations.

| Source group | Used for | Recorded terms | Release status |
|---|---|---|---|
| [RCSB PDB data and APIs](https://www1.rcsb.org/pages/usage-policy) | Experimental structures, mappings, and observed-residue evidence | PDB archive data files and programmatic API data are recorded as CC0-1.0; website content has separate terms | Terms recorded for the source inputs; do not extend CC0 to independently sourced annotations or papers |
| [UniProt](https://www.uniprot.org/api-documentation/support-data) | Reference protein sequences | CC BY 4.0 for copyrightable database parts, with attribution | Source-level terms recorded; verify each record's actual provenance |
| [DisProt](https://disprot.org/about) | Disorder annotations | CC BY 4.0, with attribution | Snapshot and date recorded; derived mapping provenance remains in the dataset |
| KnotProt 2.0 | Knot labels and candidate chain sets | No verified redistribution license in the source register | **Unresolved.** Database redistribution terms and original scored-label source lineage have not been verified; no blanket dataset permission is asserted |
| Fold-switching papers and legacy extracts | Paired experimental endpoint labels and case records | Mixed publisher terms; several sources are paywalled or their reuse terms have not been checked | **Unresolved.** No article full text is included; derived factual annotations retain their source-specific review status |
| PYP, RNase A, and BPTI papers | Condition and assembly case modules | Mixed terms, including paywalled sources; one listed PYP article is CC BY 4.0 | **Unresolved by record.** No case-specific redistribution permission is inferred from repository publication |
| ESM-2 checkpoint | Existing reference model identity | Model card records MIT; weights are excluded | No model files are bundled |

The machine-readable source table is [`data/v0.3.1/sources.tsv`](../data/v0.3.1/sources.tsv). It retains the recorded source IDs, versions, terms, and disposition notes while omitting source-local file paths.

Case-level source references are in [`data/v0.3.1/case_sources.tsv`](../data/v0.3.1/case_sources.tsv). Each case row resolves through its `source_ref` and retains the originating artifact hash plus any recorded DOI, PMCID, or citation note. These records improve traceability; they do not resolve the redistribution restrictions summarized above.

## KnotProt snapshot and terms check

The KNOT-PRESENCE scored labels and KNOT-TYPE annotations have a different lineage from the 2026-10-03 KnotProt downloads. The current package inherits `data/curated/knots.tsv`, imported by P1.04 from legacy eligibility, verdict, support, annotation, negative, and matching tables dated 2026-09-14/15. The P1.04 report says it reused those earlier topology outputs; it does not identify a complete original KnotProt file snapshot from which the current labels can be rebuilt. This gap is reflected in the source register: the P1.04 KnotProt record has no verified source snapshot or reuse terms.

The separately captured 2026-10-03 files are `knotted.txt` (7,746 bytes; SHA-256 `b61869c9f9fda9642cad0f7932671e6d6e98c520119bbd080896d8617da93eee`) and `unknotted.txt.gz` (498,035 bytes; SHA-256 `f9bc4824f379418a956c5d201970d1c03587bf78df138df2c8fc6cdccddab32e`). They were used for V.01/V.02 candidate-pool and exposure audits, not as the source snapshot for the current scored labels. The official site describes weekly updates from newly deposited PDB chains, and the downloaded files have no stable release tag.

The official [KnotProt API page](https://knotprot.cent.uw.edu.pl/api) lists the downloadable chain and knot-type data routes. Neither it nor the [KnotProt introduction](https://knotprot.cent.uw.edu.pl/introduction) stated a database reuse license when checked on 2026-10-15. The [KnotProt 2.0 methods article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6323932/) is an open-access CC BY article; its article license does not establish reuse terms for the separately served database files. The 2026-10-03 raw lists and their local project paths are not bundled with the benchmark snapshot.

The source project records the identifiers, dates, hashes, download endpoints, and candidate-audit role of the 2026-10-03 files. Those records do not repair the missing original snapshot for the scored labels or grant redistribution permission. To resolve the source-lineage gap, recover and verify the historical source inputs or create a separately reviewed, versioned annotation from versioned inputs. The current labels and scores must not be silently replaced by a new annotation.

## Included and excluded material

Included in this repository snapshot: curated sequence and label tables, group and fold assignments, source identifiers and citations, reference predictions, and evaluation provenance. Excluded: model weights, cached ESM-2 representations, raw article full text, original PDB coordinate files, and source-specific local paths.

## Version 0.4.0 additions

The separate disorder reference bundle uses the existing DisProt and PDB-observed mapping evidence. The fold-switch region control bundle adds structural-reference annotations from existing experimental pairs. Upstream source fingerprints are in [data provenance](../data/provenance/upstream_artifacts.json); data-generation and evaluation records accompany their scientific artifacts. These additions do not resolve the legacy KnotProt or literature-source redistribution entries above. The original software remains MIT licensed.
