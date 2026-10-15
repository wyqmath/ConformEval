# Source and redistribution review

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

The machine-readable source table is [`data/v0.3.1/sources.tsv`](data/v0.3.1/sources.tsv). It retains the recorded source IDs, versions, terms, and disposition notes while omitting source-local file paths.

Case-level source references are in [`data/v0.3.1/case_sources.tsv`](data/v0.3.1/case_sources.tsv). Each case row resolves through its `source_ref` and retains the originating artifact hash plus any recorded DOI, PMCID, or citation note. These records improve traceability; they do not resolve the redistribution restrictions summarized above.

## KnotProt snapshot and terms check

The KNOT-PRESENCE scored labels and KNOT-TYPE annotations have a different lineage from the 2026-10-03 KnotProt downloads. The current package inherits `data/curated/knots.tsv`, imported by P1.04 from legacy eligibility, verdict, support, annotation, negative, and matching tables dated 2026-09-14/15. The P1.04 report says it reused those earlier topology outputs; it does not identify a complete original KnotProt file snapshot from which the current labels can be rebuilt. This gap is reflected in the source register: the P1.04 KnotProt record has no verified source snapshot or reuse terms.

The separately captured 2026-10-03 files are `knotted.txt` (7,746 bytes; SHA-256 `b61869c9f9fda9642cad0f7932671e6d6e98c520119bbd080896d8617da93eee`) and `unknotted.txt.gz` (498,035 bytes; SHA-256 `f9bc4824f379418a956c5d201970d1c03587bf78df138df2c8fc6cdccddab32e`). They were used for V.01/V.02 candidate-pool and exposure audits, not as the source snapshot for the current scored labels. The official site describes weekly updates from newly deposited PDB chains, and the downloaded files have no stable release tag.

The official [KnotProt API page](https://knotprot.cent.uw.edu.pl/api) lists the downloadable chain and knot-type data routes. Neither it nor the [KnotProt introduction](https://knotprot.cent.uw.edu.pl/introduction) stated a database reuse license when checked on 2026-10-15. The [KnotProt 2.0 methods article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6323932/) is an open-access CC BY article; its article license does not establish reuse terms for the separately served database files. The 2026-10-03 raw lists and their local project paths are not bundled with the benchmark snapshot.

The source project records the identifiers, dates, hashes, download endpoints, and candidate-audit role of the 2026-10-03 files. Those records do not repair the missing original snapshot for the scored labels or grant redistribution permission. To resolve the source-lineage gap, recover and verify the historical source inputs or create a separately reviewed, versioned annotation from versioned inputs. The current labels and scores must not be silently replaced by a new annotation.

## Release decision

The dataset manifest is marked `DRAFT_NOT_CLEARED_FOR_PUBLIC_RELEASE`. The repository `LICENSE` covers the original software only; it does not grant rights to the dataset or source-derived records. The software repository is https://github.com/wyqmath/ConformEval. No DOI has been assigned. No blanket data license is asserted.

The maintainer explicitly requested publication of this complete snapshot in the public repository on 2026-10-16. That publication authorization does not resolve third-party source terms. The original data manifest status is retained, not relabeled as cleared. A future source-cleared version must document the resolved terms or excluded-source subsets and receive a separate manifest and review.

## Included and excluded material

Included in this repository snapshot: curated sequence and label tables, group and fold assignments, source identifiers and citations, reference predictions, and evaluation provenance. Excluded: model weights, cached ESM-2 representations, raw article full text, original PDB coordinate files, and source-specific local paths.

## Version 0.4.0 additions

The separate disorder reference bundle uses the existing DisProt and PDB-observed mapping evidence. The fold-switch region control bundle adds structural-reference annotations from existing experimental pairs. Source/export hashes are in SOURCE_ARTIFACTS.json; reviewed source reports are in audits/. These additions do not resolve the legacy KnotProt or literature-source redistribution entries above. The original software remains MIT licensed.
