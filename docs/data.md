# Data

## Published dataset

[Zenodo](https://zenodo.org/records/23433394) · [DOI](https://doi.org/10.5281/zenodo.23433394)
· [ConformEval-data.zip](https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1).

```bash
curl --fail --location \
  "https://zenodo.org/records/23433394/files/ConformEval-data.zip?download=1" \
  --output ConformEval-data.zip
unzip -q ConformEval-data.zip -d datasets
```

The published archive is unchanged. It extracts into `datasets/ConformEval-data/`.

| Directory | Data | Size |
|---|---|---:|
| `knot_presence/` | Sequence-level knot labels and splits | 734 samples |
| `fold_switching/` | Confirmed cases and unlabeled ranking background | 8,879 items |
| `disorder/` | Residue disorder/structure-reference labels and splits | 1,279 records; 391,718 labeled residues |

Each task includes `sequences.fasta`, `samples.tsv`, `test_splits.tsv`, and
`validation_splits.tsv`. Disorder also includes `labels.tsv`. The archive's
README defines every field; `sources.tsv` records its sources and terms.
Tables use tab-separated UTF-8 text. Read them with
`csv.DictReader(handle, delimiter="\t")` or `pandas.read_csv(path, sep="\t")`.

## Evaluator files

The checkout includes evaluator inputs under `data/` and evaluation settings under `configs/`. The Zenodo export has a separate reader-friendly layout.

| Path | Purpose |
|---|---|
| `data/base/` | Sequences, task labels, source records, and grouped folds |
| `data/disorder_reference/` | Residue labels and grouped folds |
| `data/fs_region_controls/` | Paired-structure annotations and evidence |

Compressed `.tsv.gz` files use the same fields as TSV. Open them with
`gzip.open(path, "rt", encoding="utf-8", newline="")` or use
`pandas.read_csv(path, sep="\t")`, which detects gzip automatically.
Task manifests list the required files; the included README defines fields and coordinates. The [ConformEval tutorial](../examples/notebooks/ConformEval_tutorial.ipynb) notebook provides inference and training commands.

## Label definitions

- Knot presence is a whole-chain label: 112 knotted and 622 unknotted samples.
- Fold-switch ranking has 10 confirmed cases and 8,869 unlabeled background items.
  An unlabeled item is not an experimentally verified non-switching protein.
- Disorder combines 73,766 experimental disorder labels, 35 direct ordered
  labels and 317,917 PDB-observed structure-reference labels. Unknown positions
  remain masked. Coordinates are one-based positions in the declared sequence.

The [task catalog](../data/task_catalog.tsv) also records tasks without model
scores. Group and coordinate rules are in the [evaluation guide](evaluation.md).

## Sources

Sequences and annotations come from UniProt, DisProt, RCSB PDB/SIFTS,
KnotProt-derived legacy annotations, and primary publications. The source records
retain source IDs, record versions, citations, and exclusions.
MIT covers the software; it does not grant a blanket dataset license.
UniProt and DisProt attribution and PDB source terms are recorded in the source
register. KnotProt redistribution terms and the original scored-label source
snapshot remain unresolved; literature-derived records retain their recorded
source-specific terms. See [`sources.tsv`](../data/base/sources.tsv),
[`case_sources.tsv`](../data/base/case_sources.tsv), and the
[archived source review](https://github.com/wyqmath/ConformEval/blob/3c575aad6e7a4cff3dd9cd0237948dcce17bf22d/docs/data.md#sources-and-licenses).
