import csv
import gzip
import io
from collections import defaultdict
from pathlib import Path

TASK_KNOT = "KNOT-PRESENCE"
TASK_FS = "FS-L2-grouped-PU-ranking"
TASK_DISORDER = "DISORDER-STRUCTURE-REFERENCE-CANDIDATE"
DIRECTORIES = {
    TASK_KNOT: "knot_presence",
    TASK_FS: "fold_switching",
    TASK_DISORDER: "disorder",
}


def read_table(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames) != len(
            set(reader.fieldnames)
        ):
            raise ValueError(f"missing or duplicate table columns: {path}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"inconsistent table columns: {path}")
    return rows


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    def write(stream):
        writer = csv.DictWriter(
            stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    if path.suffix == ".gz":
        with path.open("wb") as raw:
            with gzip.GzipFile(
                filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9
            ) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                    write(stream)
    else:
        with path.open("w", encoding="utf-8", newline="") as stream:
            write(stream)


def parse_fasta(path: Path) -> dict[str, str]:
    records = {}
    current = None
    with path.open(encoding="utf-8-sig") as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()
                if not name or name[0] in records:
                    raise ValueError(f"empty or duplicate FASTA identifier: {line}")
                current = name[0]
                records[current] = []
            elif current is None:
                raise ValueError("FASTA sequence precedes its header")
            else:
                records[current].append(line.upper())
    sequences = {sid: "".join(parts) for sid, parts in records.items()}
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError("FASTA is empty or contains an empty sequence")
    return sequences


def task_directory(root: Path, task: str) -> Path:
    if task not in DIRECTORIES:
        raise ValueError(f"unsupported task: {task}")
    return root if (root / "samples.tsv").is_file() else root / DIRECTORIES[task]


def sequence_path(root: Path, task: str) -> Path:
    directory = task_directory(root, task)
    path = directory / "sequences.fasta"
    return path if path.is_file() else directory.parent / "sequences.fasta"


def residue_label_path(directory: Path) -> Path:
    path = directory / "labels.tsv"
    return path if path.is_file() else directory / "labels.tsv.gz"


def task_input_paths(root: Path, task: str) -> list[Path]:
    directory = task_directory(root, task)
    paths = [sequence_path(root, task)] + [
        directory / name
        for name in ["samples.tsv", "test_splits.tsv", "validation_splits.tsv"]
    ]
    if task == TASK_DISORDER:
        paths.append(residue_label_path(directory))
    return paths


def split_assignments(
    directory: Path, samples: dict, n_outer: int
) -> tuple[dict, dict]:
    test_folds = {}
    groups = {sid: row["binding_group_id"] for sid, row in samples.items()}
    group_test = {}
    for row in read_table(directory / "test_splits.tsv"):
        sid, fold = row["sample_id"], int(row["test_fold"])
        if sid not in samples or sid in test_folds or fold not in range(n_outer):
            raise ValueError(f"invalid or duplicate test split: {sid}")
        test_folds[sid] = fold
        if group_test.setdefault(groups[sid], fold) != fold:
            raise ValueError(f"binding group crosses outer folds: {groups[sid]}")
    if set(test_folds) != set(samples) or set(test_folds.values()) != set(
        range(n_outer)
    ):
        raise ValueError("test splits must cover all samples and all outer folds")
    inner_groups = defaultdict(dict)
    inner_samples = defaultdict(set)
    for row in read_table(directory / "validation_splits.tsv"):
        outer, sid, inner = (
            int(row["test_fold"]),
            row["sample_id"],
            int(row["validation_fold"]),
        )
        if outer not in range(n_outer) or inner not in range(3) or sid not in samples:
            raise ValueError("invalid validation split identity")
        if test_folds[sid] == outer:
            raise ValueError(f"outer test sample entered inner validation: {sid}")
        if sid in inner_samples[outer]:
            raise ValueError(f"duplicate validation assignment: {outer}/{sid}")
        inner_samples[outer].add(sid)
        if inner_groups[outer].setdefault(groups[sid], inner) != inner:
            raise ValueError(f"binding group crosses inner folds: {groups[sid]}")
    for outer in range(n_outer):
        training_samples = {sid for sid, fold in test_folds.items() if fold != outer}
        if inner_samples[outer] != training_samples or set(
            inner_groups[outer].values()
        ) != {0, 1, 2}:
            raise ValueError(
                f"validation splits do not partition outer training: {outer}"
            )
    return test_folds, dict(inner_groups)


def load_samples(root: Path, task: str, example=False) -> tuple[dict, dict]:
    directory = task_directory(root, task)
    sequences = parse_fasta(sequence_path(root, task))
    samples = {}
    sequence_groups = defaultdict(set)
    for raw in read_table(directory / "samples.tsv"):
        if not {"sample_id", "group_id"} <= raw.keys():
            raise ValueError("samples require sample_id and group_id columns")
        kind = raw.get("data_kind", "benchmark")
        if kind not in {"benchmark", "synthetic"}:
            raise ValueError(f"unsupported data_kind: {kind}")
        if (kind == "synthetic") != example:
            raise ValueError(
                "synthetic samples require --example; benchmark samples must omit it"
            )
        sid, group = raw["sample_id"], raw["group_id"]
        if not sid or not group or sid in samples:
            raise ValueError(f"empty or duplicate sample identity: {sid}")
        for plural, singular in [
            ("sequence_ids", "sequence_id"),
            ("lengths", "length"),
        ]:
            if plural in raw and singular in raw and raw[plural] != raw[singular]:
                raise ValueError(f"conflicting {singular}/{plural} fields: {sid}")
        sequence_ids = raw.get("sequence_ids", raw.get("sequence_id", "")).split(";")
        lengths = raw.get("lengths", raw.get("length", "")).split(";")
        if not sequence_ids or len(sequence_ids) != len(lengths):
            raise ValueError(f"invalid sequence mapping: {sid}")
        for sequence_id, length in zip(sequence_ids, lengths):
            if sequence_id not in sequences or int(length) != len(
                sequences[sequence_id]
            ):
                raise ValueError(
                    f"sequence identity or length mismatch: {sid}/{sequence_id}"
                )
            sequence_groups[sequences[sequence_id]].add(group)
        if task != TASK_FS and len(sequence_ids) != 1:
            raise ValueError(f"task requires one sequence per sample: {sid}")
        samples[sid] = {
            **raw,
            "task_id": task,
            "binding_group_id": group,
            "sequence_ids": ";".join(sequence_ids),
            "lengths": ";".join(lengths),
        }
        if len(sequence_ids) == 1:
            samples[sid].update(sequence_id=sequence_ids[0], length=lengths[0])
    if not samples:
        raise ValueError("sample table is empty")
    if any(len(groups) != 1 for groups in sequence_groups.values()):
        raise ValueError("identical sequences cross binding groups")
    return samples, sequences


def load_pooled(root: Path, tasks=None, example=False) -> tuple[dict, dict, dict]:
    tasks = {TASK_KNOT, TASK_FS} if tasks is None else set(tasks)
    if not tasks or not tasks <= {TASK_KNOT, TASK_FS}:
        raise ValueError("unsupported pooled task selection")
    eligible, group_to_fold, inner_assignments = {}, {}, {}
    for task in sorted(tasks):
        samples, _ = load_samples(root, task, example=example)
        for sid, sample in samples.items():
            label = sample.get("label")
            mask = sample.get("mask", "1")
            if task == TASK_KNOT:
                if label not in {"0", "1"} or mask != "1":
                    raise ValueError(f"invalid knot label or mask: {sid}")
                training_label = label
                biological_target = label
            else:
                training_label = sample.get("training_label")
                if (label, mask, training_label) not in {
                    ("1", "1", "1"),
                    ("NA", "0", "0"),
                }:
                    raise ValueError(
                        f"invalid fold-switch biological label or mask: {sid}"
                    )
                if len(sample["sequence_ids"].split(";")) != (2 if label == "1" else 1):
                    raise ValueError(f"fold-switch endpoint count mismatch: {sid}")
                biological_target = "confirmed_positive" if label == "1" else "null"
            sample.update(
                label=label,
                mask=mask,
                learning_label=training_label,
                learning_mask="1",
                biological_target=biological_target,
                biological_mask=mask,
            )
            eligible[task, sid] = sample
        n_outer = 5 if task == TASK_KNOT else 10
        folds, inner = split_assignments(task_directory(root, task), samples, n_outer)
        inner_assignments[task] = inner
        for sid, sample in samples.items():
            group_to_fold[task, sample["binding_group_id"]] = str(folds[sid])
        label_key = "label" if task == TASK_KNOT else "learning_label"
        for outer in range(n_outer):
            test = [s for sid, s in samples.items() if folds[sid] == outer]
            train = [s for sid, s in samples.items() if folds[sid] != outer]
            if {s[label_key] for s in test} != {"0", "1"} or {
                s[label_key] for s in train
            } != {"0", "1"}:
                raise ValueError("outer partition lacks both operational classes")
            if task == TASK_FS and sum(s[label_key] == "1" for s in test) != 1:
                raise ValueError(
                    "fold-switch test fold must contain one confirmed case"
                )
            for validation in range(3):
                for held in (False, True):
                    labels = {
                        s[label_key]
                        for s in train
                        if (inner[outer][s["binding_group_id"]] == validation) == held
                    }
                    if labels != {"0", "1"}:
                        raise ValueError(
                            "inner partition lacks both operational classes"
                        )
    return eligible, group_to_fold, inner_assignments


def load_disorder(root: Path, example=False) -> tuple[dict, dict, dict, dict]:
    directory = task_directory(root, TASK_DISORDER)
    samples, _ = load_samples(root, TASK_DISORDER, example=example)
    expected = {sid: {} for sid in samples}
    seen = set()
    label_path = residue_label_path(directory)
    for row in read_table(label_path):
        sid, position = row["sample_id"], int(row["position"])
        if (
            sid not in samples
            or (sid, position) in seen
            or not 1 <= position <= int(samples[sid]["length"])
        ):
            raise ValueError(
                f"unknown sample or duplicate/out-of-range label coordinate: {sid}/{position}"
            )
        seen.add((sid, position))
        label, mask = row["label"], row.get("mask", "1")
        if mask == "0":
            if label != "NA":
                raise ValueError("masked residue has a class label")
            continue
        if mask != "1" or label not in {"0", "1"}:
            raise ValueError("invalid residue label or mask")
        evidence = row["evidence"]
        allowed = (
            {"synthetic_example"}
            if example
            else (
                {"direct_experimental_disorder"}
                if label == "1"
                else {"direct_experimental_order", "pdb_residue_observed"}
            )
        )
        if evidence not in allowed:
            raise ValueError(f"unsupported residue evidence: {sid}/{position}")
        expected[sid][position] = int(label)
    for sid, positions in expected.items():
        if set(positions.values()) != {0, 1}:
            raise ValueError(
                f"residue task requires both classes for every protein: {sid}"
            )
    folds, inner = split_assignments(directory, samples, 5)
    return samples, expected, folds, inner
