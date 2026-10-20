from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from .datasets import TASK_FS, TASK_KNOT, load_pooled, task_input_paths, write_tsv
from .runtime import RunRecord, json_write
from .scoring import (
    BOOTSTRAP_SEED,
    BOOTSTRAPS,
    CI_QUANTILES,
    TASKS,
    load_metadata,
    load_protocol,
    midrank,
)

PREDICTION_FIELDS = ("task_id", "sample_id", "outer_fold", "score")


def read_predictions(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or ()) != set(PREDICTION_FIELDS):
            raise ValueError(
                f"prediction table must contain exactly these fields: {', '.join(PREDICTION_FIELDS)}"
            )
        return list(reader)


def load_predictions(
    path: Path, eligible: dict, group_to_fold: dict
) -> dict[tuple[str, str], float]:
    seen: dict[tuple[str, str], float] = {}
    for line_number, row in enumerate(read_predictions(path), start=2):
        key = (row.get("task_id", ""), row.get("sample_id", ""))
        if key[0] not in TASKS:
            raise ValueError(
                f"unknown or unsupported task at prediction line {line_number}: {key[0]}"
            )
        if key not in eligible:
            raise ValueError(f"prediction sample is absent or ineligible: {key}")
        if key in seen:
            raise ValueError(f"duplicate prediction identity: {key}")
        expected_fold = group_to_fold.get((key[0], eligible[key]["binding_group_id"]))
        if row.get("outer_fold") != expected_fold:
            raise ValueError(
                f"prediction fold differs from frozen outer assignment: {key}"
            )
        try:
            score = float(row["score"])
        except (KeyError, ValueError) as exc:
            raise ValueError(f"prediction score is missing or invalid: {key}") from exc
        if not math.isfinite(score):
            raise ValueError(f"prediction score is non-finite: {key}")
        seen[key] = score
    if set(seen) != set(eligible):
        missing = sorted(set(eligible) - set(seen))[:10]
        extra = sorted(set(seen) - set(eligible))[:10]
        raise ValueError(
            f"prediction coverage differs from eligible data; missing={missing}; extra={extra}"
        )
    return seen


def rank_positive(scores: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    if int(np.sum(labels == 1)) != 1:
        raise ValueError(
            "each FS-L2 frozen outer fold must contain exactly one positive case"
        )
    rank, greater, equal = midrank(scores, int(np.flatnonzero(labels == 1)[0]))
    recall = {
        f"recall_at_{k}": max(0.0, min(equal, k - greater)) / equal
        for k in (1, 10, 100, 1000, 10000)
    }
    return {
        "rank": float(rank),
        "universe": int(len(scores)),
        "percentile": rank / len(scores),
        **recall,
    }


def group_bootstrap_auroc(
    rows: list[dict], score_map: dict, replicates: int, seed: int
) -> dict:
    groups = sorted({row["binding_group_id"] for row in rows})
    members = {
        group: [row for row in rows if row["binding_group_id"] == group]
        for group in groups
    }
    rng = np.random.default_rng(seed)
    values = []
    for draw in rng.integers(0, len(groups), size=(replicates, len(groups))):
        sampled = [row for j in draw for row in members[groups[int(j)]]]
        labels = np.asarray([int(row["label"]) for row in sampled], dtype=np.int8)
        if len(np.unique(labels)) < 2:
            continue
        scores = np.asarray(
            [score_map[row["task_id"], row["sample_id"]] for row in sampled],
            dtype=np.float64,
        )
        values.append(float(roc_auc_score(labels, scores)))
    if not values:
        raise ValueError("group bootstrap produced no valid KNOT-PRESENCE replicates")
    return {
        "unit": "binding_group_id",
        "replicates": replicates,
        "valid_replicates": len(values),
        "seed": seed,
        "conditional_on_fixed_outer_predictions": True,
        "ci95": np.quantile(values, CI_QUANTILES).tolist(),
    }


def fs_metrics(
    rows: list[dict], score_map: dict, group_to_fold: dict, replicates: int, seed: int
) -> dict:
    by_fold: dict[str, list[dict]] = {}
    for row in rows:
        fold = group_to_fold[TASK_FS, row["binding_group_id"]]
        by_fold.setdefault(fold, []).append(row)
    ranks = []
    for fold, fold_rows in sorted(by_fold.items(), key=lambda item: int(item[0])):
        ordered = sorted(fold_rows, key=lambda row: row["sample_id"])
        scores = np.asarray(
            [score_map[TASK_FS, row["sample_id"]] for row in ordered], dtype=np.float64
        )
        labels = np.asarray(
            [int(row["learning_label"]) for row in ordered], dtype=np.int8
        )
        rank = rank_positive(scores, labels)
        positive = ordered[int(np.flatnonzero(labels == 1)[0])]
        ranks.append(
            {
                "outer_fold": fold,
                "sample_id": positive["sample_id"],
                "binding_group_id": positive["binding_group_id"],
                **rank,
            }
        )
    if len(ranks) != 10 or len({row["binding_group_id"] for row in ranks}) != 10:
        raise ValueError(
            "FS-L2 positive case coverage differs from its frozen ten outer units"
        )
    values = np.asarray([row["percentile"] for row in ranks], dtype=np.float64)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(ranks), size=(replicates, len(ranks)))
    means = values[draws].mean(axis=1)
    summary = {
        "primary_metric": "mean held-positive rank percentile within its own frozen outer-fold test pool; lower is better",
        "mean_held_case_rank_percentile": float(values.mean()),
        "mean_recall_at_1": float(np.mean([row["recall_at_1"] for row in ranks])),
        "mean_recall_at_10": float(np.mean([row["recall_at_10"] for row in ranks])),
        "mean_recall_at_100": float(np.mean([row["recall_at_100"] for row in ranks])),
        "mean_recall_at_1000": float(np.mean([row["recall_at_1000"] for row in ranks])),
        "mean_recall_at_10000": float(
            np.mean([row["recall_at_10000"] for row in ranks])
        ),
        "held_positive_cases": len(ranks),
        "unlabeled_background_is_not_a_verified_negative": True,
        "group_bootstrap_ci95": np.quantile(means, CI_QUANTILES).tolist(),
        "bootstrap_unit": "positive binding_group_id",
        "bootstrap_replicates": replicates,
        "bootstrap_seed": seed,
        "bootstrap_scope": "conditional on fixed outer-fold predictions; no refitting uncertainty claim",
        "case_ranks": ranks,
    }
    return summary


def evaluate(
    data_dir: Path,
    predictions_path: Path,
    metadata_path: Path,
    protocol_path: Path,
    output_dir: Path,
    tasks: set[str] | None = None,
    example: bool = False,
) -> dict:
    tasks = TASKS if tasks is None else set(tasks)
    with RunRecord(output_dir, {"tasks": sorted(tasks), "example": example}) as record:
        protocol = load_protocol(protocol_path)
        record.data["configuration"]["scoring"] = protocol
        for path in [predictions_path, metadata_path, protocol_path]:
            record.track(path)
        metadata = load_metadata(metadata_path, example=example, require_layers=True)
        eligible, group_to_fold, _ = load_pooled(data_dir, tasks, example=example)
        for task in sorted(tasks):
            for path in task_input_paths(data_dir, task):
                record.track(path)
        prediction_map = load_predictions(predictions_path, eligible, group_to_fold)
        scored_rows = []
        metrics: dict[str, dict] = {}
        for task in sorted(tasks):
            rows = sorted(
                (row for (row_task, _), row in eligible.items() if row_task == task),
                key=lambda row: row["sample_id"],
            )
            if task == TASK_KNOT:
                labels = np.asarray([int(row["label"]) for row in rows], dtype=np.int8)
                scores = np.asarray(
                    [prediction_map[task, row["sample_id"]] for row in rows],
                    dtype=np.float64,
                )
                summary = {
                    "primary_metric": "AUROC from complete frozen outer-fold out-of-fold predictions",
                    "auroc": float(roc_auc_score(labels, scores)),
                    "average_precision_secondary": float(
                        average_precision_score(labels, scores)
                    ),
                    "n_samples": int(len(rows)),
                    "n_positive": int(labels.sum()),
                    "n_binding_groups": len({row["binding_group_id"] for row in rows}),
                    "group_bootstrap": group_bootstrap_auroc(
                        rows, prediction_map, BOOTSTRAPS, BOOTSTRAP_SEED
                    ),
                }
            else:
                summary = fs_metrics(
                    rows, prediction_map, group_to_fold, BOOTSTRAPS, BOOTSTRAP_SEED
                )
                summary["n_evaluation_units"] = len(rows)
                summary["n_binding_groups"] = len(
                    {row["binding_group_id"] for row in rows}
                )
            metrics[task] = summary
            for row in rows:
                fold = group_to_fold[task, row["binding_group_id"]]
                scored_rows.append(
                    {
                        "task_id": task,
                        "sample_id": row["sample_id"],
                        "outer_fold": fold,
                        "binding_group_id": row["binding_group_id"],
                        "score": prediction_map[task, row["sample_id"]],
                        "learning_label": row["learning_label"],
                        "learning_mask": row["learning_mask"],
                        "biological_target": row["biological_target"],
                        "biological_mask": row["biological_mask"],
                        "label": row["label"],
                        "mask": row["mask"],
                    }
                )
        scored_path = output_dir / "scored_predictions.tsv.gz"
        write_tsv(
            scored_path,
            [
                "task_id",
                "sample_id",
                "outer_fold",
                "binding_group_id",
                "score",
                "learning_label",
                "learning_mask",
                "biological_target",
                "biological_mask",
                "label",
                "mask",
            ],
            scored_rows,
        )
        json_write(output_dir / "metrics.json", {"model": metadata, "metrics": metrics})
        record.data.update(
            scope="software_example"
            if example
            else "retrospective_grouped_development",
            dataset_id="ConformEval",
            scoring_protocol_id=protocol["protocol_id"],
            model_metadata=metadata,
            eligible_sample_counts={
                task: sum(key[0] == task for key in eligible) for task in sorted(tasks)
            },
            prediction_counts={
                task: sum(key[0] == task for key in prediction_map)
                for task in sorted(tasks)
            },
            outer_folds={
                task: len(
                    {
                        group_to_fold[task, row["binding_group_id"]]
                        for key, row in eligible.items()
                        if key[0] == task
                    }
                )
                for task in sorted(tasks)
            },
            scorer_uses_outer_test_labels_for_model_selection=False,
            model_training_protocol_self_attested=True,
        )
    return record.data


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Score complete out-of-fold predictions against the frozen ConformEval labels."
    )
    parser.add_argument(
        "--example",
        action="store_true",
        help="Use synthetic data; omit input paths to run the bundled knot example",
    )
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument(
        "--task",
        action="append",
        choices=sorted(TASKS),
        help="Repeat to select pooled tasks; default: both",
    )
    parser.add_argument("--protocol", type=Path)
    parser.add_argument(
        "--predictions",
        type=Path,
        help="TSV/TSV.GZ with task_id, sample_id, outer_fold, score",
    )
    parser.add_argument("--model-metadata", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    inputs = ("data_dir", "task", "protocol", "predictions", "model_metadata")
    if args.example and not any(getattr(args, name) is not None for name in inputs):
        args.data_dir = Path("examples/knot_presence")
        args.task = [TASK_KNOT]
        args.protocol = Path("configs/scoring.json")
        args.predictions = Path("examples/predictions/knot.tsv")
        args.model_metadata = Path("examples/predictions/metadata.json")
    missing = [
        "--" + name.replace("_", "-")
        for name in inputs
        if name != "task" and getattr(args, name) is None
    ]
    if missing:
        parser.error(
            "the following input arguments are required: " + ", ".join(missing)
        )
    try:
        result = evaluate(
            args.data_dir.resolve(),
            args.predictions.resolve(),
            args.model_metadata.resolve(),
            args.protocol.resolve(),
            args.output_dir.resolve(),
            set(args.task) if args.task else None,
            example=args.example,
        )
    except (ValueError, KeyError, OSError, TypeError, json.JSONDecodeError) as error:
        print(f"Evaluation rejected: {error}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": result["status"],
                "metrics": "metrics.json",
                "output_dir": str(args.output_dir),
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
