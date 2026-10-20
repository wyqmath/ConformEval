from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .datasets import TASK_DISORDER, TASK_FS, TASK_KNOT

TASKS = {TASK_KNOT, TASK_FS}
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 2026
CONFIDENCE_LEVEL = 0.95
CI_QUANTILES = (0.025, 0.975)


def load_protocol(path: Path) -> dict:
    protocol = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(protocol, dict):
        raise ValueError("scoring configuration must be a JSON object")
    if (
        type(protocol.get("schema_version")) is not int
        or protocol["schema_version"] != 1
        or protocol.get("protocol_id") != "ConformEval-grouped-evaluation"
    ):
        raise ValueError("unsupported scoring protocol")
    tasks = protocol.get("tasks", {})
    expected = {
        TASK_KNOT: (5, "pooled_outer_fold_OOF_AUROC", "mean_inner_validation_AUROC"),
        TASK_FS: (
            10,
            "mean_held_positive_rank_percentile_within_own_outer_test_pool_lower_is_better",
            "mean_inner_held_positive_rank_percentile_within_own_validation_pool_lower_is_better",
        ),
    }
    if not isinstance(tasks, dict) or set(tasks) != TASKS | {TASK_DISORDER}:
        raise ValueError(
            "scoring protocol task set differs from the released score-eligible tasks"
        )
    expected[TASK_DISORDER] = (
        5,
        "mean_per_protein_AP",
        "mean_inner_validation_per_protein_AP",
    )
    for task, (outer, primary, inner) in expected.items():
        record = tasks[task]
        if (
            not isinstance(record, dict)
            or type(record.get("outer_folds")) is not int
            or type(record.get("inner_folds")) is not int
            or record["outer_folds"] != outer
            or record.get("inner_folds") != 3
            or record.get("primary_metric") != primary
            or (record.get("inner_selection_metric") != inner)
        ):
            raise ValueError(
                f"task metric or fold protocol differs from the scoring configuration: {task}"
            )
    uncertainty = protocol.get("uncertainty", {})
    if (
        not isinstance(uncertainty, dict)
        or type(uncertainty.get("bootstrap_replicates")) is not int
        or type(uncertainty.get("seed")) is not int
        or uncertainty["bootstrap_replicates"] != BOOTSTRAPS
        or uncertainty.get("seed") != BOOTSTRAP_SEED
        or uncertainty.get("confidence_level") != CONFIDENCE_LEVEL
    ):
        raise ValueError("uncertainty protocol differs from the scoring configuration")
    return protocol


def load_metadata(path: Path, example=False, require_layers=False) -> dict:
    metadata = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(metadata, dict):
        raise ValueError("model metadata must be a JSON object")
    for key in ["model_id", "model_revision", "prediction_method", "training_protocol"]:
        if not isinstance(metadata.get(key), str) or not metadata[key].strip():
            raise ValueError(f"missing model metadata: {key}")
    if metadata["training_protocol"] != "nested_grouped_oof":
        raise ValueError("unsupported training protocol")
    scope = "software_example" if example else "retrospective_grouped_development"
    if "scope" in metadata and metadata["scope"] != scope:
        raise ValueError("model metadata scope differs from the selected data mode")
    if require_layers:
        layers = metadata.get("candidate_layers")
        if (
            not isinstance(layers, list)
            or not layers
            or any(type(layer) is not int or layer < 1 for layer in layers)
            or len(set(layers)) != len(layers)
        ):
            raise ValueError(
                "candidate_layers must contain unique positive 1-based layer IDs"
            )
        if (
            not isinstance(metadata.get("pooling"), str)
            or not metadata["pooling"].strip()
        ):
            raise ValueError("model metadata is missing pooling")
    return metadata


def midrank(scores: np.ndarray, index: int) -> tuple[float, int, int]:
    if scores.ndim != 1 or not len(scores) or not np.isfinite(scores).all():
        raise ValueError("rank scores must be a finite nonempty vector")
    if not 0 <= index < len(scores):
        raise ValueError("rank index is outside the score vector")
    greater = int(np.sum(scores > scores[index]))
    equal = int(np.sum(scores == scores[index]))
    return greater + (equal + 1) / 2, greater, equal
