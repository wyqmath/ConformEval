from __future__ import annotations

import math
from collections import Counter

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler

from .datasets import TASK_KNOT
from .scoring import midrank

AA20 = "ACDEFGHIKLMNPQRSTVWY"


def sequence_aux(sequence: str) -> np.ndarray:
    length = len(sequence)
    counts = Counter(sequence)
    return np.asarray(
        [math.log1p(length)] + [counts[aa] / length for aa in AA20], dtype=np.float32
    )


def fold_indices(task: str, rows: list[dict], group_to_fold: dict, assignments: dict):
    n_outer = 5 if task == TASK_KNOT else 10
    inner_group = {
        str(outer): {group: str(inner) for group, inner in groups.items()}
        for outer, groups in assignments.items()
    }
    y = np.asarray(
        [
            int(row["label"] if task == TASK_KNOT else row["learning_label"])
            for row in rows
        ],
        dtype=np.int8,
    )
    groups = {row["binding_group_id"] for row in rows}
    plan = {}
    for outer_index in range(n_outer):
        outer = str(outer_index)
        test_groups = {group for group in groups if group_to_fold[task, group] == outer}
        train_groups = groups - test_groups
        inner_assignment = inner_group.get(outer, {})
        if set(inner_assignment) != train_groups:
            raise ValueError(
                f"inner validation groups do not partition outer training: {task}/{outer}"
            )
        train = np.asarray(
            [
                i
                for i, row in enumerate(rows)
                if row["binding_group_id"] in train_groups
            ],
            dtype=np.int64,
        )
        test = np.asarray(
            [i for i, row in enumerate(rows) if row["binding_group_id"] in test_groups],
            dtype=np.int64,
        )
        if set(y[train]) != {0, 1} or set(y[test]) != {0, 1}:
            raise ValueError(
                f"outer split lacks both operational classes: {task}/{outer}"
            )
        inners = []
        for inner_index in range(3):
            inner = str(inner_index)
            val_groups = {
                group for group, fold in inner_assignment.items() if fold == inner
            }
            val = np.asarray(
                [
                    i
                    for i, row in enumerate(rows)
                    if row["binding_group_id"] in val_groups
                ],
                dtype=np.int64,
            )
            inner_train = np.asarray(
                [
                    i
                    for i, row in enumerate(rows)
                    if row["binding_group_id"] in train_groups - val_groups
                ],
                dtype=np.int64,
            )
            if not len(val) or set(y[val]) != {0, 1} or set(y[inner_train]) != {0, 1}:
                raise ValueError(
                    f"inner split lacks both operational classes: {task}/{outer}/{inner}"
                )
            inners.append((inner_index, inner_train, val))
        plan[outer_index] = (train, test, inners)
    return (y, plan)


def fs_rank_percentiles(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    positives = np.flatnonzero(labels == 1)
    if not len(positives):
        raise ValueError("PU validation fold contains no held positive")
    output = []
    for index in positives:
        rank, _, _ = midrank(scores, int(index))
        output.append(rank / len(scores))
    return np.asarray(output, dtype=np.float64)


def nearest_neighbor(
    X: np.ndarray,
    y: np.ndarray,
    train: np.ndarray,
    test: np.ndarray,
    sample_ids: list[str],
):
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(X[train])
    test_scaled = scaler.transform(X[test])
    order = np.argsort(
        np.asarray([sample_ids[i] for i in train], dtype=str), kind="stable"
    )
    train_sorted = train[order]
    train_sorted_features = train_scaled[order]
    distances = cdist(test_scaled, train_sorted_features, metric="euclidean")
    if not np.isfinite(distances).all():
        raise ValueError("nonfinite nearest-neighbor distance")
    nearest = np.argmin(distances, axis=1)
    return (
        y[train_sorted[nearest]].astype(np.float64),
        {
            "reference_sample_ids": [sample_ids[i] for i in train_sorted],
            "reference_labels": y[train_sorted].astype(int).tolist(),
            "scaler_mean": scaler.mean_.tolist(),
            "scaler_scale": scaler.scale_.tolist(),
        },
    )
