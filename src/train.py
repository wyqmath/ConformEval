from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from .datasets import load_pooled, parse_fasta, sequence_path, write_tsv
from .evaluate_predictions import (
    TASK_FS,
    TASK_KNOT,
    evaluate,
    load_protocol,
    write_gzip_tsv,
)
from .evaluate_residue_predictions import (
    TASK as RESIDUE_TASK,
)
from .evaluate_residue_predictions import (
    calculate,
    load_task,
    read_scores,
)
from .model_interface import load_mean_representation, load_residue_representation
from .readers import fold_indices, fs_rank_percentiles, nearest_neighbor, sequence_aux
from .runtime import RunRecord, json_write, safe_path


def reader_config(path: Path):
    config = json.loads(path.read_text(encoding="utf-8-sig"))
    if (
        config.get("schema_version") != 1
        or config.get("preprocessing") != "fit_on_training_only"
    ):
        raise ValueError("reader configuration requires training-only preprocessing")
    predictors = config.get("predictors")
    grid = config.get("C_grid")
    if (
        not isinstance(predictors, list)
        or not predictors
        or len(set(predictors)) != len(predictors)
        or not set(predictors) <= {"logistic", "ridge"}
    ):
        raise ValueError("predictors must contain distinct logistic/ridge names")
    if (
        not isinstance(grid, list)
        or not grid
        or any(
            type(c) not in {int, float} or not np.isfinite(c) or c <= 0 for c in grid
        )
    ):
        raise ValueError("C_grid must contain finite positive values")
    if (
        type(config.get("seed")) is not int
        or type(config.get("max_iter")) is not int
        or config["max_iter"] < 1
        or config.get("class_weight") != "balanced"
    ):
        raise ValueError("invalid reader seed, iteration count or class weighting")
    if not set(config.get("baselines", [])) <= {
        "length",
        "composition",
        "nearest",
        "training-prior",
    }:
        raise ValueError("unsupported pooled baseline")
    if not set(config.get("residue_baselines", [])) <= {"traditional"}:
        raise ValueError("unsupported residue baseline")
    tolerance = config.get("metric_tie_tolerance", 1e-12)
    if (
        type(tolerance) not in {int, float}
        or not np.isfinite(tolerance)
        or tolerance < 0
    ):
        raise ValueError("metric_tie_tolerance must be finite and nonnegative")
    return config


def fit_reader(X, y, train, held, predictor, C, config, context, warning_rows):
    if set(y[train]) != {0, 1}:
        raise ValueError(f"training partition lacks both classes: {context}")
    scaler = StandardScaler()
    train_X = scaler.fit_transform(X[train])
    held_X = scaler.transform(X[held])
    if predictor == "logistic":
        model = LogisticRegression(
            C=C,
            solver="liblinear",
            class_weight="balanced",
            max_iter=config["max_iter"],
            random_state=config["seed"],
        )
    else:
        model = RidgeClassifier(
            alpha=1.0 / C,
            class_weight="balanced",
            solver="lsqr",
            max_iter=config["max_iter"],
        )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model.fit(train_X, y[train])
    for w in caught:
        warning_rows.append(
            {**context, "type": w.category.__name__, "message": str(w.message)}
        )
        if issubclass(w.category, ConvergenceWarning):
            raise ValueError("reader did not converge")
    score = (
        model.predict_proba(held_X)[:, 1]
        if predictor == "logistic"
        else model.decision_function(held_X)
    )
    if not np.isfinite(score).all():
        raise ValueError("non-finite reader predictions")
    params = {
        "scaler_mean": scaler.mean_,
        "scaler_scale": scaler.scale_,
        "coef": model.coef_,
        "intercept": model.intercept_,
        "classes": model.classes_,
        "C": np.asarray(C),
        "predictor": np.asarray(predictor),
    }
    return score, params


def load_archives(path, required, residue=False):
    manifest = json.loads(path.read_text(encoding="utf-8-sig"))
    layers = manifest.get("candidate_layers")
    if (
        not isinstance(layers, list)
        or not layers
        or len(set(layers)) != len(layers)
        or any(type(l) is not int or l < 1 for l in layers)
    ):
        raise ValueError("invalid candidate layer manifest")
    for field in ["model_id", "model_revision", "pooling"]:
        if not isinstance(manifest.get(field), str) or not manifest[field]:
            raise ValueError(f"missing representation metadata: {field}")
    if manifest["pooling"] not in {
        "mean",
        "max",
        "mean_over_residues_without_special_tokens",
    }:
        raise ValueError("unsupported declared pooling")
    entries = {}
    for row in manifest["archives"]:
        sequence_id = row["sequence_id"]
        if sequence_id in entries:
            raise ValueError("duplicate representation sequence identity")
        entries[sequence_id] = row
    missing = set(required) - set(entries)
    if missing:
        raise ValueError(
            f"missing representation sequence coverage: {sorted(missing)[:5]}"
        )
    arrays, archive_paths, dimensions = {}, [], set()
    for sequence_id, length in required.items():
        row = entries[sequence_id]
        if type(row["length"]) is not int or row["length"] != length:
            raise ValueError("representation sequence length mismatch")
        archive = safe_path(path.parent, row["file"])
        with np.load(archive, allow_pickle=False) as z:
            meta = json.loads(str(z["meta"].item()))
        if "pooling" in meta and meta["pooling"] != manifest["pooling"]:
            raise ValueError("archive pooling differs from extraction manifest")
        if (
            "extraction_config" in manifest
            and meta.get("extraction_config") != manifest["extraction_config"]
        ):
            raise ValueError("archive extraction configuration mismatch")
        if residue:
            array, positions, dimension = load_residue_representation(
                archive,
                sequence_id,
                length,
                tuple(layers),
                manifest["model_id"],
                manifest["model_revision"],
            )
            arrays[sequence_id] = (array, {int(p): i for i, p in enumerate(positions)})
        else:
            array, dimension = load_mean_representation(
                archive,
                sequence_id,
                length,
                tuple(layers),
                manifest["model_id"],
                manifest["model_revision"],
            )
            arrays[sequence_id] = array
        dimensions.add(dimension)
        archive_paths.append(str(archive))
    if len(dimensions) != 1 or (
        "output_dimension" in manifest
        and manifest["output_dimension"] != next(iter(dimensions))
    ):
        raise ValueError("representation dimensions differ across archives")
    return manifest, arrays, archive_paths


def pooled_data(data_dir, task, manifest_path, example=False):
    eligible, group_to_fold = load_pooled(data_dir, {task})
    sequences = parse_fasta(sequence_path(data_dir, task))
    rows = sorted(eligible.values(), key=lambda r: r["sample_id"])
    sample_sequences = [row["sequence_ids"].split(";") for row in rows]
    required = {sid: len(sequences[sid]) for ids in sample_sequences for sid in ids}
    manifest, arrays, archive_paths = load_archives(manifest_path, required)
    X = {
        layer: np.vstack(
            [
                np.mean([arrays[sid][i] for sid in ids], axis=0)
                for ids in sample_sequences
            ]
        )
        for i, layer in enumerate(manifest["candidate_layers"])
    }
    auxiliary = np.vstack(
        [
            np.mean([sequence_aux(sequences[sid]) for sid in ids], axis=0)
            for ids in sample_sequences
        ]
    )
    y, plan = fold_indices(task, rows, data_dir, group_to_fold)
    return dict(
        task=task,
        rows=rows,
        y=y,
        plan=plan,
        X=X,
        auxiliary=auxiliary,
        manifest=manifest,
        archive_paths=archive_paths,
        group_to_fold=group_to_fold,
        residue=False,
    )


def residue_data(data_dir, manifest_path, example):
    samples, expected, folds, inner = load_task(data_dir, example=example)
    rows = sorted(samples.values(), key=lambda r: r["sample_id"])
    required = {s["sequence_id"]: int(s["length"]) for s in rows}
    manifest, arrays, archive_paths = load_archives(
        manifest_path, required, residue=True
    )
    sequences = parse_fasta(sequence_path(data_dir, RESIDUE_TASK))
    X = {l: [] for l in manifest["candidate_layers"]}
    ys = []
    sample_index = []
    coordinates = []
    traditional = []
    aa = "ACDEFGHIKLMNPQRSTVWY"
    offsets = {}
    offset = 0
    for i, s in enumerate(rows):
        sid = s["sample_id"]
        positions = sorted(expected[sid])
        sequence_id = s["sequence_id"]
        array, index = arrays[sequence_id]
        if set(positions) - set(index):
            raise ValueError("missing labeled residue coordinates in representation")
        for j, layer in enumerate(manifest["candidate_layers"]):
            X[layer].append(array[j, [index[p] for p in positions]])
        ys.extend(expected[sid][p] for p in positions)
        sample_index.extend([i] * len(positions))
        coordinates.extend(positions)
        offsets[sid] = np.arange(offset, offset + len(positions), dtype=np.int64)
        offset += len(positions)
        seq = sequences[s["sequence_id"]]
        onehot = np.asarray(
            [[float(c == a) for a in aa] for c in seq], dtype=np.float32
        )
        # Compute each context from residues only; no dataset-fitted preprocessing.
        features = np.column_stack(
            [
                onehot,
                np.arange(len(seq), dtype=np.float32) / max(1, len(seq) - 1),
                np.vstack(
                    [
                        onehot[max(0, p - 7) : min(len(seq), p + 8)].mean(axis=0)
                        for p in range(len(seq))
                    ]
                ),
            ]
        )
        traditional.append(features[np.asarray(positions) - 1])
    X = {l: np.concatenate(v) for l, v in X.items()}

    def ix(ids):
        return np.concatenate([offsets[sid] for sid in sorted(ids)])

    plan = {}
    for outer in range(5):
        train_ids = {sid for sid in samples if folds[sid] != outer}
        test_ids = set(samples) - train_ids
        inners = []
        for k in range(3):
            validation = {
                sid
                for sid in train_ids
                if inner[outer][samples[sid]["binding_group_id"]] == k
            }
            inners.append((str(k), ix(train_ids - validation), ix(validation)))
        plan[outer] = (ix(train_ids), ix(test_ids), inners)
    return dict(
        task=RESIDUE_TASK,
        rows=rows,
        y=np.asarray(ys, dtype=np.int8),
        plan=plan,
        X=X,
        auxiliary=np.concatenate(traditional),
        manifest=manifest,
        archive_paths=archive_paths,
        residue=True,
        samples=samples,
        expected=expected,
        folds=folds,
        sample_index=np.asarray(sample_index),
        positions=np.asarray(coordinates),
    )


def selection_value(data, held, scores):
    y = data["y"][held]
    if data["residue"]:
        identity = data["sample_index"][held]
        values = [
            average_precision_score(y[identity == s], scores[identity == s])
            for s in np.unique(identity)
        ]
        return float(np.mean(values))
    if data["task"] == TASK_KNOT:
        return float(roc_auc_score(y, scores))
    return float(np.mean(fs_rank_percentiles(scores, y)))


def select_candidates(
    data, features, inners, config, outer, selection_rows, warning_rows
):
    values = []
    for layer, X in features.items():
        for predictor in config["predictors"]:
            for C in sorted(set(config["C_grid"])):
                inner_metrics = []
                inner_weights = []
                for inner, train, validation in inners:
                    context = {
                        "outer_fold": outer,
                        "inner_fold": inner,
                        "layer": str(layer),
                        "predictor": predictor,
                        "C": C,
                    }
                    scores, _ = fit_reader(
                        X,
                        data["y"],
                        train,
                        validation,
                        predictor,
                        C,
                        config,
                        context,
                        warning_rows,
                    )
                    metric = selection_value(data, validation, scores)
                    inner_metrics.append(metric)
                    proteins = (
                        len(np.unique(data["sample_index"][validation]))
                        if data["residue"]
                        else 1
                    )
                    inner_weights.append(proteins)
                    selection_rows.append(
                        {
                            **context,
                            "validation_metric": metric,
                            "validation_rows": len(validation),
                            "aggregation_weight": proteins,
                        }
                    )
                value = float(np.average(inner_metrics, weights=inner_weights))
                values.append((layer, predictor, C, value))
    return values


def choose_candidate(candidates, task, tolerance):
    """Anchor numerical ties to the exact best metric, avoiding chained ties."""
    best = (min if task == TASK_FS else max)(c[3] for c in candidates)
    tied = [c for c in candidates if np.isclose(c[3], best, rtol=0, atol=tolerance)]
    return min(
        tied, key=lambda c: (c[2], int(c[0]) if isinstance(c[0], int) else 0, c[1])
    )


def train_task(args, record, config):
    data_dir = args.data_dir.resolve()
    task = args.task
    data = (
        residue_data(data_dir, args.feature_manifest, args.example)
        if task == RESIDUE_TASK
        else pooled_data(data_dir, task, args.feature_manifest, args.example)
    )
    manifest = data["manifest"]
    layers = config.get("candidate_layers", manifest["candidate_layers"])
    if (
        not layers
        or len(set(layers)) != len(layers)
        or not set(layers) <= set(data["X"])
    ):
        raise ValueError(
            "configured candidate layers are absent from representation manifest"
        )
    data["X"] = {l: data["X"][l] for l in layers}
    for archive in data["archive_paths"]:
        record.track(Path(archive))
    for p in data_dir.rglob("*"):
        if p.is_file() and p.suffix != ".md":
            record.track(p)
    arm_scores = {}
    selection_rows = []
    outer_rows = []
    warning_rows = record.data["warnings"]
    parameters = record.output / "parameters"
    parameters.mkdir()
    for outer, (train, test, inners) in sorted(data["plan"].items()):
        candidates = select_candidates(
            data, data["X"], inners, config, outer, selection_rows, warning_rows
        )
        tolerance = config.get("metric_tie_tolerance", 1e-12)
        choices = {"PLM-nested": choose_candidate(candidates, task, tolerance)}
        if config.get("diagnostic_layers", False):
            choices.update(
                {
                    f"PLM-hidden{l}": choose_candidate(
                        [c for c in candidates if c[0] == l], task, tolerance
                    )
                    for l in layers
                }
            )
        features = dict(data["X"])
        if data["residue"]:
            if "traditional" in config.get("residue_baselines", []):
                features["traditional-AA-context"] = data["auxiliary"]
        else:
            for name, values in [
                ("length", data["auxiliary"][:, :1]),
                ("composition", data["auxiliary"][:, 1:]),
            ]:
                if name in config.get("baselines", []):
                    features["PLM-" + name] = values
        for name, X in features.items():
            if isinstance(name, int):
                continue
            baseline = select_candidates(
                data, {name: X}, inners, config, outer, selection_rows, warning_rows
            )
            choices[name] = choose_candidate(baseline, task, tolerance)
        for arm, (layer, predictor, C, metric) in choices.items():
            score, params = fit_reader(
                features[layer],
                data["y"],
                train,
                test,
                predictor,
                C,
                config,
                {"outer_fold": outer, "phase": "outer-refit", "arm": arm},
                warning_rows,
            )
            arm_scores.setdefault(arm, np.full(len(data["y"]), np.nan))[test] = score
            np.savez_compressed(
                parameters / f"{arm}_outer{outer:02d}.npz",
                **params,
                train_rows=train,
                held_rows=test,
                layer=np.asarray(str(layer)),
            )
            outer_rows.append(
                {
                    "outer_fold": outer,
                    "arm": arm,
                    "layer": layer,
                    "predictor": predictor,
                    "C": C,
                    "inner_metric": metric,
                    "train_rows": len(train),
                    "test_rows": len(test),
                }
            )
        if not data["residue"]:
            for name in config.get("baselines", []):
                if name == "training-prior":
                    score = np.full(len(test), data["y"][train].mean())
                    json_write(
                        parameters / f"PLM-training-prior_outer{outer:02d}.json",
                        {
                            "training_prior": float(score[0]),
                            "train_sample_ids": [
                                data["rows"][i]["sample_id"] for i in train
                            ],
                        },
                    )
                elif name == "nearest":
                    score, reference = nearest_neighbor(
                        data["X"][max(layers)],
                        data["y"],
                        train,
                        test,
                        [r["sample_id"] for r in data["rows"]],
                    )
                    reference["layer"] = max(layers)
                    json_write(
                        parameters / f"PLM-nearest_outer{outer:02d}.json", reference
                    )
                else:
                    continue
                arm_scores.setdefault("PLM-" + name, np.full(len(data["y"]), np.nan))[
                    test
                ] = score
        record.data["completed_outer_folds"] = int(outer) + 1
        record.save()
        print(
            f"Completed {task} outer fold {int(outer) + 1}/{len(data['plan'])}",
            flush=True,
        )
    if any(not np.isfinite(v).all() for v in arm_scores.values()):
        raise ValueError("outer predictions incomplete or non-finite")
    metadata = {
        "project": "ConformEval",
        "model_id": manifest["model_id"],
        "model_revision": manifest["model_revision"],
        "prediction_method": "nested grouped configured readers",
        "training_protocol": "nested_grouped_oof",
        "candidate_layers": list(layers),
        "pooling": manifest["pooling"],
        "arms": sorted(arm_scores),
        "extraction_config": manifest.get("extraction_config", {}),
        "reader_configuration": config,
        "scope": "software_example"
        if args.example
        else "retrospective_grouped_development",
    }
    json_write(record.output / "model_metadata.json", metadata)
    write_tsv(
        record.output / "inner_selection.tsv", list(selection_rows[0]), selection_rows
    )
    write_tsv(record.output / "outer_selections.tsv", list(outer_rows[0]), outer_rows)
    pred_rows = []
    for arm, scores in sorted(arm_scores.items()):
        for outer, (_, test, _) in sorted(data["plan"].items()):
            for index in test:
                sample = data["rows"][
                    data["sample_index"][index] if data["residue"] else index
                ]
                row = {
                    "task_id": task,
                    "sample_id": sample["sample_id"],
                    "outer_fold": outer,
                    "score": float(scores[index]),
                    "arm": arm,
                }
                if data["residue"]:
                    row["position_1based"] = int(data["positions"][index])
                pred_rows.append(row)
    prediction_path = record.output / "predictions.tsv.gz"
    write_gzip_tsv(prediction_path, pred_rows, list(pred_rows[0]))
    if data["residue"]:
        checked = read_scores(
            prediction_path,
            data["samples"],
            data["expected"],
            data["folds"],
            set(arm_scores),
            manifest["model_id"],
        )
        protein_rows, metrics = calculate(
            data["samples"], data["expected"], data["folds"], checked
        )
        write_tsv(
            record.output / "protein_metrics.tsv", list(protein_rows[0]), protein_rows
        )
        json_write(record.output / "metrics.json", metrics)
    else:
        metrics = {}
        for arm in sorted(arm_scores):
            arm_path = record.output / (arm + ".tsv.gz")
            write_gzip_tsv(
                arm_path,
                [
                    {k: v for k, v in r.items() if k != "arm"}
                    for r in pred_rows
                    if r["arm"] == arm
                ],
                ["task_id", "sample_id", "outer_fold", "score"],
            )
            evaluate(
                data_dir,
                arm_path,
                record.output / "model_metadata.json",
                args.protocol,
                record.output / "scores" / arm,
                {task},
                example=args.example,
            )
            metrics[arm] = json.loads(
                (record.output / "scores" / arm / "metrics.json").read_text(
                    encoding="utf-8-sig"
                )
            )["metrics"][task]
        json_write(record.output / "metrics.json", metrics)
    record.data.update(
        task_id=task,
        scope=metadata["scope"],
        candidate_layers=list(layers),
        model_id=manifest["model_id"],
        model_revision=manifest["model_revision"],
        eligible_samples=len(data["rows"]),
        scored_rows_per_arm=len(data["y"]),
        arms=sorted(arm_scores),
        outer_folds=len(data["plan"]),
        inner_folds=3,
        outer_labels_used_for_selection=False,
    )
    return data


def run(args):
    output = args.output_dir.resolve()
    record = RunRecord(output, {"task": getattr(args, "task", None)})
    try:
        if not args.data_dir.is_dir():
            raise FileNotFoundError(args.data_dir)
        config = reader_config(args.config)
        record.data["configuration"] = {
            "task": args.task,
            "reader": config,
            "example": args.example,
        }
        for p in [args.config, args.feature_manifest]:
            record.track(p)
        if args.protocol is None:
            args.protocol = Path("configs/scoring.json")
        load_protocol(args.protocol)
        record.track(args.protocol)
        train_task(args, record, config)
    except BaseException as error:
        record.finish(error)
        raise
    error = record.finish()
    if error:
        raise error
    print(
        json.dumps(
            {
                "status": "PASS",
                "task": args.task,
                "output_dir": str(output),
                "arms": record.data["arms"],
            }
        )
    )
    return record.data


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--task", required=True, choices=[TASK_KNOT, TASK_FS, RESIDUE_TASK])
    p.add_argument("--data-dir", required=True, type=Path)
    p.add_argument("--feature-manifest", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--protocol", type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--example", action="store_true")
    run(p.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
