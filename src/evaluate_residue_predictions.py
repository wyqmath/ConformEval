from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

from .datasets import (
    TASK_DISORDER as TASK,
)
from .datasets import (
    load_disorder as load_task,
)
from .datasets import (
    task_input_paths,
    write_tsv,
)
from .runtime import RunRecord, json_write
from .scoring import (
    BOOTSTRAP_SEED,
    BOOTSTRAPS,
    CI_QUANTILES,
    load_metadata,
    load_protocol,
)


def read_scores(path, samples, expected, folds, arms, default_arm):
    scores = {arm: {sid: {} for sid in samples} for arm in arms}
    opener = gzip.open if path.suffix == ".gz" else open
    required = {"task_id", "sample_id", "position_1based", "outer_fold", "score"}
    with opener(path, "rt", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("missing prediction columns")
        for r in reader:
            sid = r["sample_id"]
            arm = r.get("arm", default_arm)
            if r["task_id"] != TASK or sid not in samples or arm not in arms:
                raise ValueError("unexpected task/sample/arm")
            pos = int(r["position_1based"])
            if pos not in expected[sid]:
                raise ValueError(
                    "prediction coordinate is unknown, masked, or out of range"
                )
            if int(r["outer_fold"]) != folds[sid]:
                raise ValueError("prediction is outside its frozen outer test fold")
            if (
                "binding_group_id" in r
                and r["binding_group_id"] != samples[sid]["binding_group_id"]
            ):
                raise ValueError("prediction group mismatch")
            if "mask" in r and r["mask"] != "1":
                raise ValueError("masked prediction included")
            if "label" in r and r["label"] != str(expected[sid][pos]):
                raise ValueError("prediction label mismatch")
            score = float(r["score"])
            if not math.isfinite(score):
                raise ValueError("non-finite prediction")
            if pos in scores[arm][sid]:
                raise ValueError("duplicate prediction identity")
            scores[arm][sid][pos] = score
    for arm in arms:
        for sid, values in scores[arm].items():
            if set(values) != set(expected[sid]):
                raise ValueError(f"incomplete prediction coordinates: {arm}/{sid}")
    return scores


def calculate(
    samples, expected, folds, scores, replicates=BOOTSTRAPS, seed=BOOTSTRAP_SEED
):
    rows = []
    summaries = {}
    arms = sorted(scores)
    for arm in arms:
        pooled_scores = []
        pooled_labels = []
        for sid in sorted(samples):
            pos = sorted(expected[sid])
            y = np.array([expected[sid][p] for p in pos], dtype=np.int8)
            s = np.array([scores[arm][sid][p] for p in pos], dtype=np.float64)
            ap, auc = float(average_precision_score(y, s)), float(roc_auc_score(y, s))
            rows.append(
                {
                    "task_id": TASK,
                    "sample_id": sid,
                    "binding_group_id": samples[sid]["binding_group_id"],
                    "outer_fold": folds[sid],
                    "arm": arm,
                    "n_labeled_residues": len(y),
                    "n_positive": int(y.sum()),
                    "n_negative": int(len(y) - y.sum()),
                    "average_precision": ap,
                    "auroc": auc,
                }
            )
            pooled_scores.append(s)
            pooled_labels.append(y)
        rs = [r for r in rows if r["arm"] == arm]
        summaries[arm] = {
            "n_proteins": len(rs),
            "mean_per_protein_AP": float(np.mean([r["average_precision"] for r in rs])),
            "sd_per_protein_AP": float(
                np.std([r["average_precision"] for r in rs], ddof=0)
            ),
            "mean_per_protein_AUROC": float(np.mean([r["auroc"] for r in rs])),
            "pooled_residue_AP_secondary": float(
                average_precision_score(
                    np.concatenate(pooled_labels), np.concatenate(pooled_scores)
                )
            ),
        }
    groups = sorted({r["binding_group_id"] for r in rows})
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(groups), size=(replicates, len(groups)))
    boot = {}
    for arm in arms:
        by_group = {
            g: [
                r["average_precision"]
                for r in rows
                if r["arm"] == arm and r["binding_group_id"] == g
            ]
            for g in groups
        }
        totals = np.array([sum(by_group[g]) for g in groups])
        counts = np.array([len(by_group[g]) for g in groups])
        means = totals[draws].sum(axis=1) / counts[draws].sum(axis=1)
        boot[arm] = {
            "macro_ap_ci95": np.quantile(means, CI_QUANTILES).tolist(),
            "replicates": replicates,
            "bootstrap_unit": "binding_group_id",
        }
    summaries["group_bootstrap"] = boot
    return (rows, summaries)


def evaluate(
    task_dir, predictions_path, metadata_path, protocol_path, output_dir, example=False
):
    with RunRecord(output_dir, {"task": TASK, "example": example}) as record:
        protocol = load_protocol(protocol_path)
        record.data["configuration"]["scoring"] = protocol
        for path in [predictions_path, metadata_path, protocol_path] + task_input_paths(
            task_dir, TASK
        ):
            record.track(path)
        metadata = load_metadata(metadata_path, example=example)
        arms = metadata.get("arms", [metadata["model_id"]])
        if (
            not isinstance(arms, list)
            or not arms
            or any(not isinstance(arm, str) or not arm for arm in arms)
            or len(set(arms)) != len(arms)
        ):
            raise ValueError("invalid declared arms")
        metadata["arms"] = arms
        samples, expected, folds, _ = load_task(task_dir, example=example)
        scores = read_scores(
            predictions_path, samples, expected, folds, set(arms), metadata["model_id"]
        )
        uncertainty = protocol["uncertainty"]
        rows, metrics = calculate(
            samples,
            expected,
            folds,
            scores,
            replicates=uncertainty["bootstrap_replicates"],
            seed=uncertainty["seed"],
        )
        write_tsv(output_dir / "protein_metrics.tsv", list(rows[0]), rows)
        json_write(output_dir / "metrics.json", metrics)
        record.data.update(
            model_metadata=metadata,
            scoring_protocol_id=protocol["protocol_id"],
            scope="software_example"
            if example
            else "retrospective_grouped_development",
            n_predictions=sum(
                len(positions) for arm in scores.values() for positions in arm.values()
            ),
        )
    return record.data


def main():
    parser = argparse.ArgumentParser(
        description="Score out-of-fold residue disorder predictions."
    )
    parser.add_argument("--task-dir", required=True, type=Path)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--model-metadata", required=True, type=Path)
    parser.add_argument("--protocol", type=Path, default=Path("configs/scoring.json"))
    parser.add_argument(
        "--example", action="store_true", help="Score synthetic example data"
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    result = evaluate(
        args.task_dir.resolve(),
        args.predictions.resolve(),
        args.model_metadata.resolve(),
        args.protocol.resolve(),
        args.output_dir.resolve(),
        example=args.example,
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "n_predictions": result["n_predictions"],
                "output_dir": str(args.output_dir),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
