from __future__ import annotations

import argparse
import json
from pathlib import Path

from .datasets import load_pooled, parse_fasta, sequence_path, task_input_paths
from .evaluate_residue_predictions import TASK as RESIDUE_TASK
from .evaluate_residue_predictions import load_task
from .runtime import RunRecord
from .scoring import TASKS


def required_sequences(data_dir: Path, task: str, example=False):
    sequences = parse_fasta(sequence_path(data_dir, task))
    if task == RESIDUE_TASK:
        samples, _, _, _ = load_task(data_dir, example=example)
        required = {sample["sequence_id"] for sample in samples.values()}
    else:
        samples, _, _ = load_pooled(data_dir, {task}, example=example)
        required = {
            sid
            for sample in samples.values()
            for sid in sample["sequence_ids"].split(";")
        }
    return {sid: sequences[sid] for sid in sorted(required)}


def main():
    parser = argparse.ArgumentParser(
        description="Prepare exact unique FASTA inputs for one scoring task"
    )
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--task", required=True, choices=sorted(TASKS | {RESIDUE_TASK}))
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--example", action="store_true")
    args = parser.parse_args()
    record = RunRecord(
        args.output_dir.resolve(), {"task": args.task, "example": args.example}
    )
    try:
        for p in task_input_paths(args.data_dir, args.task):
            record.track(p)
        sequences = required_sequences(args.data_dir, args.task, args.example)
        (record.output / "sequences.fasta").write_bytes(
            (
                "".join(f">{sid}\n{seq}\n" for sid, seq in sorted(sequences.items()))
            ).encode("utf-8")
        )
        record.data["unique_sequences"] = len(sequences)
    except Exception as error:
        record.finish(error)
        raise
    error = record.finish()
    if error:
        raise error
    print(
        json.dumps(
            {
                "status": "PASS",
                "unique_sequences": len(sequences),
                "fasta": str(record.output / "sequences.fasta"),
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
