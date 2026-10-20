from __future__ import annotations
import argparse
import json
from pathlib import Path
from .evaluate_predictions import TASKS
from .evaluate_residue_predictions import TASK as RESIDUE_TASK, load_task
from .datasets import load_pooled, parse_fasta, sequence_path
from .runtime import RunRecord


def required_sequences(data_dir: Path, task: str, example=False):
    sequences=parse_fasta(sequence_path(data_dir,task))
    if task==RESIDUE_TASK:
        samples,_,_,_=load_task(data_dir,example=example)
        required={sample['sequence_id'] for sample in samples.values()}
    else:
        samples,_=load_pooled(data_dir,{task})
        required={sid for sample in samples.values() for sid in sample['sequence_ids'].split(';')}
    return {sid:sequences[sid] for sid in sorted(required)}


def main():
    parser = argparse.ArgumentParser(description='Prepare exact unique FASTA inputs for one scoring task')
    parser.add_argument('--data-dir', required=True, type=Path)
    parser.add_argument('--task', required=True, choices=sorted(TASKS | {RESIDUE_TASK}))
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--example', action='store_true')
    args = parser.parse_args()
    record = RunRecord(args.output_dir.resolve(), {'task':args.task,'example':args.example})
    try:
        for p in args.data_dir.rglob('*'):
            if p.is_file() and p.suffix != '.md':
                record.track(p)
        record.track(sequence_path(args.data_dir,args.task))
        sequences = required_sequences(args.data_dir, args.task, args.example)
        (record.output / 'sequences.fasta').write_bytes((''.join(f'>{sid}\n{seq}\n' for sid,seq in sorted(sequences.items()))).encode('utf-8'))
        record.data['unique_sequences'] = len(sequences)
    except Exception as error:
        record.finish(error)
        raise
    error = record.finish()
    if error:
        raise error
    print(json.dumps({'status':'PASS','unique_sequences':len(sequences),'fasta':str(record.output / 'sequences.fasta')}))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
