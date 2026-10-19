from __future__ import annotations
import argparse
import json
from pathlib import Path
from .evaluate_predictions import TASKS, eligible_frames, read_table, verify_dataset
from .evaluate_residue_predictions import TASK as RESIDUE_TASK, load_task
from .readers import parse_fasta, validate_sample_sequence_mappings
from .runtime import RunRecord


def required_sequences(data_dir: Path, task: str, example=False):
    if task == RESIDUE_TASK:
        samples, labels, folds, _ = load_task(data_dir, example=example)
        all_sequences = parse_fasta(data_dir / 'inputs/sequences.fasta')
        required = {s['sequence_id']: all_sequences[s['sequence_id']] for s in samples.values()}
    else:
        verify_dataset(data_dir)
        is_example=json.loads((data_dir/'dataset_manifest.json').read_text(encoding='utf-8-sig')).get('scope')=='software_example'
        if is_example != example:
            raise ValueError('software example requires explicit --example; benchmark data must omit it')
        eligible, _ = eligible_frames(data_dir, {task})
        all_sequences = parse_fasta(data_dir / 'sequences.fasta')
        mapping = {}
        for row in read_table(data_dir / 'sample_sequences.tsv.gz'):
            key = (row['task_id'], row['sample_id'])
            if key in eligible:
                mapping.setdefault(key, []).append(row)
        required = {}
        for key, sample in eligible.items():
            sequence_ids, lengths = validate_sample_sequence_mappings(key, sample, mapping.get(key, []))
            for sid, length in zip(sequence_ids, lengths):
                sequence = all_sequences.get(sid)
                if sequence is None or len(sequence) != length:
                    raise ValueError(f'sequence identity/length mismatch: {sid}')
                required[sid] = sequence
    return required


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
