from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
from .evaluate_predictions import eligible_frames, read_table, verify_dataset
from .model_interface import load_mean_representation
from .readers import parse_fasta, validate_sample_sequence_mappings


def required_sequences(data_dir: Path, tasks=None) -> tuple[dict[str, int], dict]:
    eligible, _ = eligible_frames(data_dir, tasks)
    sequences = parse_fasta(data_dir / 'sequences.fasta')
    sample_rows = {}
    sequence_info = {}
    for row in read_table(data_dir / 'sample_sequences.tsv.gz'):
        key = (row['task_id'], row['sample_id'])
        if key in eligible:
            sample_rows.setdefault(key, []).append(row)
    for key, sample in eligible.items():
        sequence_ids, lengths = validate_sample_sequence_mappings(key, sample, sample_rows.get(key, []))
        for sequence_id, length in zip(sequence_ids, lengths):
            sequence = sequences.get(sequence_id)
            if sequence is None or len(sequence) != length:
                raise ValueError(f'eligible input sequence is absent or has wrong length: {key}/{sequence_id}')
            if sequence_info.setdefault(sequence_id, length) != length:
                raise ValueError(f'inconsistent sequence length: {sequence_id}')
    return sequence_info, eligible


def scan_archives(archive_root: Path, sequence_info: dict[str, int], model_id: str, model_revision: str, layers: tuple[int, ...]) -> tuple[list[dict], int]:
    selected = {}
    dimensions = set()
    for path in sorted(archive_root.rglob('*.npz')):
        try:
            with np.load(path, allow_pickle=False) as archive:
                if 'meta' not in archive.files or not ({'mean_layers', 'sequence_layers'} & set(archive.files)):
                    continue
                meta = json.loads(str(archive['meta'].item()))
            sequence_id = meta.get('sequence_id')
            if sequence_id not in sequence_info:
                continue
            if meta.get('model') != model_id or meta.get('revision') != model_revision:
                continue
            length = sequence_info[sequence_id]
            _, dimension = load_mean_representation(path, sequence_id, length, layers, model_id, model_revision)
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
            continue
        if sequence_id in selected:
            raise ValueError(f'multiple archives identify the same sequence: {sequence_id}')
        selected[sequence_id] = path
        dimensions.add(dimension)
    missing = sorted(set(sequence_info) - set(selected))
    if missing:
        raise ValueError(f'representation archives missing {len(missing)} eligible sequences; first IDs: {missing[:10]}')
    if len(dimensions) != 1:
        raise ValueError(f'representation dimensions differ across archives: {sorted(dimensions)}')
    rows = [{'sequence_id': sequence_id, 'length': sequence_info[sequence_id],
             'file': path.relative_to(archive_root).as_posix()}
            for sequence_id, path in sorted(selected.items())]
    return rows, dimensions.pop()


def build(args) -> dict:
    data_dir = args.data_dir.resolve()
    archive_root = args.archive_root.resolve()
    output = args.output_manifest.resolve()
    if not data_dir.is_dir() or not archive_root.is_dir():
        raise FileNotFoundError('data directory or archive root does not exist')
    if output.exists():
        raise ValueError(f'refusing to overwrite existing manifest: {output}')
    if output.parent != archive_root:
        raise ValueError('output manifest must be directly inside --archive-root')
    for field, value in (('model_id', args.model_id), ('model_revision', args.model_revision), ('pooling', args.pooling)):
        if not value.strip():
            raise ValueError(f'{field} must be a nonempty string')
    verify_dataset(data_dir)
    sequence_info, eligible = required_sequences(data_dir, set(args.task) if args.task else None)
    archive_rows, dimension = scan_archives(archive_root, sequence_info, args.model_id, args.model_revision, tuple(args.layers))
    manifest = {'model_id': args.model_id, 'model_revision': args.model_revision,
                'pooling': args.pooling, 'candidate_layers': list(args.layers), 'archives': archive_rows}
    output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return {'status': 'PASS', 'manifest': str(output), 'eligible_samples': len(eligible),
            'unique_sequences': len(sequence_info), 'validated_archives': len(archive_rows),
            'representation_dimension': dimension}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--layers', type=int, nargs='+', required=True)
    parser.add_argument('--task', action='append', choices=['KNOT-PRESENCE', 'FS-L2-grouped-PU-ranking'])
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--archive-root', type=Path, required=True)
    parser.add_argument('--output-manifest', type=Path, required=True)
    parser.add_argument('--model-id', required=True)
    parser.add_argument('--model-revision', required=True)
    parser.add_argument('--pooling', required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(build(args), ensure_ascii=False))
    except Exception as error:
        print(f'Feature manifest rejected: {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
