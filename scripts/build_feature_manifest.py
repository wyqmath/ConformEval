from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
try:
    from .evaluate_predictions import eligible_frames, read_table, verify_dataset
    from .physinfobench_model_interface import load_mean_representation
    from .run_nested_readers import LAYERS, parse_fasta, validate_sample_sequence_mappings
except ImportError:
    from evaluate_predictions import eligible_frames, read_table, verify_dataset
    from physinfobench_model_interface import load_mean_representation
    from run_nested_readers import LAYERS, parse_fasta, validate_sample_sequence_mappings

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def required_sequences(data_dir: Path) -> tuple[dict[str, tuple[str, int]], dict]:
    eligible, _ = eligible_frames(data_dir)
    sequences = parse_fasta(data_dir / 'sequences.fasta')
    sample_rows = {}
    sequence_info: dict[str, tuple[str, int]] = {}
    for row in read_table(data_dir / 'sample_sequences.tsv.gz'):
        key = (row['task_id'], row['sample_id'])
        if key not in eligible:
            continue
        sample_rows.setdefault(key, []).append(row)
    for key, sample in eligible.items():
        hashes, lengths = validate_sample_sequence_mappings(key, sample, sample_rows.get(key, []))
        for digest, length in zip(hashes, lengths):
            sequence_id = f'SEQ_{digest}'
            sequence = sequences.get(sequence_id)
            if sequence is None or len(sequence) != length:
                raise ValueError(f'eligible input sequence is absent or has wrong length: {key}/{sequence_id}')
            actual = hashlib.sha256(sequence.encode('ascii')).hexdigest()
            if actual != digest:
                raise ValueError(f'eligible input sequence hash differs from FASTA: {key}/{sequence_id}')
            old = sequence_info.setdefault(digest, (sequence_id, length))
            if old != (sequence_id, length):
                raise ValueError(f'sequence hash maps to inconsistent identity or length: {digest}')
    return (sequence_info, eligible)

def scan_archives(archive_root: Path, sequence_info: dict[str, tuple[str, int]], model_id: str, model_revision: str) -> tuple[list[dict], int]:
    selected: dict[str, tuple[Path, str]] = {}
    dimensions = set()
    required = set(sequence_info)
    for path in sorted(archive_root.rglob('*.npz')):
        try:
            with np.load(path, allow_pickle=False) as archive:
                if not {'meta', 'mean_layers'}.issubset(archive.files):
                    continue
                meta = json.loads(str(archive['meta'].item()))
            digest = meta.get('seq_sha256')
            if digest not in required:
                continue
            if meta.get('model') != model_id or meta.get('revision') != model_revision:
                continue
            _, length = sequence_info[digest]
            if meta.get('len_full') != length or meta.get('len_used') != length or meta.get('truncated') not in (False, 0, 'false'):
                continue
            expected_hash = sha256(path)
            _, dimension = load_mean_representation(path, digest, length, expected_hash, LAYERS, model_id, model_revision)
            dimensions.add(dimension)
            old = selected.get(digest)
            if old is not None and old[1] != expected_hash:
                raise ValueError(f'multiple different archives map to required sequence {digest}; keep one reviewed archive or make the duplicate byte-identical')
            selected.setdefault(digest, (path, expected_hash))
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
            if isinstance(error, ValueError) and str(error).startswith('multiple different'):
                raise
            continue
    missing = sorted(required - set(selected))
    if missing:
        raise ValueError(f'representation archives missing {len(missing)} eligible sequences; first IDs: {missing[:10]}')
    if len(dimensions) != 1:
        raise ValueError(f'representation dimensions differ across archives: {sorted(dimensions)}')
    rows = [{'sequence_sha256': digest, 'length': sequence_info[digest][1], 'file': path.relative_to(archive_root).as_posix(), 'archive_sha256': archive_hash} for digest, (path, archive_hash) in sorted(selected.items())]
    return (rows, dimensions.pop())

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
    if len(args.extractor_code_sha256) != 64:
        raise ValueError('extractor-code-sha256 must be a 64-character SHA-256')
    try:
        int(args.extractor_code_sha256, 16)
    except ValueError as exc:
        raise ValueError('extractor-code-sha256 must be hexadecimal') from exc
    verify_dataset(data_dir)
    sequence_info, eligible = required_sequences(data_dir)
    archive_rows, dimension = scan_archives(archive_root, sequence_info, args.model_id, args.model_revision)
    manifest = {'model_id': args.model_id, 'model_revision': args.model_revision, 'pooling': args.pooling, 'extractor_code_sha256': args.extractor_code_sha256.lower(), 'candidate_layers': list(LAYERS), 'archives': archive_rows}
    output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return {'status': 'PASS', 'manifest': str(output), 'eligible_samples': len(eligible), 'unique_sequences': len(sequence_info), 'validated_archives': len(archive_rows), 'representation_dimension': dimension}

def main() -> int:
    parser = argparse.ArgumentParser(description='Build a checked manifest for exact-sequence mean-representation archives.')
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--archive-root', type=Path, required=True, help='directory containing one exact-sequence .npz archive per required sequence')
    parser.add_argument('--output-manifest', type=Path, required=True, help='new JSON manifest path directly inside --archive-root')
    parser.add_argument('--model-id', required=True)
    parser.add_argument('--model-revision', required=True)
    parser.add_argument('--pooling', required=True)
    parser.add_argument('--extractor-code-sha256', required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(build(args), ensure_ascii=False))
    except Exception as error:
        print(f'Feature manifest rejected: {error}', file=sys.stderr)
        return 2
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
