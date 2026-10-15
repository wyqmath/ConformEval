from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
from collections import Counter
from pathlib import Path
EXPECTED_COUNTS = {'KNOT-PRESENCE': 734, 'FS-L2-grouped-PU-ranking': 8879}
EXPECTED_ARMS = {'PLM-nested', 'PLM-hidden17', 'PLM-hidden33', 'PLM-length', 'PLM-composition', 'PLM-nearest', 'PLM-training-prior'}
TEXT_PATH_MARKERS = re.compile('(?:/Users/|/tmp/|/private/var/|/home/|[A-Za-z]:\\\\)')

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()

def read_json(path: Path) -> dict:
    with path.open(encoding='utf-8') as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f'expected a JSON object: {path}')
    return value

def review_passed(path: Path) -> bool:
    text = path.read_text(encoding='utf-8')
    return bool(re.search('结论：\\*\\*PASS\\b|\\bPASS\\b[^\\n]{0,80}无阻塞项', text))

def verify_prediction_rows(path: Path) -> dict:
    expected_fields = {'task_id', 'outer_fold', 'arm', 'sample_id', 'binding_group_id', 'learning_label', 'learning_mask', 'biological_target', 'biological_mask', 'label', 'mask', 'score'}
    counts: Counter[str] = Counter()
    arms: set[str] = set()
    seen: set[tuple[str, str, str]] = set()
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if set(reader.fieldnames or ()) != expected_fields:
            raise ValueError('predictions.tsv schema differs from the audited runner')
        for row in reader:
            task = row['task_id']
            if task not in EXPECTED_COUNTS:
                raise ValueError(f'unexpected scored task: {task}')
            if row['arm'] not in EXPECTED_ARMS:
                raise ValueError(f"unexpected prediction arm: {row['arm']}")
            key = (task, row['sample_id'], row['arm'])
            if key in seen:
                raise ValueError(f'duplicate prediction key: {key}')
            seen.add(key)
            try:
                score = float(row['score'])
            except ValueError as exc:
                raise ValueError(f'non-numeric prediction score for {key}') from exc
            if not (score == score and abs(score) != float('inf')):
                raise ValueError(f'non-finite prediction score for {key}')
            arms.add(row['arm'])
            counts[task] += 1
    if arms != EXPECTED_ARMS:
        raise ValueError(f'prediction arms differ: {sorted(arms)}')
    expected_rows = {task: n * len(EXPECTED_ARMS) for task, n in EXPECTED_COUNTS.items()}
    if dict(counts) != expected_rows:
        raise ValueError(f'prediction row counts differ: {dict(counts)} != {expected_rows}')
    return {'rows_by_task': dict(counts), 'samples_by_task': EXPECTED_COUNTS, 'arms': sorted(arms)}

def gzip_copy(source: Path, destination: Path) -> None:
    with source.open('rb') as src, destination.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as dst:
            shutil.copyfileobj(src, dst)

def output_inventory(root: Path) -> dict:
    files = {}
    for path in sorted((p for p in root.rglob('*') if p.is_file() and p.name != 'manifest.json')):
        rel = path.relative_to(root).as_posix()
        files[rel] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
    return files

def scan_for_machine_paths(root: Path) -> None:
    for path in root.rglob('*'):
        if not path.is_file() or path.suffix.lower() == '.npz':
            continue
        opener = gzip.open if path.suffix.lower() == '.gz' else open
        try:
            with opener(path, 'rt', encoding='utf-8') as handle:
                if any((TEXT_PATH_MARKERS.search(line) for line in handle)):
                    raise ValueError(f'machine-specific absolute path found in exported text: {path.name}')
        except UnicodeDecodeError:
            continue

def export(source_results: Path, protocol: Path, data_manifest: Path, evaluator_script: Path, review_log: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise ValueError(f'refusing to overwrite existing output directory: {output_dir}')
    for path in (source_results, protocol, data_manifest, evaluator_script, review_log):
        if not path.exists():
            raise FileNotFoundError(path)
    run = read_json(source_results / 'run.json')
    if run.get('status') != 'COMPLETE_PENDING_INDEPENDENT_REVIEW':
        raise ValueError(f"unexpected source run status: {run.get('status')}")
    if not review_passed(review_log):
        raise ValueError('independent review log does not contain a PASS conclusion')
    source_hashes = run.get('outputs_sha256', {})
    selected = ['case_ranks.tsv', 'configuration.json', 'fold_metrics.tsv', 'inner_selection.tsv', 'metrics_FS-L2-grouped-PU-ranking.json', 'metrics_KNOT-PRESENCE.json', 'nearest_references_FS-L2-grouped-PU-ranking.json', 'nearest_references_KNOT-PRESENCE.json', 'outer_selections.tsv', 'predictions.tsv', 'warnings.tsv']
    selected.extend(sorted((k for k in source_hashes if k.startswith('reader_parameters/'))))
    selected.append('fold_roles.tsv')
    for rel in selected:
        source = source_results / rel
        expected = source_hashes.get(rel)
        if not expected or not source.is_file() or sha256(source) != expected:
            raise ValueError(f'missing or changed reviewed output: {rel}')
    row_audit = verify_prediction_rows(source_results / 'predictions.tsv')
    metrics = {task: read_json(source_results / f'metrics_{task}.json') for task in EXPECTED_COUNTS}
    for task, record in metrics.items():
        if record.get('task_id') != task or record.get('sample_count') != EXPECTED_COUNTS[task]:
            raise ValueError(f'metric file disagrees with the fixed task denominator: {task}')
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix=f'.{output_dir.name}.', dir=output_dir.parent))
    try:
        for rel in selected:
            source = source_results / rel
            if rel == 'predictions.tsv':
                gzip_copy(source, temp_dir / 'predictions.tsv.gz')
            elif rel == 'fold_roles.tsv':
                gzip_copy(source, temp_dir / 'fold_roles.tsv.gz')
            else:
                destination = temp_dir / rel
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
        (temp_dir / 'README.md').write_text('# Existing ESM-2 reference results\n\nAudited retrospective grouped-development results for KNOT-PRESENCE and FS-L2-grouped-PU-ranking. This directory contains predictions, metrics, fold roles, selections, nearest-reference records, and fitted reader parameters. It does not contain model weights or embedding caches. The result scope and limitations are documented in the benchmark-level data card.\n', encoding='utf-8')
        provenance = {'artifact_id': 'esm2-650m-existing-model-reference-v0.3', 'benchmark_version': '0.3-existing-model-grouped-development', 'evaluation_scope': 'retrospective_development_grouped_evaluation', 'independent_review': 'PASS', 'model_testing_new_models': False, 'source_run_status': run['status'], 'started_at': run['started_at'], 'ended_at': run['ended_at'], 'encoder': {'id': 'facebook/esm2_t33_650M_UR50D', 'revision': '08e4846e537177426273712802403f7ba8261b6c', 'layers': [5, 11, 17, 23, 29, 33]}, 'data_lock_sha256': sha256(data_manifest), 'protocol_sha256': sha256(protocol), 'evaluator_sha256': sha256(evaluator_script), 'review_record_sha256': sha256(review_log), 'source_input_hashes': {'data_v06_manifest': run['input_sha256']['results/benchmark_fixed/20261014/data_v06/manifest.json'], 'data_v06_protocol': run['input_sha256']['results/benchmark_fixed/20261014/data_v06/protocol.json'], 'accepted_model_protocol': run['input_sha256']['configs/physinfobench_existing_20261012.json'], 'representation_overlay': run['input_sha256']['results/benchmark_fixed/20261014/cache_completion_v01/knot_presence/feature_coverage_overlay.tsv.gz']}, 'versions': run.get('dependencies', run.get('versions', {})), 'prediction_coverage': row_audit, 'metric_files': {task: {'sample_count': record['sample_count'], 'outer_folds': record['outer_folds']} for task, record in metrics.items()}, 'limitations': ['Only KNOT-PRESENCE and FS-L2 have scorable frozen folds in this run.', 'FS-L2 unlabeled background is not a verified biological negative class.', 'Intervals are conditional on fixed outer-fold predictions and do not include refitting uncertainty.', 'This retrospective grouped development evaluation does not certify family-independent or prospective generalization.']}
        (temp_dir / 'provenance.json').write_text(json.dumps(provenance, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        scan_for_machine_paths(temp_dir)
        manifest = {'schema_version': '1.0', 'artifact_id': provenance['artifact_id'], 'files': output_inventory(temp_dir), 'model_weights_included': False, 'embedding_cache_included': False}
        (temp_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        os.rename(temp_dir, output_dir)
        return {'output_dir': str(output_dir), 'files': len(manifest['files']), 'prediction_rows': row_audit['rows_by_task']}
    except Exception:
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise

def main() -> int:
    parser = argparse.ArgumentParser(description='Export the audited ESM-2 reference run into a portable benchmark artifact.')
    parser.add_argument('--source-results', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--data-manifest', type=Path, required=True)
    parser.add_argument('--evaluator-script', type=Path, required=True)
    parser.add_argument('--review-log', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    result = export(args.source_results, args.protocol, args.data_manifest, args.evaluator_script, args.review_log, args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
