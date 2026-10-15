from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path
PATH_MARKERS = re.compile('(?:/Users/|/tmp/|/private/var/|/home/|[A-Za-z]:\\\\)')
EXPECTED = {'KNOT-PRESENCE': 734 * 7, 'FS-L2-grouped-PU-ranking': 8879 * 7}
REFERENCE_ARMS = {'PLM-composition', 'PLM-hidden17', 'PLM-hidden33', 'PLM-length', 'PLM-nearest', 'PLM-nested', 'PLM-training-prior'}
PREDICTION_IDENTITY_FIELDS = ('binding_group_id', 'learning_label', 'learning_mask', 'biological_target', 'biological_mask', 'label', 'mask')

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()

def verify_manifest(root: Path, manifest_path: Path) -> int:
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    records = manifest.get('files', {})
    if not records:
        raise ValueError(f'empty file manifest: {manifest_path}')
    for name, record in records.items():
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError(f'unsafe path in manifest: {name}')
        path = root / rel
        expected_hash = record['sha256'] if isinstance(record, dict) else record
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f'release file hash mismatch: {name}')
        if isinstance(record, dict) and path.stat().st_size != record.get('bytes'):
            raise ValueError(f'release file size mismatch: {name}')
    return len(records)

def scan_text_for_paths(root: Path) -> int:
    scanned = 0
    for path in root.rglob('*'):
        if any((part in {'.venv', '.git', 'build', 'dist', 'runs', 'wheelhouse', '__pycache__'} or part.endswith('.egg-info') or part.startswith('.venv_') for part in path.relative_to(root).parts)):
            continue
        if not path.is_file() or path.suffix.lower() in {'.py', '.pyc', '.npz', '.pdf', '.png', '.jpg'}:
            continue
        opener = gzip.open if path.suffix.lower() == '.gz' else open
        try:
            with opener(path, 'rt', encoding='utf-8') as handle:
                for line in handle:
                    if PATH_MARKERS.search(line):
                        raise ValueError(f'machine-specific path found in {path.relative_to(root)}')
            scanned += 1
        except UnicodeDecodeError:
            continue
    return scanned

def verify_predictions(path: Path) -> dict:
    counts = {task: 0 for task in EXPECTED}
    seen = set()
    with gzip.open(path, 'rt', newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'task_id', 'outer_fold', 'arm', 'sample_id', 'learning_label', 'biological_target', 'biological_mask', 'label', 'mask', 'score'}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError('prediction table is missing required columns')
        for row in reader:
            task = row['task_id']
            if task not in EXPECTED:
                raise ValueError(f'unexpected scored task in predictions: {task}')
            key = (task, row['outer_fold'], row['arm'], row['sample_id'])
            if key in seen:
                raise ValueError(f'duplicate prediction identity: {key}')
            seen.add(key)
            try:
                score = float(row['score'])
            except ValueError as exc:
                raise ValueError(f'invalid prediction score: {key}') from exc
            if not math.isfinite(score):
                raise ValueError(f'non-finite prediction score: {key}')
            if task == 'FS-L2-grouped-PU-ranking' and row['learning_label'] == '0':
                if row['biological_target'] != 'null' or row['biological_mask'] != '0' or row['mask'] != '0' or (row['label'] != 'NA'):
                    raise ValueError(f'FS-L2 unlabeled background was converted to a biological negative: {key}')
            counts[task] += 1
    if counts != EXPECTED:
        raise ValueError(f'prediction denominators differ: {counts} != {EXPECTED}')
    return counts

def verify_prediction_assignments(predictions_path: Path, data_dir: Path, expected_tasks: set[str] | None=None, expected_arms: set[str] | None=None) -> dict:
    tasks = set(EXPECTED) if expected_tasks is None else set(expected_tasks)
    arms = set(REFERENCE_ARMS) if expected_arms is None else set(expected_arms)
    if not tasks or not arms:
        raise ValueError('prediction audit requires nonempty task and arm sets')
    sample_path = data_dir / 'task_samples.tsv.gz'
    with gzip.open(sample_path, 'rt', newline='', encoding='utf-8') as handle:
        sample_rows = list(csv.DictReader(handle, delimiter='\t'))
    eligible = {}
    for row in sample_rows:
        if row['task_id'] not in tasks or row['data_eligible'] != '1':
            continue
        key = (row['task_id'], row['sample_id'])
        if key in eligible:
            raise ValueError(f'duplicate eligible sample identity in data: {key}')
        eligible[key] = row
    fold_path = data_dir / 'folds/outer_group_holdout.tsv'
    group_fold = {}
    with fold_path.open(newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['task_id'] not in tasks:
                continue
            if row['role'] != 'test':
                raise ValueError(f'unexpected outer-fold role in released test map: {row}')
            key = (row['task_id'], row['group_id'])
            if key in group_fold:
                raise ValueError(f'binding group has multiple outer test folds: {key}')
            group_fold[key] = row['outer_fold']
    observed: dict[tuple[str, str], set[str]] = {}
    actual_samples = set()
    with gzip.open(predictions_path, 'rt', newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            task = row['task_id']
            if task not in tasks:
                raise ValueError(f'unexpected task in prediction assignment audit: {task}')
            sample_key = (task, row['sample_id'])
            sample = eligible.get(sample_key)
            if sample is None:
                raise ValueError(f'prediction sample is absent or ineligible in data: {sample_key}')
            expected_fold = group_fold.get((task, row['binding_group_id']))
            if expected_fold is None or row['outer_fold'] != expected_fold:
                raise ValueError(f'prediction is not assigned to its frozen outer test fold: {sample_key}')
            for field in PREDICTION_IDENTITY_FIELDS:
                if row[field] != sample[field]:
                    raise ValueError(f'prediction/data {field} mismatch: {sample_key}')
            arm = row['arm']
            if arm not in arms:
                raise ValueError(f'unexpected prediction arm: {arm}')
            arm_set = observed.setdefault(sample_key, set())
            if arm in arm_set:
                raise ValueError(f'duplicate prediction for sample and arm: {sample_key}/{arm}')
            arm_set.add(arm)
            actual_samples.add(sample_key)
    expected_samples = set(eligible)
    if actual_samples != expected_samples:
        missing = sorted(expected_samples - actual_samples)[:10]
        extra = sorted(actual_samples - expected_samples)[:10]
        raise ValueError(f'prediction sample coverage differs from eligible data; missing={missing}; extra={extra}')
    for sample_key, sample_arms in observed.items():
        if sample_arms != arms:
            raise ValueError(f'prediction arm coverage differs for {sample_key}: {sorted(sample_arms)}')
    return {task: {'eligible_samples': sum((1 for sample_task, _ in expected_samples if sample_task == task)), 'arms_per_sample': len(arms)} for task in sorted(tasks)}

def verify_task_catalog(root: Path) -> dict:
    catalog_path = root / 'task_catalog.tsv'
    with catalog_path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    count_fields = ('candidate_count_unit', 'eligible_count_unit', 'scored_count_unit')
    if not rows or any((field not in rows[0] for field in count_fields)):
        raise ValueError('task catalog must declare a unit for every denominator')
    for row in rows:
        if any((not row.get(field, '').strip() for field in count_fields)):
            raise ValueError(f"task catalog has an unspecified denominator unit: {row['task_id']}")
    task_ids = [row['task_id'] for row in rows]
    if len(rows) != 13 or len(task_ids) != len(set(task_ids)):
        raise ValueError('task catalog must contain 13 unique registered tasks')
    scored = {row['task_id']: row for row in rows if int(row['scored_units']) > 0}
    if set(scored) != set(EXPECTED) | {'DISORDER-STRUCTURE-REFERENCE-CANDIDATE'}:
        raise ValueError(f'scored task set differs from included reference results: {sorted(scored)}')
    for task, row in scored.items():
        if row['status'] != 'evaluated_existing_ESM2':
            raise ValueError(f"unexpected scored task status: {task}/{row['status']}")
    for row in rows:
        rel = row['result_artifact'].split('#', 1)[0]
        path = Path(rel)
        if path.is_absolute() or '..' in path.parts or (not (root / path).is_file()):
            raise ValueError(f"task catalog artifact path is invalid: {row['task_id']}/{rel}")
    by_task = {row['task_id']: row for row in rows}
    expected_case_counts = {'FS-L1-PAIRED': ('16', 'candidate_case', '10', 'accepted_case'), 'FS-REGION': ('10', 'case', '3', 'case_with_supported_region')}
    for task, expected in expected_case_counts.items():
        row = by_task[task]
        actual = (row['candidate_units'], row['candidate_count_unit'], row['eligible_units'], row['eligible_count_unit'])
        if actual != expected:
            raise ValueError(f'case denominator or count unit changed: {task}: {actual} != {expected}')
    return {'registered_tasks': len(rows), 'scored_tasks': sorted(scored)}

def verify_knot_representation_coverage(data_dir: Path) -> dict:
    coverage_path = data_dir / 'representation_coverage.tsv.gz'
    with gzip.open(coverage_path, 'rt', newline='', encoding='utf-8') as handle:
        coverage_rows = list(csv.DictReader(handle, delimiter='\t'))
    coverage = {(row['task_id'], row['sample_id'], row['sequence_sha256']): row for row in coverage_rows}
    ready_counts = {}
    sample_path = data_dir / 'task_samples.tsv.gz'
    with gzip.open(sample_path, 'rt', newline='', encoding='utf-8') as handle:
        samples = list(csv.DictReader(handle, delimiter='\t'))
    for task in ('KNOT-PRESENCE', 'KNOT-TYPE'):
        eligible = [row for row in samples if row['task_id'] == task and row['data_eligible'] == '1']
        if not eligible:
            raise ValueError(f'no eligible rows in released task data: {task}')
        for row in eligible:
            hashes = row['sequence_sha256s'].split(';')
            if len(hashes) != 1:
                raise ValueError(f"knot sample does not map to exactly one sequence: {task}/{row['sample_id']}")
            found = coverage.get((task, row['sample_id'], hashes[0]))
            if found is None or found['feature_ready'] != '1':
                raise ValueError(f"eligible knot sample lacks representation coverage: {task}/{row['sample_id']}")
        ready_counts[task] = len(eligible)
    if ready_counts['KNOT-PRESENCE'] != 734 or ready_counts['KNOT-TYPE'] != 65:
        raise ValueError(f'knot eligible denominators changed: {ready_counts}')
    return ready_counts

def audit_group_splits(eligible_groups: dict[str, set[str]], outer_rows: list[dict[str, str]], inner_rows: list[dict[str, str]], fold_counts: dict[str, int], inner_fold_count: int=3) -> dict:
    outer_by_task_fold: dict[tuple[str, str], set[str]] = {}
    outer_seen: dict[str, set[str]] = {task: set() for task in eligible_groups}
    for row in outer_rows:
        task, fold, group = (row['task_id'], row['outer_fold'], row['group_id'])
        if task not in eligible_groups or row['role'] != 'test':
            raise ValueError(f'unexpected outer split row: {task}/{fold}/{group}')
        bucket = outer_by_task_fold.setdefault((task, fold), set())
        if group in bucket or group in outer_seen[task]:
            raise ValueError(f'group assigned to more than one outer test fold: {task}/{group}')
        bucket.add(group)
        outer_seen[task].add(group)
    for task, expected in eligible_groups.items():
        if outer_seen[task] != expected:
            raise ValueError(f'outer split groups differ from eligible groups: {task}')
        folds = sorted({fold for row_task, fold in outer_by_task_fold if row_task == task}, key=int)
        if len(folds) != fold_counts[task] or folds != [str(i) for i in range(fold_counts[task])]:
            raise ValueError(f'outer fold IDs/count differ from catalog: {task}')
    inner_validation: dict[tuple[str, str, str], set[str]] = {}
    for row in inner_rows:
        task, outer_fold, inner_fold, group = (row['task_id'], row['outer_fold'], row['inner_fold'], row['group_id'])
        if task not in eligible_groups or row['role'] != 'validation':
            raise ValueError(f'unexpected inner split row: {task}/{outer_fold}/{inner_fold}/{group}')
        key = (task, outer_fold, inner_fold)
        bucket = inner_validation.setdefault(key, set())
        if group in bucket:
            raise ValueError(f'duplicate inner validation group: {key}/{group}')
        bucket.add(group)
    audit = {}
    for task, groups in eligible_groups.items():
        for outer_fold in (str(i) for i in range(fold_counts[task])):
            outer_test = outer_by_task_fold[task, outer_fold]
            outer_train = groups - outer_test
            validation_sets = []
            for inner_fold in (str(i) for i in range(inner_fold_count)):
                key = (task, outer_fold, inner_fold)
                validation = inner_validation.get(key)
                if not validation:
                    raise ValueError(f'empty or missing inner validation fold: {key}')
                if validation & outer_test:
                    raise ValueError(f'outer test group leaked into inner validation: {key}')
                if not validation.issubset(outer_train):
                    raise ValueError(f'unknown group in inner validation: {key}')
                validation_sets.append(validation)
            for i, left in enumerate(validation_sets):
                if any((left & right for right in validation_sets[i + 1:])):
                    raise ValueError(f'inner validation groups overlap: {task}/{outer_fold}')
            if set.union(*validation_sets) != outer_train:
                raise ValueError(f'inner validation folds do not partition outer training groups: {task}/{outer_fold}')
        audit[task] = {'eligible_groups': len(groups), 'outer_folds': fold_counts[task], 'inner_folds_per_outer': inner_fold_count}
    return audit

def verify_group_splits(data_dir: Path, scored_tasks: dict[str, dict[str, str]]) -> dict:
    with gzip.open(data_dir / 'task_samples.tsv.gz', 'rt', newline='', encoding='utf-8') as handle:
        samples = list(csv.DictReader(handle, delimiter='\t'))
    eligible_groups = {task: {row['binding_group_id'] for row in samples if row['task_id'] == task and row['data_eligible'] == '1'} for task in scored_tasks}
    outer_path = data_dir / 'folds/outer_group_holdout.tsv'
    with outer_path.open(newline='', encoding='utf-8') as handle:
        outer_rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['task_id'] in scored_tasks]
    inner_path = data_dir / 'folds/inner_validation_groups.tsv.gz'
    with gzip.open(inner_path, 'rt', newline='', encoding='utf-8') as handle:
        inner_rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['task_id'] in scored_tasks]
    fold_counts = {task: int(row['outer_folds']) for task, row in scored_tasks.items()}
    return audit_group_splits(eligible_groups, outer_rows, inner_rows, fold_counts)

def verify_case_source_references(data_dir: Path) -> dict:
    with (data_dir / 'case_modules.tsv').open(newline='', encoding='utf-8') as handle:
        modules = list(csv.DictReader(handle, delimiter='\t'))
    with (data_dir / 'case_sources.tsv').open(newline='', encoding='utf-8') as handle:
        sources = list(csv.DictReader(handle, delimiter='\t'))
    if not modules or not sources:
        raise ValueError('case module or source registry is empty')
    if 'source_path' in (modules[0].keys() if modules else ()):
        raise ValueError('project-local source_path field remains in case modules')
    source_by_ref = {}
    for row in sources:
        ref = row.get('source_ref', '')
        if not ref or ref in source_by_ref:
            raise ValueError(f'missing or duplicate case source reference: {ref}')
        if not re.fullmatch('CS\\d{4}', ref):
            raise ValueError(f'invalid portable case source reference: {ref}')
        if not re.fullmatch('[0-9a-f]{64}', row.get('source_artifact_sha256', '')):
            raise ValueError(f'invalid source artifact fingerprint: {ref}')
        source_by_ref[ref] = row
    seen_refs = set()
    for row in modules:
        ref = row.get('source_ref', '')
        source = source_by_ref.get(ref)
        if source is None:
            raise ValueError(f'case module source reference does not resolve: {ref}')
        if (source['module'], source['case_id']) != (row['module'], row['case_id']):
            raise ValueError(f'case source reference identifies a different record: {ref}')
        if ref in seen_refs:
            raise ValueError(f'case source reference is reused: {ref}')
        seen_refs.add(ref)
    if seen_refs != set(source_by_ref):
        raise ValueError('case source registry contains unreferenced entries')
    return {'case_rows': len(modules), 'resolvable_source_records': len(sources)}

def verify(root: Path) -> dict:
    data = root / 'data/v0.3.1'
    results = root / 'results/esm2_existing_v0.3'
    data_manifest = json.loads((data / 'dataset_manifest.json').read_text(encoding='utf-8'))
    if data_manifest.get('distribution_status') != 'DRAFT_NOT_CLEARED_FOR_PUBLIC_RELEASE':
        raise ValueError('unexpected data redistribution status; requires maintainer review')
    data_files = verify_manifest(data, data / 'dataset_manifest.json')
    result_manifest = json.loads((results / 'manifest.json').read_text(encoding='utf-8'))
    result_files = verify_manifest(results, results / 'manifest.json')
    prediction_counts = verify_predictions(results / 'predictions.tsv.gz')
    prediction_assignments = verify_prediction_assignments(results / 'predictions.tsv.gz', data)
    catalog = verify_task_catalog(root)
    knot_coverage = verify_knot_representation_coverage(data)
    scored_rows = {}
    with (root / 'task_catalog.tsv').open(newline='', encoding='utf-8') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['task_id'] in EXPECTED:
                scored_rows[row['task_id']] = row
    split_audit = verify_group_splits(data, scored_rows)
    case_sources = verify_case_source_references(data)
    for task, expected_n in (('KNOT-PRESENCE', 734), ('FS-L2-grouped-PU-ranking', 8879)):
        metrics = json.loads((results / f'metrics_{task}.json').read_text(encoding='utf-8'))
        if metrics.get('sample_count') != expected_n:
            raise ValueError(f'metric denominator mismatch: {task}')
    path_scan = scan_text_for_paths(root)
    return {'status': 'PASS_LOCAL_DRAFT', 'distribution_status': data_manifest['distribution_status'], 'dataset_files_verified': data_files, 'result_files_verified': result_files, 'prediction_rows': prediction_counts, 'prediction_sample_arm_fold_audit': prediction_assignments, 'task_catalog': catalog, 'group_split_audit': split_audit, 'case_source_references': case_sources, 'eligible_knot_representations': knot_coverage, 'text_files_scanned_for_local_paths': path_scan, 'result_model_weights_included': result_manifest.get('model_weights_included'), 'result_embedding_cache_included': result_manifest.get('embedding_cache_included')}

def main() -> int:
    parser = argparse.ArgumentParser(description='Verify file integrity and core label semantics in a PhysInfoBench release tree.')
    parser.add_argument('--benchmark-root', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.benchmark_root.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
