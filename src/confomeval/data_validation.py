"""Scientific integrity checks for versioned data and reference predictions.

These functions validate scoring inputs. They do not register repository files,
restrict source comments, or approve a delivery.
"""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path
from .evaluate_residue_predictions import check_manifest

EXPECTED = {'KNOT-PRESENCE': 734 * 7, 'FS-L2-grouped-PU-ranking': 8879 * 7}
REFERENCE_ARMS = {'PLM-composition', 'PLM-hidden17', 'PLM-hidden33', 'PLM-length', 'PLM-nearest', 'PLM-nested', 'PLM-training-prior'}
PREDICTION_IDENTITY_FIELDS = ('binding_group_id', 'learning_label', 'learning_mask', 'biological_target', 'biological_mask', 'label', 'mask')

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()

def read_tsv(path: Path, compressed: bool=False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def write_gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as zipped:
            with __import__('io').TextIOWrapper(zipped, encoding='utf-8', newline='') as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)

def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = None
    sequence: list[str] = []
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if line.startswith('>'):
            if identifier is not None:
                if identifier in records:
                    raise ValueError(f'duplicate FASTA ID: {identifier}')
                records[identifier] = ''.join(sequence)
            identifier, sequence = (line[1:].split()[0], [])
        elif line.strip():
            sequence.append(line.strip())
    if identifier is not None:
        if identifier in records:
            raise ValueError(f'duplicate FASTA ID: {identifier}')
        records[identifier] = ''.join(sequence)
    if not records or any((not seq for seq in records.values())):
        raise ValueError('FASTA is empty or contains an empty sequence')
    return records

def verify_sequences(data_dir: Path) -> dict:
    fasta = fasta_records(data_dir / 'sequences.fasta')
    rows = read_tsv(data_dir / 'sample_sequences.tsv.gz', compressed=True)
    for row in rows:
        sequence = fasta.get(row['sequence_id'])
        if sequence is None:
            raise ValueError(f"sample sequence absent from FASTA: {row['sequence_id']}")
        if len(sequence) != int(row['length']):
            raise ValueError(f"sequence length mismatch: {row['sequence_id']}")
        digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        if digest != row['sequence_sha256']:
            raise ValueError(f"sequence hash mismatch: {row['sequence_id']}")
    return {'fasta_sequences': len(fasta), 'sample_sequence_rows': len(rows)}

def verify_masks(data_dir: Path) -> dict:
    counts = {}
    for filename in ('position_labels_strict.tsv.gz', 'position_labels_structure_reference_candidate.tsv.gz'):
        rows = read_tsv(data_dir / filename, compressed=True)
        for row in rows:
            label, mask = (row['label'], row['mask'])
            if mask not in {'0', '1'}:
                raise ValueError(f'invalid mask in {filename}: {mask}')
            if mask == '0' and label != 'NA':
                raise ValueError(f'unknown residue label is not NA/masked in {filename}')
            if mask == '1' and label not in {'0', '1'}:
                raise ValueError(f'active residue label is invalid in {filename}')
        counts[filename] = {'rows': len(rows), 'masked_unknown': sum((r['mask'] == '0' for r in rows))}
    return counts

def verify_coordinates(data_dir: Path) -> dict:
    fasta = fasta_records(data_dir / 'sequences.fasta')
    lengths_by_sha = {identifier.removeprefix('SEQ_'): len(sequence) for identifier, sequence in fasta.items()}
    audited = {}
    for filename in ('position_labels_strict.tsv.gz', 'position_labels_structure_reference_candidate.tsv.gz'):
        rows = read_tsv(data_dir / filename, compressed=True)
        seen = set()
        for row in rows:
            digest = row['sequence_sha256']
            if digest not in lengths_by_sha:
                raise ValueError(f'residue label sequence is absent from FASTA in {filename}: {digest}')
            try:
                position = int(row['position'])
            except ValueError as exc:
                raise ValueError(f"non-integer residue coordinate in {filename}: {row['position']}") from exc
            if not 1 <= position <= lengths_by_sha[digest]:
                raise ValueError(f'residue coordinate outside sequence in {filename}: {digest}:{position}')
            key = (row['sample_id'], digest, position)
            if key in seen:
                raise ValueError(f'duplicate residue label coordinate in {filename}: {key}')
            seen.add(key)
        audited[filename] = {'rows': len(rows), 'unique_coordinates': len(seen)}
    sequence_rows = read_tsv(data_dir / 'sample_sequences.tsv.gz', compressed=True)
    lengths_by_sample: dict[str, set[int]] = {}
    for row in sequence_rows:
        lengths_by_sample.setdefault(row['sample_id'], set()).add(int(row['length']))
    fs_rows = read_tsv(data_dir / 'fs_region_full_length_labels.tsv.gz', compressed=True)
    for row in fs_rows:
        sample_id = row['endpoint_sample_id']
        if sample_id not in lengths_by_sample or len(lengths_by_sample[sample_id]) != 1:
            raise ValueError(f'FS-region endpoint identity is missing or ambiguous: {sample_id}')
        length = next(iter(lengths_by_sample[sample_id]))
        for column in ('sequence_position', 'start', 'end'):
            try:
                coordinate = int(row[column])
            except ValueError as exc:
                raise ValueError(f'non-integer FS-region coordinate: {sample_id}:{column}') from exc
            if not 1 <= coordinate <= length:
                raise ValueError(f'FS-region coordinate outside endpoint sequence: {sample_id}:{column}={coordinate}')
        if row['mask'] not in {'0', '1'}:
            raise ValueError(f"invalid FS-region mask: {row['mask']}")
        if row['mask'] == '0' and row['target_value'] != 'NA':
            raise ValueError('unknown FS-region value is not NA/masked')
        if row['mask'] == '1' and row['target_value'] not in {'0', '1'}:
            raise ValueError('active FS-region target is not binary')
    audited['fs_region_full_length_labels.tsv.gz'] = {'rows': len(fs_rows), 'coordinate_bounds_checked': len(fs_rows)}
    return audited

def verify_predictions(path: Path) -> dict:
    counts = {task: 0 for task in EXPECTED}
    seen = set()
    with gzip.open(path, 'rt', newline='', encoding='utf-8-sig') as handle:
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
    with gzip.open(sample_path, 'rt', newline='', encoding='utf-8-sig') as handle:
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
    with fold_path.open(newline='', encoding='utf-8-sig') as handle:
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
    with gzip.open(predictions_path, 'rt', newline='', encoding='utf-8-sig') as handle:
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
    catalog_path = root / 'data/task_catalog.tsv'
    with catalog_path.open(newline='', encoding='utf-8-sig') as handle:
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
    with gzip.open(coverage_path, 'rt', newline='', encoding='utf-8-sig') as handle:
        coverage_rows = list(csv.DictReader(handle, delimiter='\t'))
    coverage = {(row['task_id'], row['sample_id'], row['sequence_sha256']): row for row in coverage_rows}
    ready_counts = {}
    sample_path = data_dir / 'task_samples.tsv.gz'
    with gzip.open(sample_path, 'rt', newline='', encoding='utf-8-sig') as handle:
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
    with gzip.open(data_dir / 'task_samples.tsv.gz', 'rt', newline='', encoding='utf-8-sig') as handle:
        samples = list(csv.DictReader(handle, delimiter='\t'))
    eligible_groups = {task: {row['binding_group_id'] for row in samples if row['task_id'] == task and row['data_eligible'] == '1'} for task in scored_tasks}
    outer_path = data_dir / 'folds/outer_group_holdout.tsv'
    with outer_path.open(newline='', encoding='utf-8-sig') as handle:
        outer_rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['task_id'] in scored_tasks]
    inner_path = data_dir / 'folds/inner_validation_groups.tsv.gz'
    with gzip.open(inner_path, 'rt', newline='', encoding='utf-8-sig') as handle:
        inner_rows = [row for row in csv.DictReader(handle, delimiter='\t') if row['task_id'] in scored_tasks]
    fold_counts = {task: int(row['outer_folds']) for task, row in scored_tasks.items()}
    return audit_group_splits(eligible_groups, outer_rows, inner_rows, fold_counts)

def verify_case_source_references(data_dir: Path) -> dict:
    with (data_dir / 'case_modules.tsv').open(newline='', encoding='utf-8-sig') as handle:
        modules = list(csv.DictReader(handle, delimiter='\t'))
    with (data_dir / 'case_sources.tsv').open(newline='', encoding='utf-8-sig') as handle:
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

def verify_stable_controls(root):
    check_manifest(root)
    with (root / 'endpoint_provenance.tsv').open(newline='', encoding='utf-8-sig') as handle:
        endpoints = {r['endpoint_sample_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    counts = {}
    tables = {}
    for name in ['labels_primary.tsv.gz', 'labels_strict_sensitivity.tsv.gz']:
        seen = set()
        c = Counter()
        values = {}
        with gzip.open(root / name, 'rt', encoding='utf-8-sig', newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                sid = row['endpoint_sample_id']
                p = int(row['position_1based'])
                key = (sid, p)
                if sid not in endpoints or key in seen or (not 1 <= p <= int(endpoints[sid]['length'])):
                    raise ValueError('duplicate/out-of-range FS control identity')
                ep = endpoints[sid]
                for field in ['sequence_id', 'sequence_sha256', 'binding_group_id', 'case_id']:
                    if row[field] != ep[field]:
                        raise ValueError('FS control endpoint identity mismatch')
                if (row['mask'], row['label']) not in {('0', 'NA'), ('1', '0'), ('1', '1')}:
                    raise ValueError('unknown FS label not masked')
                seen.add(key)
                c[row['label']] += 1
                values[key] = row['label']
        for sid, ep in endpoints.items():
            if {p for s, p in seen if s == sid} != set(range(1, int(ep['length']) + 1)):
                raise ValueError('FS control endpoint coordinate coverage differs')
        counts[name] = dict(c)
        tables[name] = values
    primary = tables['labels_primary.tsv.gz']
    strict = tables['labels_strict_sensitivity.tsv.gz']
    if counts['labels_primary.tsv.gz'] != {'1': 124, '0': 848, 'NA': 4746}:
        raise ValueError('FS primary counts changed')
    if counts['labels_strict_sensitivity.tsv.gz'] != {'1': 124, '0': 580, 'NA': 5014}:
        raise ValueError('FS strict counts changed')
    if any((primary[k] != '0' for k, v in strict.items() if v == '0')) or {k for k, v in primary.items() if v == '1'} != {k for k, v in strict.items() if v == '1'}:
        raise ValueError('FS strict subset or inherited positives changed')
    return {'cases': 10, 'endpoints': len(endpoints), 'primary': counts['labels_primary.tsv.gz'], 'strict': counts['labels_strict_sensitivity.tsv.gz']}
