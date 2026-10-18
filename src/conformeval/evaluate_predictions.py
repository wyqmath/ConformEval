from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import platform
import re
import sys
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score
TASK_KNOT = 'KNOT-PRESENCE'
TASK_FS = 'FS-L2-grouped-PU-ranking'
TASKS = {TASK_KNOT, TASK_FS}
PREDICTION_FIELDS = ('task_id', 'sample_id', 'outer_fold', 'score')
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 2026
HASH_RE = re.compile('^[0-9a-f]{64}$')
ACCEPTED_PROTOCOL_SHA256 = '46cb535c76590af8fc29b13b1a9d5d4dbdd492ad9640470e11c0769d05ecdff5'

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()

def now() -> str:
    return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')

def read_table(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_gzip_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0, compresslevel=9) as zipped:
            with io.TextIOWrapper(zipped, encoding='utf-8', newline='') as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)

def verify_dataset(data_dir: Path) -> dict[str, str]:
    manifest_path = data_dir / 'dataset_manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    records = manifest.get('files', {})
    if not records:
        raise ValueError('dataset manifest has no file records')
    hashes = {'data/base/dataset_manifest.json': sha256(manifest_path)}
    for relative, record in records.items():
        rel = Path(relative)
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError(f'unsafe dataset manifest path: {relative}')
        path = data_dir / rel
        expected = record['sha256'] if isinstance(record, dict) else record
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f'dataset file hash mismatch: {relative}')
        if isinstance(record, dict) and path.stat().st_size != record.get('bytes'):
            raise ValueError(f'dataset file size mismatch: {relative}')
        hashes[f'data/base/{rel.as_posix()}'] = expected
    return hashes

def load_protocol(path: Path) -> dict:
    protocol = json.loads(path.read_text(encoding='utf-8-sig'))
    if protocol.get('schema_version') != 1 or protocol.get('protocol_id') != 'ConformEval-0.3-frozen-grouped-development':
        raise ValueError('unsupported scoring protocol')
    if protocol.get('source_accepted_protocol_sha256') != ACCEPTED_PROTOCOL_SHA256:
        raise ValueError('scoring protocol does not cite the accepted source protocol')
    tasks = protocol.get('tasks', {})
    expected = {TASK_KNOT: (5, 'pooled_outer_fold_OOF_AUROC', 'mean_inner_validation_AUROC'), TASK_FS: (10, 'mean_held_positive_rank_percentile_within_own_outer_test_pool_lower_is_better', 'mean_inner_held_positive_rank_percentile_within_own_validation_pool_lower_is_better')}
    if set(tasks) != TASKS:
        raise ValueError('scoring protocol task set differs from the released score-eligible tasks')
    for task, (outer, primary, inner) in expected.items():
        record = tasks[task]
        if record.get('outer_folds') != outer or record.get('inner_folds') != 3 or record.get('primary_metric') != primary or (record.get('inner_selection_metric') != inner):
            raise ValueError(f'task metric or fold protocol differs from the accepted lock: {task}')
    uncertainty = protocol.get('uncertainty', {})
    if uncertainty.get('bootstrap_replicates') != BOOTSTRAPS or uncertainty.get('seed') != BOOTSTRAP_SEED or uncertainty.get('confidence_level') != 0.95:
        raise ValueError('uncertainty protocol differs from the accepted lock')
    return protocol

def eligible_frames(data_dir: Path, tasks: set[str] | None = None) -> tuple[dict[tuple[str, str], dict], dict[tuple[str, str], str]]:
    tasks = TASKS if tasks is None else set(tasks)
    if not tasks or not tasks <= TASKS:
        raise ValueError('unsupported or empty task selection')
    samples = read_table(data_dir / 'task_samples.tsv.gz')
    eligible: dict[tuple[str, str], dict] = {}
    for row in samples:
        task = row['task_id']
        if task not in tasks or row['data_eligible'] != '1' or row['use_policy'] != 'development_only':
            continue
        if task == TASK_KNOT:
            if row['mask'] != '1' or row['label'] not in {'0', '1'}:
                continue
            if row['biological_mask'] != '1' or row['biological_target'] != row['label']:
                raise ValueError(f"KNOT-PRESENCE biological label/mask mismatch: {row['sample_id']}")
        else:
            if row['learning_mask'] != '1' or row['learning_label'] not in {'0', '1'}:
                continue
            if row['learning_label'] == '0':
                if row['biological_target'] != 'null' or row['biological_mask'] != '0' or row['mask'] != '0' or (row['label'] != 'NA'):
                    raise ValueError(f"FS-L2 unlabeled row was encoded as a biological negative: {row['sample_id']}")
            elif row['biological_target'] != 'confirmed_positive' or row['biological_mask'] != '1':
                raise ValueError(f"FS-L2 positive label/mask mismatch: {row['sample_id']}")
        key = (task, row['sample_id'])
        if key in eligible:
            raise ValueError(f'duplicate benchmark sample identity: {key}')
        eligible[key] = row
    fold_rows = read_table(data_dir / 'folds/outer_group_holdout.tsv')
    group_to_fold: dict[tuple[str, str], str] = {}
    for row in fold_rows:
        if row['task_id'] not in tasks:
            continue
        if row['role'] != 'test':
            raise ValueError(f'outer fold map contains a non-test row: {row}')
        key = (row['task_id'], row['group_id'])
        if key in group_to_fold:
            raise ValueError(f'binding group assigned to multiple outer folds: {key}')
        group_to_fold[key] = row['outer_fold']
    for task in tasks:
        frame = [row for key, row in eligible.items() if key[0] == task]
        if not frame:
            raise ValueError(f'no eligible samples for scored task {task}')
        groups = {row['binding_group_id'] for row in frame}
        if groups != {group for t, group in group_to_fold if t == task}:
            raise ValueError(f'outer-fold group coverage does not match eligible data: {task}')
        folds = {group_to_fold[task, group] for group in groups}
        expected_count = 5 if task == TASK_KNOT else 10
        if len(folds) != expected_count:
            raise ValueError(f'wrong frozen outer-fold count for {task}: {len(folds)}')
    audit_nested_group_maps(data_dir, eligible, group_to_fold)
    return (eligible, group_to_fold)

def audit_nested_group_maps(data_dir: Path, eligible: dict, group_to_fold: dict) -> None:
    tasks = {task for task, _ in eligible}
    inner_rows = read_table(data_dir / 'folds/inner_validation_groups.tsv.gz')
    qc_rows = read_table(data_dir / 'folds/fold_qc.tsv')
    inner_map: dict[tuple[str, str, str], set[str]] = {}
    for row in inner_rows:
        task, outer, inner, group = (row['task_id'], row['outer_fold'], row['inner_fold'], row['group_id'])
        if task not in tasks:
            continue
        if row['role'] != 'validation':
            raise ValueError(f'unexpected frozen inner-fold row: {row}')
        key = (task, outer, inner)
        bucket = inner_map.setdefault(key, set())
        if group in bucket:
            raise ValueError(f'duplicate group in inner validation: {key}/{group}')
        bucket.add(group)
    qc = {}
    for row in qc_rows:
        if row['task_id'] not in tasks:
            continue
        key = (row['task_id'], row['stage'], row['outer_fold'], row['inner_fold'])
        if key in qc or row['valid'] != '1':
            raise ValueError(f'invalid or duplicate frozen fold-QC row: {key}')
        qc[key] = row

    def compare_qc(task: str, stage: str, outer: str, inner: str, train_rows: list[dict], held_rows: list[dict], train_groups: set[str], held_groups: set[str]) -> None:
        key = (task, stage, outer, inner)
        record = qc.get(key)
        if record is None:
            raise ValueError(f'missing frozen fold-QC record: {key}')
        label_field = 'label' if task == TASK_KNOT else 'learning_label'
        train_classes = {label: sum((row[label_field] == label for row in train_rows)) for label in ('0', '1')}
        held_classes = {label: sum((row[label_field] == label for row in held_rows)) for label in ('0', '1')}
        actual = (len(train_rows), len(held_rows), len(train_groups), len(held_groups), train_classes, held_classes)
        recorded = (int(record['train_samples']), int(record['held_samples']), int(record['train_groups']), int(record['held_groups']), json.loads(record['train_classes']), json.loads(record['held_classes']))
        if actual != recorded:
            raise ValueError(f'fold-QC counts/classes differ from joined data: {key}')
        if set(train_classes) - {'0', '1'} or set(held_classes) - {'0', '1'} or (not all(train_classes.values())) or (not all(held_classes.values())):
            raise ValueError(f'frozen fold lacks both operational classes: {key}')
    expected_qc = set()
    for task in tasks:
        task_rows = [row for (row_task, _), row in eligible.items() if row_task == task]
        group_set = {row['binding_group_id'] for row in task_rows}
        n_outer = 5 if task == TASK_KNOT else 10
        n_inner = 3
        for fold_num in range(n_outer):
            outer = str(fold_num)
            test_groups = {group for (row_task, group), fold in group_to_fold.items() if row_task == task and fold == outer}
            train_groups = group_set - test_groups
            if not test_groups or test_groups & train_groups:
                raise ValueError(f'invalid outer group partition: {task}/{outer}')
            test_rows = [row for row in task_rows if row['binding_group_id'] in test_groups]
            train_rows = [row for row in task_rows if row['binding_group_id'] in train_groups]
            expected_qc.add((task, 'outer', outer, 'NA'))
            compare_qc(task, 'outer', outer, 'NA', train_rows, test_rows, train_groups, test_groups)
            validation_sets = []
            for inner_num in range(n_inner):
                inner = str(inner_num)
                key = (task, outer, inner)
                validation = inner_map.get(key)
                if not validation or validation & test_groups or (not validation.issubset(train_groups)):
                    raise ValueError(f'invalid inner validation group set: {key}')
                validation_sets.append(validation)
                val_rows = [row for row in train_rows if row['binding_group_id'] in validation]
                inner_train_groups = train_groups - validation
                inner_train_rows = [row for row in train_rows if row['binding_group_id'] in inner_train_groups]
                expected_qc.add((task, 'inner', outer, inner))
                compare_qc(task, 'inner', outer, inner, inner_train_rows, val_rows, inner_train_groups, validation)
            if any((left & right for i, left in enumerate(validation_sets) for right in validation_sets[i + 1:])):
                raise ValueError(f'inner validation groups overlap: {task}/{outer}')
            if set.union(*validation_sets) != train_groups:
                raise ValueError(f'inner validation groups do not partition outer training: {task}/{outer}')
    if set(qc) != expected_qc or set(inner_map) != {(task, str(outer), str(inner)) for task in tasks for outer in range(5 if task == TASK_KNOT else 10) for inner in range(3)}:
        raise ValueError('frozen fold-QC or inner-group map has unexpected or missing records')

def read_predictions(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if set(reader.fieldnames or ()) != set(PREDICTION_FIELDS):
            raise ValueError(f"prediction table must contain exactly these fields: {', '.join(PREDICTION_FIELDS)}")
        return list(reader)

def load_predictions(path: Path, eligible: dict, group_to_fold: dict) -> dict[tuple[str, str], float]:
    seen: dict[tuple[str, str], float] = {}
    for line_number, row in enumerate(read_predictions(path), start=2):
        key = (row.get('task_id', ''), row.get('sample_id', ''))
        if key[0] not in TASKS:
            raise ValueError(f'unknown or unsupported task at prediction line {line_number}: {key[0]}')
        if key not in eligible:
            raise ValueError(f'prediction sample is absent or ineligible: {key}')
        if key in seen:
            raise ValueError(f'duplicate prediction identity: {key}')
        expected_fold = group_to_fold.get((key[0], eligible[key]['binding_group_id']))
        if row.get('outer_fold') != expected_fold:
            raise ValueError(f'prediction fold differs from frozen outer assignment: {key}')
        try:
            score = float(row['score'])
        except (KeyError, ValueError) as exc:
            raise ValueError(f'prediction score is missing or invalid: {key}') from exc
        if not math.isfinite(score):
            raise ValueError(f'prediction score is non-finite: {key}')
        seen[key] = score
    if set(seen) != set(eligible):
        missing = sorted(set(eligible) - set(seen))[:10]
        extra = sorted(set(seen) - set(eligible))[:10]
        raise ValueError(f'prediction coverage differs from eligible data; missing={missing}; extra={extra}')
    return seen

def rank_positive(scores: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    if int(np.sum(labels == 1)) != 1:
        raise ValueError('each FS-L2 frozen outer fold must contain exactly one positive case')
    positive = float(scores[np.flatnonzero(labels == 1)[0]])
    greater = int(np.sum(scores > positive))
    equal = int(np.sum(scores == positive))
    rank = greater + (equal + 1) / 2
    recall = {f'recall_at_{k}': max(0.0, min(equal, k - greater)) / equal for k in (1, 10, 100, 1000, 10000)}
    return {'rank': float(rank), 'universe': int(len(scores)), 'percentile': rank / len(scores), **recall}

def group_bootstrap_auroc(rows: list[dict], score_map: dict, replicates: int, seed: int) -> dict:
    groups = sorted({row['binding_group_id'] for row in rows})
    members = {group: [row for row in rows if row['binding_group_id'] == group] for group in groups}
    rng = np.random.default_rng(seed)
    values = []
    for draw in rng.integers(0, len(groups), size=(replicates, len(groups))):
        sampled = [row for j in draw for row in members[groups[int(j)]]]
        labels = np.asarray([int(row['label']) for row in sampled], dtype=np.int8)
        if len(np.unique(labels)) < 2:
            continue
        scores = np.asarray([score_map[row['task_id'], row['sample_id']] for row in sampled], dtype=np.float64)
        values.append(float(roc_auc_score(labels, scores)))
    if not values:
        raise ValueError('group bootstrap produced no valid KNOT-PRESENCE replicates')
    return {'unit': 'binding_group_id', 'replicates': replicates, 'valid_replicates': len(values), 'seed': seed, 'conditional_on_fixed_outer_predictions': True, 'ci95': np.quantile(values, [0.025, 0.975]).tolist()}

def fs_metrics(rows: list[dict], score_map: dict, group_to_fold: dict, replicates: int, seed: int) -> dict:
    by_fold: dict[str, list[dict]] = {}
    for row in rows:
        fold = group_to_fold[TASK_FS, row['binding_group_id']]
        by_fold.setdefault(fold, []).append(row)
    ranks = []
    for fold, fold_rows in sorted(by_fold.items(), key=lambda item: int(item[0])):
        ordered = sorted(fold_rows, key=lambda row: row['sample_id'])
        scores = np.asarray([score_map[TASK_FS, row['sample_id']] for row in ordered], dtype=np.float64)
        labels = np.asarray([int(row['learning_label']) for row in ordered], dtype=np.int8)
        rank = rank_positive(scores, labels)
        positive = ordered[int(np.flatnonzero(labels == 1)[0])]
        ranks.append({'outer_fold': fold, 'sample_id': positive['sample_id'], 'binding_group_id': positive['binding_group_id'], **rank})
    if len(ranks) != 10 or len({row['binding_group_id'] for row in ranks}) != 10:
        raise ValueError('FS-L2 positive case coverage differs from its frozen ten outer units')
    values = np.asarray([row['percentile'] for row in ranks], dtype=np.float64)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(ranks), size=(replicates, len(ranks)))
    means = values[draws].mean(axis=1)
    summary = {'primary_metric': 'mean held-positive rank percentile within its own frozen outer-fold test pool; lower is better', 'mean_held_case_rank_percentile': float(values.mean()), 'mean_recall_at_1': float(np.mean([row['recall_at_1'] for row in ranks])), 'mean_recall_at_10': float(np.mean([row['recall_at_10'] for row in ranks])), 'mean_recall_at_100': float(np.mean([row['recall_at_100'] for row in ranks])), 'mean_recall_at_1000': float(np.mean([row['recall_at_1000'] for row in ranks])), 'mean_recall_at_10000': float(np.mean([row['recall_at_10000'] for row in ranks])), 'held_positive_cases': len(ranks), 'unlabeled_background_is_not_a_verified_negative': True, 'group_bootstrap_ci95': np.quantile(means, [0.025, 0.975]).tolist(), 'bootstrap_unit': 'positive binding_group_id', 'bootstrap_replicates': replicates, 'bootstrap_seed': seed, 'bootstrap_scope': 'conditional on fixed outer-fold predictions; no refitting uncertainty claim', 'case_ranks': ranks}
    return summary

def evaluate(data_dir: Path, predictions_path: Path, metadata_path: Path, protocol_path: Path, output_dir: Path, tasks: set[str] | None = None) -> dict:
    tasks = TASKS if tasks is None else set(tasks)
    if output_dir.exists():
        raise ValueError(f'refusing to overwrite existing output directory: {output_dir}')
    for path in (data_dir, predictions_path, metadata_path, protocol_path):
        if not path.exists():
            raise FileNotFoundError(path)
    protocol_hash = sha256(protocol_path)
    protocol = load_protocol(protocol_path)
    metadata_hash = sha256(metadata_path)
    predictions_hash = sha256(predictions_path)
    metadata = json.loads(metadata_path.read_text(encoding='utf-8-sig'))
    for key in ('model_id', 'model_revision', 'prediction_method', 'prediction_code_sha256', 'training_protocol'):
        if not isinstance(metadata.get(key), str) or not metadata[key].strip():
            raise ValueError(f'model metadata is missing required field: {key}')
    if not HASH_RE.fullmatch(metadata['prediction_code_sha256']):
        raise ValueError('prediction_code_sha256 must be a lowercase SHA-256 digest')
    candidate_layers = metadata.get('candidate_layers')
    if not isinstance(candidate_layers, list) or not candidate_layers or any((type(layer) is not int or layer < 1 for layer in candidate_layers)) or (len(candidate_layers) != len(set(candidate_layers))):
        raise ValueError('candidate_layers must contain unique positive 1-based layer IDs')
    if not isinstance(metadata.get('pooling'), str) or not metadata['pooling'].strip():
        raise ValueError('model metadata is missing pooling')
    if metadata['training_protocol'] != 'frozen_nested_grouped_oof_v1':
        raise ValueError('predictions must declare the frozen nested grouped out-of-fold protocol')
    started = now()
    code_path = Path(__file__).resolve()
    code_hash = sha256(code_path)
    input_hashes = verify_dataset(data_dir)
    dataset_metadata = json.loads((data_dir / 'dataset_manifest.json').read_text(encoding='utf-8-sig'))
    input_hashes['prediction_file'] = predictions_hash
    input_hashes['model_metadata'] = metadata_hash
    input_hashes['scoring_protocol'] = protocol_hash
    eligible, group_to_fold = eligible_frames(data_dir, tasks)
    prediction_map = load_predictions(predictions_path, eligible, group_to_fold)
    scored_rows = []
    metrics: dict[str, dict] = {}
    for task in sorted(tasks):
        rows = sorted((row for (row_task, _), row in eligible.items() if row_task == task), key=lambda row: row['sample_id'])
        if task == TASK_KNOT:
            labels = np.asarray([int(row['label']) for row in rows], dtype=np.int8)
            scores = np.asarray([prediction_map[task, row['sample_id']] for row in rows], dtype=np.float64)
            summary = {'primary_metric': 'AUROC from complete frozen outer-fold out-of-fold predictions', 'auroc': float(roc_auc_score(labels, scores)), 'average_precision_secondary': float(average_precision_score(labels, scores)), 'n_samples': int(len(rows)), 'n_positive': int(labels.sum()), 'n_binding_groups': len({row['binding_group_id'] for row in rows}), 'group_bootstrap': group_bootstrap_auroc(rows, prediction_map, BOOTSTRAPS, BOOTSTRAP_SEED)}
        else:
            summary = fs_metrics(rows, prediction_map, group_to_fold, BOOTSTRAPS, BOOTSTRAP_SEED)
            summary['n_evaluation_units'] = len(rows)
            summary['n_binding_groups'] = len({row['binding_group_id'] for row in rows})
        metrics[task] = summary
        for row in rows:
            fold = group_to_fold[task, row['binding_group_id']]
            scored_rows.append({'task_id': task, 'sample_id': row['sample_id'], 'outer_fold': fold, 'binding_group_id': row['binding_group_id'], 'score': prediction_map[task, row['sample_id']], 'learning_label': row['learning_label'], 'learning_mask': row['learning_mask'], 'biological_target': row['biological_target'], 'biological_mask': row['biological_mask'], 'label': row['label'], 'mask': row['mask']})
    if verify_dataset(data_dir) != {key: value for key, value in input_hashes.items() if key.startswith('data/')}:
        raise ValueError('benchmark data changed during evaluation')
    if sha256(predictions_path) != predictions_hash or sha256(metadata_path) != metadata_hash or sha256(protocol_path) != protocol_hash:
        raise ValueError('prediction, model metadata, or scoring protocol changed during evaluation')
    output_dir.mkdir(parents=True)
    scored_path = output_dir / 'scored_predictions.tsv.gz'
    write_gzip_tsv(scored_path, scored_rows, ['task_id', 'sample_id', 'outer_fold', 'binding_group_id', 'score', 'learning_label', 'learning_mask', 'biological_target', 'biological_mask', 'label', 'mask'])
    metrics_path = output_dir / 'metrics.json'
    metrics_path.write_text(json.dumps({'model': metadata, 'metrics': metrics}, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')
    if sha256(code_path) != code_hash:
        raise ValueError('scoring code changed during evaluation')
    output_hashes = {path.name: sha256(path) for path in (scored_path, metrics_path)}
    run = {'status': 'PASS', 'started_at': started, 'completed_at': now(), 'scope': 'frozen retrospective grouped-development predictions; model training and selection history are self-reported and not verified by this scorer', 'benchmark_data_version': '2026-10-14-data-lock-v01+FD.03', 'scoring_protocol_id': protocol['protocol_id'], 'input_sha256': input_hashes, 'scorer_sha256': code_hash, 'model_metadata': metadata, 'eligible_sample_counts': {task: sum((key[0] == task for key in eligible)) for task in sorted(tasks)}, 'prediction_counts': {task: sum((key[0] == task for key in prediction_map)) for task in sorted(tasks)}, 'outer_folds': {task: len({group_to_fold[task, row['binding_group_id']] for key, row in eligible.items() if key[0] == task}) for task in sorted(tasks)}, 'scorer_uses_outer_test_labels_for_model_selection': False, 'model_training_protocol_self_attested': True, 'warnings': [], 'dependencies': {'python': platform.python_version(), 'numpy': np.__version__, 'scikit_learn': version('scikit-learn')}, 'output_sha256': output_hashes}
    run_path = output_dir / 'run.json'
    run['project'] = 'ConformEval'
    if dataset_metadata.get('scope') == 'software_example':
        run['scope'] = 'software_example'
        run['benchmark_data_version'] = dataset_metadata['artifact_id']
    run_path.write_text(json.dumps(run, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')
    return run

def main() -> int:
    parser = argparse.ArgumentParser(description='Score complete out-of-fold predictions against the frozen ConformEval labels.')
    parser.add_argument('--example', action='store_true', help='Run the bundled synthetic knot example from the repository root')
    parser.add_argument('--data-dir', type=Path)
    parser.add_argument('--task', action='append', choices=sorted(TASKS), help='Repeat to select pooled tasks; default: both')
    parser.add_argument('--protocol', type=Path)
    parser.add_argument('--predictions', type=Path, help='TSV/TSV.GZ with task_id, sample_id, outer_fold, score')
    parser.add_argument('--model-metadata', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    inputs = ('data_dir', 'task', 'protocol', 'predictions', 'model_metadata')
    if args.example:
        supplied = ['--' + name.replace('_', '-') for name in inputs if getattr(args, name) is not None]
        if supplied:
            parser.error('--example cannot be combined with ' + ', '.join(supplied))
        args.data_dir = Path('examples/mini/base')
        args.task = [TASK_KNOT]
        args.protocol = Path('configs/scoring.json')
        args.predictions = Path('examples/predictions/knot.tsv')
        args.model_metadata = Path('examples/predictions/metadata.json')
        if not all(path.exists() for path in (args.data_dir, args.protocol, args.predictions, args.model_metadata)):
            parser.error('--example requires the bundled files; run it from the ConformEval repository root')
    else:
        missing = ['--' + name.replace('_', '-') for name in inputs if name != 'task' and getattr(args, name) is None]
        if missing:
            parser.error('the following arguments are required without --example: ' + ', '.join(missing))
    try:
        result = evaluate(args.data_dir.resolve(), args.predictions.resolve(), args.model_metadata.resolve(), args.protocol.resolve(), args.output_dir.resolve(), set(args.task) if args.task else None)
    except (ValueError, KeyError, OSError, TypeError, json.JSONDecodeError) as error:
        print(f'Evaluation rejected: {error}', file=sys.stderr)
        return 2
    print(json.dumps({'status': result['status'], 'metrics': 'metrics.json', 'output_dir': str(args.output_dir)}, ensure_ascii=False))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
