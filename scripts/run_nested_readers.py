from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import sys
import warnings
from collections import Counter, defaultdict
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import scipy
from scipy.spatial.distance import cdist
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
try:
    from .evaluate_predictions import TASK_FS, TASK_KNOT, eligible_frames, evaluate, load_protocol, read_table, verify_dataset, sha256, write_gzip_tsv
    from .physinfobench_model_interface import load_mean_representation
except ImportError:
    from evaluate_predictions import TASK_FS, TASK_KNOT, eligible_frames, evaluate, load_protocol, read_table, verify_dataset, sha256, write_gzip_tsv
    from physinfobench_model_interface import load_mean_representation
LAYERS = (17, 33)
C_GRID = (0.1, 1.0)
SEED = 2026
MAX_ITER = 5000
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
ARMS = ('PLM-nested', 'PLM-hidden17', 'PLM-hidden33', 'PLM-length', 'PLM-composition', 'PLM-nearest', 'PLM-training-prior')

def stamp() -> str:
    return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')

def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    for raw in path.read_text(encoding='utf-8').splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            current = line[1:].split()[0]
            if current in records:
                raise ValueError(f'duplicate FASTA identifier: {current}')
            records[current] = []
        elif current is None:
            raise ValueError('FASTA sequence precedes its header')
        else:
            records[current].append(line.upper())
    sequences = {key: ''.join(parts) for key, parts in records.items()}
    if not sequences or any((not value for value in sequences.values())):
        raise ValueError('FASTA is empty or contains an empty sequence')
    return sequences

def sequence_aux(sequence: str) -> np.ndarray:
    length = len(sequence)
    counts = Counter(sequence)
    return np.asarray([math.log1p(length)] + [counts[aa] / length for aa in AA20], dtype=np.float32)

def safe_archive_path(root: Path, relative: str) -> Path:
    rel = Path(relative)
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError(f'unsafe representation archive path: {relative}')
    path = (root / rel).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError(f'representation archive is missing or escapes its root: {relative}')
    return path

def load_representation_manifest(path: Path, data_dir: Path, eligible: dict):
    root = path.parent.resolve()
    manifest = json.loads(path.read_text(encoding='utf-8'))
    for field in ('model_id', 'model_revision', 'pooling', 'extractor_code_sha256', 'candidate_layers', 'archives'):
        if field not in manifest:
            raise ValueError(f'representation manifest lacks {field}')
    if any((not isinstance(manifest[field], str) or not manifest[field].strip() for field in ('model_id', 'model_revision', 'pooling'))):
        raise ValueError('representation model identity and pooling must be nonempty strings')
    if not isinstance(manifest['extractor_code_sha256'], str) or len(manifest['extractor_code_sha256']) != 64:
        raise ValueError('extractor_code_sha256 must be a SHA-256 string')
    try:
        int(manifest['extractor_code_sha256'], 16)
    except ValueError as exc:
        raise ValueError('extractor_code_sha256 must be hexadecimal') from exc
    if manifest['candidate_layers'] != list(LAYERS):
        raise ValueError('current frozen protocol requires 1-based candidate layers [17, 33]')
    if not isinstance(manifest['archives'], list) or not manifest['archives']:
        raise ValueError('representation manifest archives must be a nonempty list')
    sequences = parse_fasta(data_dir / 'sequences.fasta')
    sample_sequences = read_table(data_dir / 'sample_sequences.tsv.gz')
    used_keys = set(eligible)
    sequence_rows: dict[tuple[str, str], list[dict]] = defaultdict(list)
    required: dict[str, tuple[str, int]] = {}
    for row in sample_sequences:
        key = (row['task_id'], row['sample_id'])
        if key not in used_keys:
            continue
        digest = row['sequence_sha256']
        length = int(row['length'])
        sequence = sequences.get(row['sequence_id'])
        if sequence is None or len(sequence) != length:
            raise ValueError(f"sample sequence is absent or has wrong length: {key}/{row['sequence_id']}")
        actual_digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        if actual_digest != digest or row['sequence_id'] != f'SEQ_{digest}':
            raise ValueError(f"sample sequence identity differs from FASTA: {key}/{row['sequence_id']}")
        old = required.setdefault(digest, (row['sequence_id'], length))
        if old != (row['sequence_id'], length):
            raise ValueError(f'one sequence hash maps to inconsistent IDs or lengths: {digest}')
        sequence_rows[key].append(row)
    archive_rows = {}
    for row in manifest['archives']:
        digest = row.get('sequence_sha256')
        length = row.get('length')
        file_name = row.get('file')
        expected_hash = row.get('archive_sha256')
        if not isinstance(digest, str) or digest in archive_rows or digest not in required:
            raise ValueError(f'duplicate or unexpected sequence in representation manifest: {digest}')
        if type(length) is not int or length != required[digest][1]:
            raise ValueError(f'representation manifest length differs for {digest}')
        if not isinstance(expected_hash, str) or len(expected_hash) != 64:
            raise ValueError(f'invalid archive SHA-256 in representation manifest: {digest}')
        try:
            int(expected_hash, 16)
        except ValueError as exc:
            raise ValueError(f'invalid archive SHA-256 in representation manifest: {digest}') from exc
        archive_path = safe_archive_path(root, file_name)
        if sha256(archive_path) != expected_hash:
            raise ValueError(f'representation archive hash mismatch: {file_name}')
        selected, dimension = load_mean_representation(archive_path, digest, length, expected_hash, LAYERS, manifest['model_id'], manifest['model_revision'])
        if selected.shape != (len(LAYERS), dimension) or not np.isfinite(selected).all():
            raise ValueError(f'invalid selected representation array: {file_name}')
        archive_rows[digest] = {'path': archive_path, 'sha256': expected_hash, 'features': selected, 'dimension': dimension, 'length': length}
    if set(archive_rows) != set(required):
        missing = sorted(set(required) - set(archive_rows))[:10]
        raise ValueError(f'representation coverage differs from eligible sequence universe; missing={missing}')
    feature_by_sample = {TASK_KNOT: {}, TASK_FS: {}}
    auxiliary_by_sample = {TASK_KNOT: {}, TASK_FS: {}}
    task_rows = {TASK_KNOT: [], TASK_FS: []}
    for key, sample in eligible.items():
        task, sample_id = key
        mappings = sequence_rows.get(key, [])
        hashes, lengths = validate_sample_sequence_mappings(key, sample, mappings)
        vectors, auxiliary = ([], [])
        for digest, length in zip(hashes, lengths):
            source = archive_rows[digest]
            vectors.append(source['features'])
            seq_id = required[digest][0]
            auxiliary.append(sequence_aux(sequences[seq_id]))
        feature_by_sample[task][sample_id] = np.mean(np.stack(vectors), axis=0, dtype=np.float32)
        auxiliary_by_sample[task][sample_id] = np.mean(np.stack(auxiliary), axis=0, dtype=np.float32)
        task_rows[task].append(sample)
    for task in task_rows:
        task_rows[task].sort(key=lambda row: row['sample_id'])
        if any((feature_by_sample[task][row['sample_id']].shape[0] != len(LAYERS) for row in task_rows[task])):
            raise ValueError(f'candidate layer coverage differs within {task}')
    input_hashes = {'feature_manifest': sha256(path)}
    input_hashes.update({f'representation/{digest}': item['sha256'] for digest, item in sorted(archive_rows.items())})
    return (manifest, archive_rows, feature_by_sample, auxiliary_by_sample, task_rows, input_hashes)

def validate_sample_sequence_mappings(key: tuple[str, str], sample: dict, mappings: list[dict]) -> tuple[list[str], list[int]]:
    try:
        hashes = sample['sequence_sha256s'].split(';')
        sequence_ids = sample['sequence_ids'].split(';')
        lengths = [int(value) for value in sample['lengths'].split(';')]
    except (KeyError, ValueError) as exc:
        raise ValueError(f'eligible sample has invalid endpoint identity fields: {key}') from exc
    if not hashes or len(hashes) != len(lengths) or len(hashes) != len(sequence_ids):
        raise ValueError(f'eligible sample has invalid endpoint identity fields: {key}')
    expected = Counter(zip(sequence_ids, hashes, lengths))
    observed = Counter(((row.get('sequence_id'), row.get('sequence_sha256'), int(row['length'])) for row in mappings))
    if expected != observed:
        raise ValueError(f'eligible sample sequence mapping does not close: {key}')
    if key[0] == TASK_FS and sample.get('unit_type') == 'positive_endpoint_pair':
        roles = Counter((row.get('endpoint_role') for row in mappings))
        if roles != Counter({'endpoint_A': 1, 'endpoint_B': 1}):
            raise ValueError(f'fold-switch pair must map one endpoint_A and one endpoint_B: {key}')
    return (hashes, lengths)

def fold_indices(task: str, rows: list[dict], data_dir: Path, group_to_fold: dict):
    inner_rows = [row for row in read_table(data_dir / 'folds/inner_validation_groups.tsv.gz') if row['task_id'] == task]
    inner_group = defaultdict(dict)
    for row in inner_rows:
        outer, inner, group = (row['outer_fold'], row['inner_fold'], row['group_id'])
        if group in inner_group[outer]:
            raise ValueError(f'duplicate inner group assignment: {task}/{outer}/{group}')
        inner_group[outer][group] = inner
    y = np.asarray([int(row['label'] if task == TASK_KNOT else row['learning_label']) for row in rows], dtype=np.int8)
    groups = {row['binding_group_id'] for row in rows}
    n_outer = 5 if task == TASK_KNOT else 10
    plan = {}
    for outer_index in range(n_outer):
        outer = str(outer_index)
        test_groups = {group for group in groups if group_to_fold[task, group] == outer}
        train_groups = groups - test_groups
        inner_assignment = inner_group.get(outer, {})
        if set(inner_assignment) != train_groups:
            raise ValueError(f'inner validation groups do not partition outer training: {task}/{outer}')
        train = np.asarray([i for i, row in enumerate(rows) if row['binding_group_id'] in train_groups], dtype=np.int64)
        test = np.asarray([i for i, row in enumerate(rows) if row['binding_group_id'] in test_groups], dtype=np.int64)
        if set(y[train]) != {0, 1} or set(y[test]) != {0, 1}:
            raise ValueError(f'outer split lacks both operational classes: {task}/{outer}')
        inners = []
        for inner_index in range(3):
            inner = str(inner_index)
            val_groups = {group for group, fold in inner_assignment.items() if fold == inner}
            val = np.asarray([i for i, row in enumerate(rows) if row['binding_group_id'] in val_groups], dtype=np.int64)
            inner_train = np.asarray([i for i, row in enumerate(rows) if row['binding_group_id'] in train_groups - val_groups], dtype=np.int64)
            if not len(val) or set(y[val]) != {0, 1} or set(y[inner_train]) != {0, 1}:
                raise ValueError(f'inner split lacks both operational classes: {task}/{outer}/{inner}')
            inners.append((inner_index, inner_train, val))
        plan[outer_index] = (train, test, inners)
    return (y, plan)

def fs_rank_percentiles(scores: np.ndarray, labels: np.ndarray) -> np.ndarray:
    positives = np.flatnonzero(labels == 1)
    if not len(positives):
        raise ValueError('PU validation fold contains no held positive')
    output = []
    for index in positives:
        score = scores[index]
        greater = int(np.sum(scores > score))
        equal = int(np.sum(scores == score))
        output.append((greater + (equal + 1) / 2) / len(scores))
    return np.asarray(output, dtype=np.float64)

def selection_metric(task: str, y: np.ndarray, indices: np.ndarray, scores: np.ndarray) -> float:
    if task == TASK_KNOT:
        return float(roc_auc_score(y[indices], scores))
    return float(fs_rank_percentiles(scores, y[indices]).mean())

def fit_logistic(features: np.ndarray, labels: np.ndarray, train: np.ndarray, test: np.ndarray, C: float, context: dict, warnings_out: list[dict]):
    if set(labels[train]) != {0, 1}:
        raise ValueError(f'reader training partition lacks both classes: {context}')
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(features[train])
    test_scaled = scaler.transform(features[test])
    model = LogisticRegression(C=C, solver='liblinear', class_weight='balanced', max_iter=MAX_ITER, random_state=SEED)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.fit(train_scaled, labels[train])
    for warning in caught:
        warnings_out.append({**context, 'type': warning.category.__name__, 'message': str(warning.message)})
        if issubclass(warning.category, ConvergenceWarning):
            raise ValueError(f'reader did not converge: {context}')
    if list(model.classes_) != [0, 1]:
        raise ValueError(f'unexpected fitted class order: {context}')
    scores = model.predict_proba(test_scaled)[:, 1]
    if not np.isfinite(scores).all():
        raise ValueError(f'nonfinite reader prediction: {context}')
    params = {'scaler_mean': scaler.mean_.copy(), 'scaler_scale': scaler.scale_.copy(), 'coef': model.coef_.copy(), 'intercept': model.intercept_.copy(), 'classes': model.classes_.copy(), 'C': np.asarray(C)}
    return (scores, params)

def choose_hyperparameters(task: str, layer_features: dict[int, np.ndarray], length: np.ndarray, composition: np.ndarray, y: np.ndarray, inners: list, outer_fold: int, warnings_out: list, selection_rows: list) -> dict:
    candidates = [(f'PLM-hidden{layer}', layer, C, layer_features[layer]) for layer in LAYERS for C in C_GRID]
    candidates.extend(((arm, 0, C, X) for arm, X in (('PLM-length', length), ('PLM-composition', composition)) for C in C_GRID))
    values = defaultdict(list)
    for arm, layer, C, X in candidates:
        fold_values = []
        for inner, train, val in inners:
            context = {'task_id': task, 'outer_fold': outer_fold, 'inner_fold': inner, 'arm': arm, 'layer': layer, 'C': C, 'phase': 'inner'}
            scores, _ = fit_logistic(X, y, train, val, C, context, warnings_out)
            metric = selection_metric(task, y, val, scores)
            fold_values.append(metric)
            selection_rows.append({**context, 'selection_metric': metric})
        values[arm, layer, C] = float(np.mean(fold_values))
    maximize = task == TASK_KNOT
    key = lambda score, C, layer=0: (-score if maximize else score, C, layer)
    plm_options = [(score, C, layer) for (arm, layer, C), score in values.items() if arm.startswith('PLM-hidden')]
    nested = sorted(plm_options, key=lambda item: key(*item))[0]
    selected = {'PLM-nested': (nested[2], nested[1], nested[0])}
    for layer in LAYERS:
        options = [(score, C) for (arm, candidate_layer, C), score in values.items() if arm == f'PLM-hidden{layer}' and candidate_layer == layer]
        best = sorted(options, key=lambda item: key(item[0], item[1]))[0]
        selected[f'PLM-hidden{layer}'] = (layer, best[1], best[0])
    for arm in ('PLM-length', 'PLM-composition'):
        options = [(score, C) for (candidate_arm, _, C), score in values.items() if candidate_arm == arm]
        best = sorted(options, key=lambda item: key(item[0], item[1]))[0]
        selected[arm] = (None, best[1], best[0])
    return selected

def nearest_neighbor(X: np.ndarray, y: np.ndarray, train: np.ndarray, test: np.ndarray, sample_ids: list[str]):
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(X[train])
    test_scaled = scaler.transform(X[test])
    order = np.argsort(np.asarray([sample_ids[i] for i in train], dtype=str), kind='stable')
    train_sorted = train[order]
    train_sorted_features = train_scaled[order]
    distances = cdist(test_scaled, train_sorted_features, metric='euclidean')
    if not np.isfinite(distances).all():
        raise ValueError('nonfinite nearest-neighbor distance')
    nearest = np.argmin(distances, axis=1)
    return (y[train_sorted[nearest]].astype(np.float64), {'reference_sample_ids': [sample_ids[i] for i in train_sorted], 'reference_labels': y[train_sorted].astype(int).tolist(), 'scaler_mean': scaler.mean_.tolist(), 'scaler_scale': scaler.scale_.tolist()})

def run_task(task: str, rows: list[dict], features: dict[str, np.ndarray], auxiliary: dict[str, np.ndarray], fold_plan: dict, output_dir: Path, warnings_out: list, selection_rows: list, outer_rows: list, fold_rows: list, nearest_rows: list) -> dict[str, list[dict]]:
    sample_ids = [row['sample_id'] for row in rows]
    y = np.asarray([int(row['label'] if task == TASK_KNOT else row['learning_label']) for row in rows], dtype=np.int8)
    layers = {layer: np.vstack([features[sample][LAYERS.index(layer)] for sample in sample_ids]) for layer in LAYERS}
    aux = np.vstack([auxiliary[sample] for sample in sample_ids])
    length = aux[:, :1]
    composition = aux
    predictions = {arm: np.full(len(rows), np.nan, dtype=np.float64) for arm in ARMS}
    parameter_dir = output_dir / 'parameters' / task
    parameter_dir.mkdir(parents=True, exist_ok=True)
    fold_results = {}
    for outer_fold, (train, test, inners) in fold_plan.items():
        selected = choose_hyperparameters(task, layers, length, composition, y, inners, outer_fold, warnings_out, selection_rows)
        fold_results[outer_fold] = {'train_samples': len(train), 'test_samples': len(test), 'selected': selected}
        arm_features = {'PLM-nested': layers[selected['PLM-nested'][0]], 'PLM-hidden17': layers[17], 'PLM-hidden33': layers[33], 'PLM-length': length, 'PLM-composition': composition}
        for arm, X in arm_features.items():
            layer, C, metric = selected[arm]
            context = {'task_id': task, 'outer_fold': outer_fold, 'arm': arm, 'layer': layer if layer is not None else 'NA', 'C': C, 'phase': 'outer'}
            scores, params = fit_logistic(X, y, train, test, C, context, warnings_out)
            predictions[arm][test] = scores
            np.savez_compressed(parameter_dir / f'fold_{outer_fold:02d}_{arm}.npz', **params)
            outer_rows.append({'task_id': task, 'outer_fold': outer_fold, 'arm': arm, 'selected_layer': layer if layer is not None else 'NA', 'selected_C': C, 'inner_selection_metric': metric})
            fold_rows.append({'task_id': task, 'outer_fold': outer_fold, 'arm': arm, 'train_samples': len(train), 'test_samples': len(test)})
        near, reference = nearest_neighbor(layers[33], y, train, test, sample_ids)
        predictions['PLM-nearest'][test] = near
        prior = float(y[train].mean())
        predictions['PLM-training-prior'][test] = prior
        nearest_rows.append({'task_id': task, 'outer_fold': outer_fold, **reference})
        for arm in ('PLM-nearest', 'PLM-training-prior'):
            fold_rows.append({'task_id': task, 'outer_fold': outer_fold, 'arm': arm, 'train_samples': len(train), 'test_samples': len(test)})
    for arm, values in predictions.items():
        if not np.isfinite(values).all():
            raise ValueError(f'out-of-fold predictions are incomplete: {task}/{arm}')
    return {arm: [{'task_id': task, 'sample_id': sample_id, 'outer_fold': str(next((fold for fold, (_, test, _) in fold_plan.items() if index in test))), 'score': float(predictions[arm][index])} for index, sample_id in enumerate(sample_ids)] for arm in ARMS}

def _execute(args, output_dir: Path, started_at: str) -> dict:
    data_dir = args.data_dir.resolve()
    protocol_path = args.protocol.resolve()
    manifest_path = args.feature_manifest.resolve()
    if not data_dir.is_dir():
        raise FileNotFoundError(data_dir)
    protocol = load_protocol(protocol_path)
    eligible, group_to_fold = eligible_frames(data_dir)
    manifest, archive_rows, feature_by_sample, auxiliary_by_sample, task_rows, feature_hashes = load_representation_manifest(manifest_path, data_dir, eligible)
    script_path = Path(__file__).resolve()
    code_hash = sha256(script_path)
    input_hashes = {f'data/{key}': value for key, value in verify_dataset(data_dir).items()}
    input_hashes.update(feature_hashes)
    input_hashes['scoring_protocol'] = sha256(protocol_path)
    input_hashes['nested_reader_code'] = code_hash
    input_hashes['prediction_scorer_code'] = sha256(Path(__file__).with_name('evaluate_predictions.py'))
    input_hashes['model_interface_code'] = sha256(Path(__file__).with_name('physinfobench_model_interface.py'))
    all_predictions = {arm: [] for arm in ARMS}
    selection_rows, outer_rows, fold_rows, nearest_rows, warnings_out = ([], [], [], [], [])
    for task in (TASK_KNOT, TASK_FS):
        rows = task_rows[task]
        _, plan = fold_indices(task, rows, data_dir, group_to_fold)
        predictions = run_task(task, rows, feature_by_sample[task], auxiliary_by_sample[task], plan, output_dir, warnings_out, selection_rows, outer_rows, fold_rows, nearest_rows)
        for arm in ARMS:
            all_predictions[arm].extend(predictions[arm])
    write_fields = ['task_id', 'sample_id', 'outer_fold', 'score']
    scorer_results = {}
    model_manifest_copy = {key: value for key, value in manifest.items() if key != 'archives'}
    for arm in ARMS:
        prediction_path = output_dir / f'predictions_{arm}.tsv.gz'
        rows = sorted(all_predictions[arm], key=lambda row: (row['task_id'], row['sample_id']))
        write_gzip_tsv(prediction_path, rows, write_fields)
        method = 'nested_balanced_logistic_regression' if arm.startswith('PLM-') and arm not in {'PLM-nearest', 'PLM-training-prior'} else arm
        metadata = {'model_id': manifest['model_id'], 'model_revision': manifest['model_revision'], 'prediction_method': method, 'prediction_code_sha256': code_hash, 'training_protocol': 'frozen_nested_grouped_oof_v1', 'candidate_layers': list(LAYERS), 'pooling': manifest['pooling'], 'arm': arm, 'feature_extractor_code_sha256': manifest['extractor_code_sha256']}
        metadata_path = output_dir / f'metadata_{arm}.json'
        metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        scorer_out = output_dir / 'scores' / arm
        scorer_results[arm] = evaluate(data_dir, prediction_path, metadata_path, protocol_path, scorer_out)
    selection_fields = ['task_id', 'outer_fold', 'inner_fold', 'arm', 'layer', 'C', 'phase', 'selection_metric']
    outer_fields = ['task_id', 'outer_fold', 'arm', 'selected_layer', 'selected_C', 'inner_selection_metric']
    fold_fields = ['task_id', 'outer_fold', 'arm', 'train_samples', 'test_samples']
    write_gzip_tsv(output_dir / 'inner_selection.tsv.gz', selection_rows, selection_fields)
    write_gzip_tsv(output_dir / 'outer_selection.tsv.gz', outer_rows, outer_fields)
    write_gzip_tsv(output_dir / 'fold_metrics.tsv.gz', fold_rows, fold_fields)
    (output_dir / 'nearest_references.json').write_text(json.dumps(nearest_rows, indent=2) + '\n', encoding='utf-8')
    (output_dir / 'warnings.json').write_text(json.dumps(warnings_out, indent=2) + '\n', encoding='utf-8')
    expected_data = {f'data/{key}': value for key, value in verify_dataset(data_dir).items()}
    if {key: input_hashes[key] for key in expected_data} != expected_data:
        raise ValueError('benchmark data changed during nested evaluation')
    if sha256(manifest_path) != feature_hashes['feature_manifest'] or sha256(protocol_path) != input_hashes['scoring_protocol']:
        raise ValueError('feature manifest or scoring protocol changed during nested evaluation')
    for digest, item in archive_rows.items():
        if sha256(item['path']) != item['sha256']:
            raise ValueError(f'representation archive changed during evaluation: {digest}')
    if sha256(script_path) != code_hash:
        raise ValueError('nested reader code changed during evaluation')
    output_hashes = {str(path.relative_to(output_dir)): sha256(path) for path in sorted(output_dir.rglob('*')) if path.is_file() and path != output_dir / 'run.json'}
    result = {'status': 'PASS', 'started_at': started_at, 'completed_at': stamp(), 'scope': 'frozen retrospective grouped-development nested readers; no prospective or family-independent claim', 'protocol_id': protocol['protocol_id'], 'model': model_manifest_copy, 'candidate_layers': list(LAYERS), 'C_grid': list(C_GRID), 'outer_folds': {TASK_KNOT: 5, TASK_FS: 10}, 'inner_folds': 3, 'eligible_samples': {task: len(task_rows[task]) for task in task_rows}, 'input_sha256': input_hashes, 'code_sha256': code_hash, 'dependencies': {'python': platform.python_version(), 'numpy': np.__version__, 'scikit_learn': version('scikit-learn'), 'scipy': scipy.__version__}, 'warnings': warnings_out, 'output_sha256_scope': 'all output files except this self-referential run.json', 'arms': {arm: {'scorer_status': scorer_results[arm]['status'], 'prediction_rows': sum((row['task_id'] in {TASK_KNOT, TASK_FS} for row in all_predictions[arm])), 'scored_result_dir': f'scores/{arm}'} for arm in ARMS}, 'output_sha256': output_hashes, 'training_history_audit': 'Each scorer run records model-owner metadata. Prediction identity and metrics are verified; model code/weight use is recorded from manifests and hashes but cannot prove absence of unreported data leakage.'}
    return result

def run(args) -> dict:
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise ValueError(f'refusing to overwrite existing output directory: {output_dir}')
    output_dir.mkdir(parents=True)
    started_at = stamp()
    run_record = {'status': 'RUNNING', 'started_at': started_at, 'output_sha256_scope': 'all output files except this self-referential run.json'}
    run_path = output_dir / 'run.json'
    run_path.write_text(json.dumps(run_record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    try:
        result = _execute(args, output_dir, started_at)
    except Exception as error:
        run_record.update({'status': 'FAILED', 'completed_at': stamp(), 'error_type': type(error).__name__, 'error': str(error), 'data_dir': str(args.data_dir.resolve()), 'protocol': str(args.protocol.resolve()), 'feature_manifest': str(args.feature_manifest.resolve())})
        run_path.write_text(json.dumps(run_record, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        raise
    run_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return result

def main() -> int:
    parser = argparse.ArgumentParser(description='Fit the frozen global readers on supplied per-sequence representations.')
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--feature-manifest', type=Path, required=True, help='JSON manifest adjacent to one archive per exact input sequence')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = run(args)
    except Exception as error:
        print(f'Nested evaluation rejected: {error}', file=sys.stderr)
        return 2
    print(json.dumps({'status': result['status'], 'arms': list(result['arms']), 'output_dir': str(args.output_dir)}, ensure_ascii=False))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
