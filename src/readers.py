from __future__ import annotations
import json
import math
import warnings
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
from .evaluate_predictions import TASK_KNOT, TASK_FS, read_table
MAX_ITER = 5000
SEED = 2026
AA20 = 'ACDEFGHIKLMNPQRSTVWY'

def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    for raw in path.read_text(encoding='utf-8-sig').splitlines():
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

def validate_sample_sequence_mappings(key: tuple[str, str], sample: dict, mappings: list[dict]) -> tuple[list[str], list[int]]:
    try:
        sequence_ids = sample['sequence_ids'].split(';')
        lengths = [int(value) for value in sample['lengths'].split(';')]
        expected = Counter(zip(sequence_ids, lengths))
        observed = Counter((row.get('sequence_id'), int(row['length'])) for row in mappings)
    except (KeyError, ValueError) as exc:
        raise ValueError(f'eligible sample has invalid endpoint identity fields: {key}') from exc
    if not sequence_ids or any(not sid for sid in sequence_ids) or len(sequence_ids) != len(lengths) or any(n < 1 for n in lengths):
        raise ValueError(f'eligible sample has invalid endpoint identity fields: {key}')
    if expected != observed:
        raise ValueError(f'eligible sample sequence mapping does not close: {key}')
    if key[0] == TASK_FS and sample.get('unit_type') == 'positive_endpoint_pair':
        roles = Counter(row.get('endpoint_role') for row in mappings)
        if roles != Counter({'endpoint_A': 1, 'endpoint_B': 1}):
            raise ValueError(f'fold-switch pair must map one endpoint_A and one endpoint_B: {key}')
    return sequence_ids, lengths

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
