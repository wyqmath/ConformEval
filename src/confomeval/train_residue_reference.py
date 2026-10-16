from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import time
import traceback
import warnings
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import scipy
import sklearn
import torch
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from .model_interface import load_residue_representation
from .residue_scoring_core import validate_reference_config
MODEL = 'facebook/esm2_t33_650M_UR50D'
REVISION = '08e4846e537177426273712802403f7ba8261b6c'
LAYERS = (11, 23, 33)
C_GRID = (0.1, 1.0)
AA = 'ACDEFGHIKLMNPQRSTVWY'
SEEDS = (13, 42, 2026)
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 2026
MLP_WIDTH = 64
MLP_BATCH = 8192
MLP_MAX_EPOCHS = 100
MLP_PATIENCE = 10
MLP_LR = 0.001
MLP_DROPOUT = 0.1
DETECTED_WARNINGS: list[dict[str, str]] = []

def stamp() -> str:
    return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode('ascii')).hexdigest()

def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def load_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    with path.open(encoding='utf-8') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                current = line[1:].split()[0]
                if not current or current in records:
                    raise ValueError(f'duplicate/empty FASTA identifier: {current!r}')
                records[current] = []
            elif current is None:
                raise ValueError('FASTA sequence before header')
            else:
                records[current].append(line.upper())
    return {key: ''.join(parts) for key, parts in records.items()}

def make_traditional_features(sequence: str, positions_1based: np.ndarray) -> np.ndarray:
    aa_index = {aa: i for i, aa in enumerate(AA)}
    length = len(sequence)
    out = np.zeros((len(positions_1based), 41), dtype=np.float32)
    for row, pos1 in enumerate(positions_1based):
        pos = int(pos1) - 1
        residue = sequence[pos]
        if residue in aa_index:
            out[row, aa_index[residue]] = 1.0
        out[row, 20] = pos / (length - 1) if length > 1 else 0.0
        window = sequence[max(0, pos - 7):min(length, pos + 8)]
        denom = float(len(window))
        for aa, j in aa_index.items():
            out[row, 21 + j] = window.count(aa) / denom
    return out

def validate_label_rows(samples: dict[str, dict], labels_by_sample: dict[str, list[dict]], sequences: dict[str, str]) -> None:
    if not samples:
        raise ValueError('candidate sample table is empty')
    if set(labels_by_sample) - set(samples):
        raise ValueError('position labels contain samples absent from the frozen sample table')
    for sid, sample in samples.items():
        seq_id = sample['sequence_id']
        seq = sequences.get(seq_id)
        if seq is None:
            raise ValueError(f'sample sequence is absent from FASTA: {sid}')
        seq_sha = sample['sequence_sha256']
        if seq_id != f'SEQ_{seq_sha}' or len(seq) != int(sample['length']) or sha256_text(seq) != seq_sha:
            raise ValueError(f'sequence identity/length/hash mismatch: {sid}')
        if sample['use_policy'] != 'development_only' or sample['label_view'] != 'reference':
            raise ValueError(f'sample is outside the locked development reference task: {sid}')
        rows = labels_by_sample.get(sid, [])
        seen = set()
        counts = Counter()
        for row in rows:
            if row['sequence_id'] != seq_id or row['sequence_sha256'] != seq_sha:
                raise ValueError(f'label sequence identity mismatch: {sid}')
            pos = int(row['position'])
            if pos < 1 or pos > len(seq) or pos in seen:
                raise ValueError(f'duplicate/out-of-range label coordinate: {sid}/{pos}')
            seen.add(pos)
            label = row['label']
            if label not in {'0', '1'}:
                raise ValueError(f'unknown label was included in scoring rows: {sid}/{pos}')
            counts[label] += 1
            if label == '1':
                if not (row['direct_label'] == '1' and row['direct_mask'] == '1' and (row['label_source'] == 'direct_experimental_disorder')):
                    raise ValueError(f'positive lacks direct experimental disorder evidence: {sid}/{pos}')
            elif row['label_source'] == 'direct_experimental_order':
                if not (row['direct_label'] == '0' and row['direct_mask'] == '1'):
                    raise ValueError(f'direct ordered reference label is inconsistent: {sid}/{pos}')
            elif row['label_source'] == 'pdb_residue_observed':
                if row['pdb_observed'] != '1':
                    raise ValueError(f'PDB reference negative is not observed: {sid}/{pos}')
            else:
                raise ValueError(f'unsupported negative label source: {sid}/{pos}')
        expected = {'0': int(sample['negative_positions']), '1': int(sample['positive_positions'])}
        if counts != Counter(expected) or min(counts.get('0', 0), counts.get('1', 0)) == 0:
            raise ValueError(f'sample class counts do not match frozen eligibility: {sid}')
        if counts['0'] + counts['1'] + int(sample['unknown_positions']) != len(seq):
            raise ValueError(f'unknown mask count does not reconcile with sequence length: {sid}')

def validate_fold_roles(samples: dict[str, dict], outer_rows: list[dict], inner_rows: list[dict], expected_outer: int=5, expected_inner: int=3) -> tuple[dict, dict]:
    sample_ids = set(samples)
    groups = {sid: row['binding_group_id'] for sid, row in samples.items()}
    outer: dict[int, dict[str, str]] = defaultdict(dict)
    for row in outer_rows:
        fold, sid = (int(row['outer_fold']), row['sample_id'])
        if sid not in sample_ids or groups[sid] != row['binding_group_id']:
            raise ValueError(f'outer fold identity/group mismatch: {sid}')
        if row['role'] not in {'train', 'test'} or sid in outer[fold]:
            raise ValueError(f'invalid/duplicate outer role: {fold}/{sid}')
        outer[fold][sid] = row['role']
    if set(outer) != set(range(expected_outer)):
        raise ValueError('outer fold ids differ from the frozen five-fold task')
    test_fold_by_group = {}
    for fold, roles in outer.items():
        if set(roles) != sample_ids:
            raise ValueError(f'outer fold {fold} does not cover every sample exactly once')
        role_by_group: dict[str, set[str]] = defaultdict(set)
        for sid, role in roles.items():
            role_by_group[groups[sid]].add(role)
        if any((len(value) != 1 for value in role_by_group.values())):
            raise ValueError(f'binding group crosses outer roles in fold {fold}')
        test_groups = {g for g, value in role_by_group.items() if value == {'test'}}
        if not test_groups:
            raise ValueError(f'outer fold {fold} has no held groups')
        for group in test_groups:
            if group in test_fold_by_group:
                raise ValueError(f'binding group is held out in multiple outer folds: {group}')
            test_fold_by_group[group] = fold
    if set(test_fold_by_group) != set(groups.values()):
        raise ValueError('outer test folds do not hold out every binding group exactly once')
    inner: dict[tuple[int, int], dict[str, str]] = defaultdict(dict)
    for row in inner_rows:
        outer_fold, inner_fold, sid = (int(row['outer_fold']), int(row['inner_fold']), row['sample_id'])
        if outer_fold not in outer or sid not in outer[outer_fold]:
            raise ValueError(f'inner role references unknown outer sample: {outer_fold}/{sid}')
        if outer[outer_fold][sid] != 'train':
            raise ValueError(f'outer test sample entered inner selection: {outer_fold}/{sid}')
        if groups[sid] != row['binding_group_id']:
            raise ValueError(f'inner fold group identity mismatch: {sid}')
        if row['role'] not in {'train', 'validation'} or sid in inner[outer_fold, inner_fold]:
            raise ValueError(f'invalid/duplicate inner role: {outer_fold}/{inner_fold}/{sid}')
        inner[outer_fold, inner_fold][sid] = row['role']
    expected_keys = {(o, i) for o in range(expected_outer) for i in range(expected_inner)}
    if set(inner) != expected_keys:
        raise ValueError('inner fold ids differ from the frozen 5x3 task')
    for outer_fold in range(expected_outer):
        outer_train = {sid for sid, role in outer[outer_fold].items() if role == 'train'}
        valid_fold_by_group = {}
        for inner_fold in range(expected_inner):
            roles = inner[outer_fold, inner_fold]
            if set(roles) != outer_train:
                raise ValueError(f'inner fold {outer_fold}/{inner_fold} does not cover outer train')
            role_by_group: dict[str, set[str]] = defaultdict(set)
            for sid, role in roles.items():
                role_by_group[groups[sid]].add(role)
            if any((len(value) != 1 for value in role_by_group.values())):
                raise ValueError(f'binding group crosses inner roles: {outer_fold}/{inner_fold}')
            val_groups = {g for g, value in role_by_group.items() if value == {'validation'}}
            if not val_groups:
                raise ValueError(f'inner fold {outer_fold}/{inner_fold} has no validation groups')
            for group in val_groups:
                if group in valid_fold_by_group:
                    raise ValueError(f'inner binding group validated more than once: {group}')
                valid_fold_by_group[group] = inner_fold
        if set(valid_fold_by_group) != {groups[sid] for sid in outer_train}:
            raise ValueError(f'inner validations do not cover all outer-training groups: {outer_fold}')
    return (dict(outer), dict(inner))

def protein_ap(scores: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    if len(scores) != len(labels) or len(np.unique(labels)) != 2:
        raise ValueError('per-protein scoring requires paired scores and both classes')
    return (float(average_precision_score(labels, scores)), float(roc_auc_score(labels, scores)))

def make_matrix_by_sequence(task_dir: Path, cache_dir: Path, samples: dict, labels_by_sample: dict, sequences: dict, layers: tuple[int, ...]=LAYERS):
    resid_index_rows = read_tsv(task_dir / 'inputs' / 'resid_index.tsv')
    resid_index = {}
    for row in resid_index_rows:
        name = row['name']
        indices = [int(x) for x in row['resid_indices'].split(',') if x]
        if name in resid_index or not indices or len(set(indices)) != len(indices):
            raise ValueError(f'duplicate/empty/duplicate-coordinate residue index: {name}')
        if indices != sorted(indices):
            raise ValueError(f'residue index is not strictly ordered: {name}')
        resid_index[name] = indices
    wanted_positions: dict[str, set[int]] = defaultdict(set)
    for sid, rows in labels_by_sample.items():
        seq_id = samples[sid]['sequence_id']
        wanted_positions[seq_id].update((int(row['position']) for row in rows))
    if {k: set(v) for k, v in resid_index.items()} != dict(wanted_positions):
        raise ValueError('resid_index coordinate union differs from valid selected labels')
    manifest_path = cache_dir / 'extract_manifest.tsv'
    manifest_rows = read_tsv(manifest_path)
    manifest = {row['name']: row for row in manifest_rows}
    if len(manifest) != len(manifest_rows) or set(manifest) != set(resid_index):
        raise ValueError('representation manifest sequence set differs from locked residue index')
    n_coords = sum(map(len, resid_index.values()))
    x_layers = [np.empty((n_coords, 1280), dtype=np.float16) for _ in layers]
    traditional = np.empty((n_coords, 41), dtype=np.float32)
    coord_rows = {}
    seq_offsets = {}
    cursor = 0
    archive_hashes = {}
    for seq_id in sorted(resid_index):
        positions = resid_index[seq_id]
        sequence = sequences[seq_id]
        if len(sequence) != int(manifest[seq_id]['len_full']):
            raise ValueError(f'cache manifest length mismatch: {seq_id}')
        mrow = manifest[seq_id]
        if mrow['revision'] != REVISION or mrow['seq_sha256'] != seq_id.removeprefix('SEQ_') or mrow['truncated'] not in {'False', 'false', '0'} or (int(mrow['len_used']) != len(sequence)):
            raise ValueError(f'cache manifest identity/revision/truncation mismatch: {seq_id}')
        archive = cache_dir / f"{mrow['key']}.npz"
        if not archive.is_file():
            raise ValueError(f'representation archive missing: {seq_id}')
        archive_sha = sha256_file(archive)
        archive_hashes[archive.name] = archive_sha
        values, indices, dim = load_residue_representation(archive, seq_id.removeprefix('SEQ_'), len(sequence), archive_sha, layers, MODEL, REVISION)
        if values.shape != (len(layers), len(positions), 1280) or dim != 1280:
            raise ValueError(f'cache matrix shape differs for {seq_id}: {values.shape}')
        if indices.tolist() != positions or int(mrow['n_resid_stored']) != len(positions):
            raise ValueError(f'cached residue coordinates differ from locked selection: {seq_id}')
        with np.load(archive, allow_pickle=False) as z:
            meta = json.loads(str(z['meta'].item()))
        if meta.get('name') != seq_id or meta.get('device') != 'cuda':
            raise ValueError(f'cache metadata name/device mismatch: {seq_id}')
        if not np.isfinite(values).all():
            raise ValueError(f'cache contains non-finite vectors: {seq_id}')
        seq_rows = np.arange(cursor, cursor + len(positions), dtype=np.int64)
        coord_rows[seq_id] = dict(zip(positions, seq_rows.tolist()))
        seq_offsets[seq_id] = seq_rows
        values = np.moveaxis(values, 0, 1)
        for li, x in enumerate(layers):
            x_layers[li][seq_rows] = values[:, li, :].astype(np.float16, copy=False)
        traditional[seq_rows] = make_traditional_features(sequence, np.asarray(positions))
        cursor += len(positions)
    sample_indices, sample_y = ({}, {})
    for sid, rows in labels_by_sample.items():
        seq_id = samples[sid]['sequence_id']
        sorted_rows = sorted(rows, key=lambda row: int(row['position']))
        sample_indices[sid] = np.asarray([coord_rows[seq_id][int(row['position'])] for row in sorted_rows], dtype=np.int64)
        sample_y[sid] = np.asarray([int(row['label']) for row in sorted_rows], dtype=np.int8)
    return (x_layers, traditional, sample_indices, sample_y, archive_hashes, manifest_path)

def rows_for_samples(sample_ids: list[str], sample_indices: dict, sample_y: dict):
    if not sample_ids:
        raise ValueError('empty protein set for model fit or evaluation')
    indices = np.concatenate([sample_indices[sid] for sid in sample_ids])
    targets = np.concatenate([sample_y[sid] for sid in sample_ids])
    if len(indices) != len(targets) or len(np.unique(targets)) != 2:
        raise ValueError('residue training/validation domain lacks paired rows or both classes')
    return (indices, targets)

def fit_linear(matrix: np.ndarray, train_ids: list[str], C: float, sample_indices: dict, sample_y: dict):
    idx, y = rows_for_samples(train_ids, sample_indices, sample_y)
    X = np.asarray(matrix[idx], dtype=np.float32)
    if not np.isfinite(X).all():
        raise ValueError('linear training matrix contains non-finite values')
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    model = LogisticRegression(C=C, class_weight='balanced', solver='liblinear', max_iter=5000, random_state=2026)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        model.fit(Xs, y)
    for item in caught:
        DETECTED_WARNINGS.append({'type': item.category.__name__, 'message': str(item.message)})
    if any((issubclass(item.category, ConvergenceWarning) for item in caught)):
        raise ValueError('linear residue reader did not converge')
    if list(model.classes_) != [0, 1]:
        raise ValueError('linear reader class order differs from [0,1]')
    return (scaler, model, len(y))

def predict_linear(matrix: np.ndarray, protein_ids: list[str], scaler, model, sample_indices: dict) -> dict[str, np.ndarray]:
    out = {}
    for sid in protein_ids:
        idx = sample_indices[sid]
        X = np.asarray(matrix[idx], dtype=np.float32)
        score = model.predict_proba(scaler.transform(X))[:, 1]
        if len(score) != len(idx) or not np.isfinite(score).all():
            raise ValueError(f'invalid residue-level linear scores: {sid}')
        out[sid] = score.astype(np.float32)
    return out

def choose_linear_C(arm: str, matrix: np.ndarray, outer_fold: int, train_ids: list[str], inner_by_fold: dict, samples: dict, sample_indices: dict, sample_y: dict, selection_rows: list[dict]) -> float:
    metrics_by_C = {}
    for C in C_GRID:
        protein_aps = []
        for inner_fold in range(3):
            roles = inner_by_fold[outer_fold, inner_fold]
            inner_train = [sid for sid in train_ids if roles[sid] == 'train']
            validation = [sid for sid in train_ids if roles[sid] == 'validation']
            scaler, model, ntrain = fit_linear(matrix, inner_train, C, sample_indices, sample_y)
            scores = predict_linear(matrix, validation, scaler, model, sample_indices)
            ap_vals = [protein_ap(scores[sid], sample_y[sid])[0] for sid in validation]
            if not ap_vals:
                raise ValueError(f'no scored proteins in inner validation {outer_fold}/{inner_fold}')
            protein_aps.extend(ap_vals)
            selection_rows.append({'task_id': 'DISORDER-STRUCTURE-REFERENCE-CANDIDATE', 'outer_fold': outer_fold, 'inner_fold': inner_fold, 'arm': arm, 'C': C, 'n_train_residues': ntrain, 'n_validation_proteins': len(validation), 'mean_validation_protein_AP': float(np.mean(ap_vals)), 'phase': 'inner'})
        value = float(np.mean(protein_aps))
        metrics_by_C[C] = value
        selection_rows.append({'task_id': 'DISORDER-STRUCTURE-REFERENCE-CANDIDATE', 'outer_fold': outer_fold, 'inner_fold': 'all_mean', 'arm': arm, 'C': C, 'n_validation_proteins': len(protein_aps), 'mean_validation_protein_AP': value, 'phase': 'inner_mean'})
    return min((C for C, metric in metrics_by_C.items() if np.isclose(metric, max(metrics_by_C.values()), rtol=0, atol=1e-12)))

class Width64Reader(torch.nn.Module):

    def __init__(self, input_dim: int):
        super().__init__()
        self.net = torch.nn.Sequential(torch.nn.Linear(input_dim, MLP_WIDTH), torch.nn.ReLU(), torch.nn.Dropout(MLP_DROPOUT), torch.nn.Linear(MLP_WIDTH, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)

def fit_mlp_epochs(X: np.ndarray, y: np.ndarray, epochs: int, seed: int, device: torch.device):
    if epochs < 1 or len(np.unique(y)) != 2:
        raise ValueError('invalid MLP fit epoch count or class coverage')
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = Width64Reader(X.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=MLP_LR)
    n0, n1 = np.bincount(y.astype(np.int64), minlength=2)
    weights = torch.tensor(np.where(y == 1, len(y) / (2.0 * n1), len(y) / (2.0 * n0)), dtype=torch.float32, device=device)
    xt = torch.as_tensor(X, dtype=torch.float32, device=device)
    yt = torch.as_tensor(y.astype(np.float32), dtype=torch.float32, device=device)
    generator = torch.Generator().manual_seed(seed)
    loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(xt, yt, weights), batch_size=MLP_BATCH, shuffle=True, generator=generator, drop_last=False)
    loss_fn = torch.nn.BCEWithLogitsLoss(reduction='none')
    last_loss = None
    for epoch in range(epochs):
        model.train()
        epoch_losses = []
        for xb, yb, wb in loader:
            optimizer.zero_grad(set_to_none=True)
            loss = (loss_fn(model(xb), yb) * wb).mean()
            if not torch.isfinite(loss):
                raise ValueError(f'non-finite MLP loss at epoch {epoch + 1}')
            loss.backward()
            optimizer.step()
            epoch_losses.append(float(loss.detach().cpu()))
        last_loss = float(np.mean(epoch_losses))
    return (model, last_loss)

def predict_mlp(model, scaler: StandardScaler, matrix: np.ndarray, protein_ids: list[str], sample_indices: dict, device: torch.device) -> dict[str, np.ndarray]:
    model.eval()
    out = {}
    with torch.no_grad():
        for sid in protein_ids:
            idx = sample_indices[sid]
            X = scaler.transform(np.asarray(matrix[idx], dtype=np.float32))
            xt = torch.as_tensor(X, dtype=torch.float32, device=device)
            scores = torch.sigmoid(model(xt)).detach().cpu().numpy()
            if len(scores) != len(idx) or not np.isfinite(scores).all():
                raise ValueError(f'invalid MLP residue scores: {sid}')
            out[sid] = scores.astype(np.float32)
    return out

def early_stop_mlp(matrix: np.ndarray, train_ids: list[str], validation_ids: list[str], sample_indices: dict, sample_y: dict, seed: int, device: torch.device) -> tuple[int, float, int]:
    tr_idx, ytr = rows_for_samples(train_ids, sample_indices, sample_y)
    va_idx, yva = rows_for_samples(validation_ids, sample_indices, sample_y)
    scaler = StandardScaler()
    Xtr = scaler.fit_transform(np.asarray(matrix[tr_idx], dtype=np.float32))
    Xva = scaler.transform(np.asarray(matrix[va_idx], dtype=np.float32))
    if not np.isfinite(Xtr).all() or not np.isfinite(Xva).all():
        raise ValueError('MLP standardized inputs contain non-finite values')
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = Width64Reader(Xtr.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=MLP_LR)
    n0, n1 = np.bincount(ytr.astype(np.int64), minlength=2)
    sample_weight = np.where(ytr == 1, len(ytr) / (2.0 * n1), len(ytr) / (2.0 * n0)).astype(np.float32)
    xt = torch.as_tensor(Xtr, dtype=torch.float32, device=device)
    yt = torch.as_tensor(ytr.astype(np.float32), dtype=torch.float32, device=device)
    wt = torch.as_tensor(sample_weight, dtype=torch.float32, device=device)
    generator = torch.Generator().manual_seed(seed)
    loader = torch.utils.data.DataLoader(torch.utils.data.TensorDataset(xt, yt, wt), batch_size=MLP_BATCH, shuffle=True, generator=generator, drop_last=False)
    xval = torch.as_tensor(Xva, dtype=torch.float32, device=device)
    loss_fn = torch.nn.BCEWithLogitsLoss(reduction='none')
    best_ap, best_epoch, wait = (-math.inf, 0, 0)
    best_group_scores = None
    for epoch in range(1, MLP_MAX_EPOCHS + 1):
        model.train()
        for xb, yb, wb in loader:
            optimizer.zero_grad(set_to_none=True)
            loss = (loss_fn(model(xb), yb) * wb).mean()
            if not torch.isfinite(loss):
                raise ValueError(f'non-finite MLP loss at epoch {epoch}')
            loss.backward()
            optimizer.step()
        model.eval()
        with torch.no_grad():
            val_score = torch.sigmoid(model(xval)).detach().cpu().numpy()
        offset, per_protein = (0, [])
        for sid in validation_ids:
            n = len(sample_indices[sid])
            ap, _ = protein_ap(val_score[offset:offset + n], sample_y[sid])
            per_protein.append(ap)
            offset += n
        value = float(np.mean(per_protein))
        if value > best_ap + 1e-12:
            best_ap, best_epoch, wait = (value, epoch, 0)
        else:
            wait += 1
            if wait >= MLP_PATIENCE:
                break
    if best_epoch == 0:
        raise ValueError('MLP early stopping did not find a valid epoch')
    return (best_epoch, best_ap, len(ytr))

def group_bootstrap(metrics_rows: list[dict], arms: list[str]) -> dict:
    proteins = sorted({row['sample_id'] for row in metrics_rows})
    groups = sorted({row['binding_group_id'] for row in metrics_rows})
    group_proteins = {group: [sid for sid in proteins if next((r['binding_group_id'] for r in metrics_rows if r['sample_id'] == sid)) == group] for group in groups}
    by_arm = {arm: {row['sample_id']: float(row['average_precision']) for row in metrics_rows if row['arm'] == arm} for arm in arms}
    if any((set(values) != set(proteins) for values in by_arm.values())):
        raise ValueError('arm protein-metric coverage differs')
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(0, len(groups), size=(BOOTSTRAPS, len(groups)))
    out = {}
    reference = 'ESM2-L33'
    for arm in arms:
        means, deltas = ([], [])
        for draw in draws:
            sampled = [sid for j in draw for sid in group_proteins[groups[int(j)]]]
            means.append(float(np.mean([by_arm[arm][sid] for sid in sampled])))
            if arm != reference:
                deltas.append(float(np.mean([by_arm[arm][sid] - by_arm[reference][sid] for sid in sampled])))
        record = {'macro_ap_ci95': np.quantile(means, [0.025, 0.975]).tolist(), 'replicates': BOOTSTRAPS, 'bootstrap_unit': 'binding_group_id'}
        if arm != reference:
            record['delta_vs_ESM2_L33'] = float(np.mean([by_arm[arm][sid] - by_arm[reference][sid] for sid in proteins]))
            record['delta_vs_ESM2_L33_ci95'] = np.quantile(deltas, [0.025, 0.975]).tolist()
        out[arm] = record
    return out

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--task-dir', type=Path, required=True)
    ap.add_argument('--config-dir', type=Path, required=True, help='Directory containing the residue task protocols and extraction request')
    ap.add_argument('--cache-dir', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--device', choices=['cuda', 'cpu'], default='cuda')
    args = ap.parse_args()
    task_dir, cache_dir, out_dir = (args.task_dir.resolve(), args.cache_dir.resolve(), args.output_dir.resolve())
    if out_dir.exists():
        raise ValueError(f'output directory already exists: {out_dir}')
    if not task_dir.is_dir() or not cache_dir.is_dir():
        raise ValueError('task or representation cache directory is missing')
    from .evaluate_residue_predictions import load_task
    load_task(task_dir)
    config_dir = args.config_dir.resolve()
    validate_reference_config(config_dir)
    protocol_path = config_dir / 'evaluation_protocol.json'
    protocol = json.loads(protocol_path.read_text(encoding='utf-8'))
    if protocol.get('linear_arms', {}).get('revision') != REVISION:
        raise ValueError('evaluation protocol encoder revision mismatch')
    device = torch.device('cuda' if args.device == 'cuda' and torch.cuda.is_available() else 'cpu')
    started_at, start_clock = (stamp(), time.time())
    out_dir.mkdir(parents=True, exist_ok=False)
    (out_dir / 'parameters').mkdir()
    run = {'status': 'RUNNING', 'started_at': started_at, 'device': str(device), 'encoder': MODEL, 'revision': REVISION, 'layers': list(LAYERS), 'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'scikit_learn': sklearn.__version__, 'torch': torch.__version__, 'cuda': torch.version.cuda}, 'warnings': [], 'input_sha256': {}, 'cache_archive_sha256': {}, 'code_sha256_at_start': {}, 'output_sha256': {}}
    (out_dir / 'run.json').write_text(json.dumps(run, indent=2, sort_keys=True) + '\n')
    try:
        input_files = [protocol_path, config_dir / 'data_protocol.json', task_dir / 'run_manifest.json', config_dir / 'extraction_request.json', task_dir / 'selected_samples.tsv', task_dir / 'selected_labels.tsv.gz', task_dir / 'outer_fold_roles.tsv', task_dir / 'inner_fold_roles.tsv', task_dir / 'fold_qc.tsv', task_dir / 'inputs' / 'resid_index.tsv', task_dir / 'inputs' / 'sequences.fasta', cache_dir / 'extract_manifest.tsv']
        for path in input_files:
            if not path.is_file():
                raise ValueError(f'required input file missing: {path}')
            run['input_sha256'][str(path)] = sha256_file(path)
        code_files = [Path(__file__).resolve(), Path(__file__).resolve().parent / 'model_interface.py', Path(__file__).resolve().parent / 'residue_scoring_core.py']
        run['code_sha256_at_start'] = {str(path): sha256_file(path) for path in code_files}
        samples_rows = read_tsv(task_dir / 'selected_samples.tsv')
        samples = {row['sample_id']: row for row in samples_rows}
        if len(samples) != len(samples_rows) or len(samples) != 1279:
            raise ValueError('candidate sample denominator differs from frozen 1279 proteins')
        label_rows = read_tsv(task_dir / 'selected_labels.tsv.gz')
        labels_by_sample: dict[str, list[dict]] = defaultdict(list)
        for row in label_rows:
            labels_by_sample[row['sample_id']].append(row)
        sequences = load_fasta(task_dir / 'inputs' / 'sequences.fasta')
        validate_label_rows(samples, labels_by_sample, sequences)
        if len(sequences) != 1275 or sum(map(len, labels_by_sample.values())) != 391718:
            raise ValueError('candidate sequence or valid-residue denominator differs from frozen inputs')
        outer, inner = validate_fold_roles(samples, read_tsv(task_dir / 'outer_fold_roles.tsv'), read_tsv(task_dir / 'inner_fold_roles.tsv'))
        print('PRECHECK labels, sequence identities, masks, and 5x3 grouped folds: PASS', flush=True)
        x_layers, x_trad, sample_indices, sample_y, cache_hashes, manifest_path = make_matrix_by_sequence(task_dir, cache_dir, samples, labels_by_sample, sequences)
        run['cache_archive_sha256'] = cache_hashes
        run['n_cache_archives'] = len(cache_hashes)
        if len(cache_hashes) != 1275:
            raise ValueError('representation archive count differs from 1275 unique sequences')
        run['n_unique_labeled_coordinates'] = int(x_layers[0].shape[0])
        print(f'PRECHECK representations: {len(cache_hashes)} archives, {x_layers[0].shape[0]} unique labeled coordinates; PASS', flush=True)
        groups_by_sample = {sid: row['binding_group_id'] for sid, row in samples.items()}
        arms = [*(f'ESM2-L{layer}' for layer in LAYERS), 'traditional-AA-context', *(f'MLP-h33-seed{seed}' for seed in SEEDS)]
        predictions: dict[str, dict[str, np.ndarray]] = {arm: {} for arm in arms}
        selection_rows, fit_rows = ([], [])
        matrices = {**{f'ESM2-L{layer}': x_layers[i] for i, layer in enumerate(LAYERS)}, 'traditional-AA-context': x_trad}
        for outer_fold in range(5):
            train_ids = sorted((sid for sid, role in outer[outer_fold].items() if role == 'train'))
            test_ids = sorted((sid for sid, role in outer[outer_fold].items() if role == 'test'))
            for arm, matrix in matrices.items():
                C = choose_linear_C(arm, matrix, outer_fold, train_ids, inner, samples, sample_indices, sample_y, selection_rows)
                scaler, model, ntrain = fit_linear(matrix, train_ids, C, sample_indices, sample_y)
                scores = predict_linear(matrix, test_ids, scaler, model, sample_indices)
                predictions[arm].update(scores)
                fit_rows.append({'outer_fold': outer_fold, 'arm': arm, 'selected_C': C, 'n_train_proteins': len(train_ids), 'n_test_proteins': len(test_ids), 'n_train_residues': ntrain})
                param_name = f'linear_{arm}_outer{outer_fold:02d}.npz'
                np.savez_compressed(out_dir / 'parameters' / param_name, coef=model.coef_, intercept=model.intercept_, scaler_mean=scaler.mean_, scaler_scale=scaler.scale_, selected_C=np.asarray(C), classes=model.classes_, arm=np.asarray(arm), outer_fold=np.asarray(outer_fold))
                print(f'LINEAR outer={outer_fold} arm={arm} C={C:g} proteins={len(train_ids)}/{len(test_ids)}', flush=True)
        if device.type != 'cuda':
            print('WARNING MLP capacity diagnostic running on CPU', flush=True)
            run['warnings'].append('CUDA unavailable; MLP diagnostic ran on CPU')
        h33 = matrices['ESM2-L33']
        for outer_fold in range(5):
            train_ids = sorted((sid for sid, role in outer[outer_fold].items() if role == 'train'))
            test_ids = sorted((sid for sid, role in outer[outer_fold].items() if role == 'test'))
            for seed in SEEDS:
                best_epochs, inner_records = ([], [])
                for inner_fold in range(3):
                    roles = inner[outer_fold, inner_fold]
                    inner_train = [sid for sid in train_ids if roles[sid] == 'train']
                    validation = [sid for sid in train_ids if roles[sid] == 'validation']
                    epoch, metric, ntrain = early_stop_mlp(h33, inner_train, validation, sample_indices, sample_y, seed, device)
                    best_epochs.append(epoch)
                    inner_records.append({'outer_fold': outer_fold, 'inner_fold': inner_fold, 'arm': f'MLP-h33-seed{seed}', 'seed': seed, 'best_epoch': epoch, 'best_validation_mean_protein_AP': metric, 'n_train_residues': ntrain, 'n_validation_proteins': len(validation)})
                selection_rows.extend(inner_records)
                final_epochs = int(np.median(best_epochs))
                tr_idx, ytr = rows_for_samples(train_ids, sample_indices, sample_y)
                scaler = StandardScaler()
                Xtr = scaler.fit_transform(np.asarray(h33[tr_idx], dtype=np.float32))
                model, final_loss = fit_mlp_epochs(Xtr, ytr, final_epochs, seed, device)
                scores = predict_mlp(model, scaler, h33, test_ids, sample_indices, device)
                arm = f'MLP-h33-seed{seed}'
                predictions[arm].update(scores)
                fit_rows.append({'outer_fold': outer_fold, 'arm': arm, 'seed': seed, 'selected_epochs': final_epochs, 'inner_best_epochs': ','.join(map(str, best_epochs)), 'n_train_proteins': len(train_ids), 'n_test_proteins': len(test_ids), 'n_train_residues': len(ytr), 'final_train_loss': final_loss})
                state = {f"state_{k.replace('.', '_')}": v.detach().cpu().numpy() for k, v in model.state_dict().items()}
                state.update({'scaler_mean': scaler.mean_, 'scaler_scale': scaler.scale_, 'seed': np.asarray(seed), 'final_epochs': np.asarray(final_epochs), 'inner_best_epochs': np.asarray(best_epochs, dtype=np.int32), 'revision': np.asarray(REVISION)})
                np.savez_compressed(out_dir / 'parameters' / f'mlp_h33_seed{seed}_outer{outer_fold:02d}.npz', **state)
                print(f"MLP outer={outer_fold} seed={seed} epochs={final_epochs} AP-inner={np.mean([r['best_validation_mean_protein_AP'] for r in inner_records]):.5f}", flush=True)
        for arm, scores in predictions.items():
            if set(scores) != set(samples):
                raise ValueError(f'outer test predictions do not cover every protein exactly once: {arm}')
        protein_rows = []
        for arm in arms:
            for sid in sorted(samples):
                ap_value, auc_value = protein_ap(predictions[arm][sid], sample_y[sid])
                fold_id = next((f for f, roles in outer.items() if roles[sid] == 'test'))
                protein_rows.append({'task_id': 'DISORDER-STRUCTURE-REFERENCE-CANDIDATE', 'sample_id': sid, 'binding_group_id': groups_by_sample[sid], 'outer_fold': fold_id, 'arm': arm, 'n_labeled_residues': len(sample_y[sid]), 'n_positive': int(sample_y[sid].sum()), 'n_negative': int(len(sample_y[sid]) - sample_y[sid].sum()), 'average_precision': ap_value, 'auroc': auc_value})
        write_tsv(out_dir / 'protein_metrics.tsv', list(protein_rows[0]), protein_rows)
        pred_fields = ['task_id', 'sample_id', 'binding_group_id', 'outer_fold', 'arm', 'position_1based', 'label', 'mask', 'score']
        pred_rows = []
        for arm in arms:
            for sid in sorted(samples):
                ordered = sorted(labels_by_sample[sid], key=lambda row: int(row['position']))
                score_values = predictions[arm][sid]
                for row, score in zip(ordered, score_values):
                    pred_rows.append({'task_id': 'DISORDER-STRUCTURE-REFERENCE-CANDIDATE', 'sample_id': sid, 'binding_group_id': groups_by_sample[sid], 'outer_fold': next((f for f, roles in outer.items() if roles[sid] == 'test')), 'arm': arm, 'position_1based': row['position'], 'label': row['label'], 'mask': 1, 'score': float(score)})
        write_tsv(out_dir / 'residue_predictions.tsv.gz', pred_fields, pred_rows)
        selection_fields = sorted({key for row in selection_rows for key in row})
        write_tsv(out_dir / 'inner_selection.tsv', selection_fields, selection_rows)
        write_tsv(out_dir / 'outer_fits.tsv', sorted({k for row in fit_rows for k in row}), fit_rows)
        summaries = {}
        for arm in arms:
            rows = [r for r in protein_rows if r['arm'] == arm]
            all_scores = np.concatenate([predictions[arm][sid] for sid in sorted(samples)])
            all_labels = np.concatenate([sample_y[sid] for sid in sorted(samples)])
            summaries[arm] = {'n_proteins': len(rows), 'mean_per_protein_AP': float(np.mean([r['average_precision'] for r in rows])), 'sd_per_protein_AP': float(np.std([r['average_precision'] for r in rows])), 'mean_per_protein_AUROC': float(np.mean([r['auroc'] for r in rows])), 'pooled_residue_AP_secondary': float(average_precision_score(all_labels, all_scores))}
        summaries['group_bootstrap'] = group_bootstrap(protein_rows, arms)
        metrics_path = out_dir / 'metrics.json'
        metrics_path.write_text(json.dumps(summaries, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        report = ['# DISORDER-STRUCTURE-REFERENCE-CANDIDATE nested residue evaluation', '', f'Completed: {stamp()} (Asia/Shanghai).', '', 'This is a separately versioned development evaluation. Label 1 is direct DisProt disorder; label 0 combines direct DisProt order and amino-acid-matched PDB-observed reference positions. It is not strict direct ordered/disordered truth and does not establish family-independent generalization.', '', f'The frozen task contains {len(samples)} proteins in {len(set(groups_by_sample.values()))} binding groups. It scores 391,718 labeled protein-residue rows over 1,279 proteins; unknown positions remain excluded. All outer/inner roles are consumed from the saved grouped fold tables.', '', '| Arm | Mean protein AP (primary) | Mean protein AUROC | Pooled residue AP (secondary) |', '|---|---:|---:|---:|']
        for arm in arms:
            m = summaries[arm]
            report.append(f"| {arm} | {m['mean_per_protein_AP']:.6f} | {m['mean_per_protein_AUROC']:.6f} | {m['pooled_residue_AP_secondary']:.6f} |")
        report.extend(['', '95% intervals are binding-group bootstrap intervals conditional on the fixed outer predictions. MLP seeds are all reported; no seed is selected by outer performance.', '', 'The traditional baseline uses amino-acid one-hot, relative position, and centered 15-residue composition. All scaling is fit on training residues only. The MLP capacity diagnostic uses hidden layer 33 and inner-fold early stopping.', ''])
        (out_dir / 'REPORT.md').write_text('\n'.join(report), encoding='utf-8')
        run['status'] = 'COMPLETE_PENDING_INDEPENDENT_REVIEW'
        run['ended_at'] = stamp()
        run['duration_seconds'] = time.time() - start_clock
        run['arms'] = arms
        run['n_predictions'] = len(pred_rows)
        run['warnings'].extend(DETECTED_WARNINGS)
        run['code_sha256_at_end'] = {str(path): sha256_file(path) for path in code_files}
        if run['code_sha256_at_start'] != run['code_sha256_at_end']:
            raise ValueError('runner or model-interface code changed during the run')
        for path in sorted((p for p in out_dir.rglob('*') if p.is_file() and p.name != 'run.json')):
            run['output_sha256'][str(path.relative_to(out_dir))] = sha256_file(path)
        run['results_summary'] = summaries
    except Exception:
        run['status'] = 'FAILED'
        run['ended_at'] = stamp()
        run['duration_seconds'] = time.time() - start_clock
        run['error'] = traceback.format_exc()
        run['warnings'].extend(DETECTED_WARNINGS)
        (out_dir / 'run.json').write_text(json.dumps(run, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        raise
    (out_dir / 'run.json').write_text(json.dumps(run, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
