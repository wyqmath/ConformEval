from __future__ import annotations
import csv, gzip, hashlib
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 2026

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

def validate_label_rows(samples: dict[str, dict], labels_by_sample: dict[str, list[dict]], sequences: dict[str, str], example: bool = False) -> None:
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
        if sample['use_policy'] != 'development_only' or sample['label_view'] != ('synthetic_example' if example else 'reference'):
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
            if example:
                if row['label_source'] != 'synthetic_example':
                    raise ValueError('software example requires synthetic labels')
                continue
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
