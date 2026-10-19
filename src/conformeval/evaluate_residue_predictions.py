from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
from collections import defaultdict
from pathlib import Path
import numpy as np
import scipy
import sklearn
from .residue_scoring_core import read_tsv, load_fasta, validate_label_rows, validate_fold_roles, protein_ap, group_bootstrap, write_tsv, stamp
TASK = 'DISORDER-STRUCTURE-REFERENCE-CANDIDATE'
REFERENCE_ARMS = {'ESM2-L11', 'ESM2-L23', 'ESM2-L33', 'traditional-AA-context', 'MLP-h33-seed13', 'MLP-h33-seed42', 'MLP-h33-seed2026'}

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

def check_manifest(root):
    manifest = json.loads((root / 'manifest.json').read_text(encoding='utf-8-sig'))
    hashes = {'manifest.json': sha256(root / 'manifest.json')}
    for name, record in manifest['files'].items():
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError('unsafe manifest path')
        p = root / rel
        if not p.is_file() or sha256(p) != record['sha256'] or p.stat().st_size != record['bytes']:
            raise ValueError(f'dataset hash/size mismatch: {name}')
        hashes[name] = record['sha256']
    return hashes

def load_task(root, example=False):
    hashes = check_manifest(root)
    manifest = json.loads((root / 'manifest.json').read_text(encoding='utf-8-sig'))
    if example and (manifest.get('artifact_id') != 'confomeval-software-example-residue-v1' or manifest.get('scope') != 'software_example'):
        raise ValueError('example mode requires the bundled software example')
    if not example and hashes['manifest.json'] != 'd8527744b6c8b8997f39dca56e49725a815207237ed4f97547d3ae0c2df7907a':
        raise ValueError('unsupported disorder reference data version')
    rows = read_tsv(root / 'selected_samples.tsv')
    samples = {r['sample_id']: r for r in rows}
    if len(samples) != len(rows):
        raise ValueError('duplicate sample identity')
    labels = defaultdict(list)
    for r in read_tsv(root / 'selected_labels.tsv.gz'):
        labels[r['sample_id']].append(r)
    validate_label_rows(samples, labels, load_fasta(root / 'inputs/sequences.fasta'), example=example)
    sequence_groups = defaultdict(set)
    for s in samples.values():
        sequence_groups[s['sequence_sha256']].add(s['binding_group_id'])
    if any((len(gs) != 1 for gs in sequence_groups.values())):
        raise ValueError('exact sequence crosses binding groups')
    outer, inner = validate_fold_roles(samples, read_tsv(root / 'outer_fold_roles.tsv'), read_tsv(root / 'inner_fold_roles.tsv'))
    expected = {sid: {int(r['position']): int(r['label']) for r in rs} for sid, rs in labels.items()}
    folds = {sid: f for f, roles in outer.items() for sid, role in roles.items() if role == 'test'}
    return (samples, expected, folds, hashes)

def read_scores(path, samples, expected, folds, arms, default_arm):
    scores = {arm: {sid: {} for sid in samples} for arm in arms}
    opener = gzip.open if path.suffix == '.gz' else open
    required = {'task_id', 'sample_id', 'position_1based', 'outer_fold', 'score'}
    with opener(path, 'rt', newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('missing prediction columns')
        for r in reader:
            sid = r['sample_id']
            arm = r.get('arm', default_arm)
            if r['task_id'] != TASK or sid not in samples or arm not in arms:
                raise ValueError('unexpected task/sample/arm')
            pos = int(r['position_1based'])
            if pos not in expected[sid]:
                raise ValueError('prediction coordinate is unknown, masked, or out of range')
            if int(r['outer_fold']) != folds[sid]:
                raise ValueError('prediction is outside its frozen outer test fold')
            if 'binding_group_id' in r and r['binding_group_id'] != samples[sid]['binding_group_id']:
                raise ValueError('prediction group mismatch')
            if 'mask' in r and r['mask'] != '1':
                raise ValueError('masked prediction included')
            if 'label' in r and r['label'] != str(expected[sid][pos]):
                raise ValueError('prediction label mismatch')
            score = float(r['score'])
            if not math.isfinite(score):
                raise ValueError('non-finite prediction')
            if pos in scores[arm][sid]:
                raise ValueError('duplicate prediction identity')
            scores[arm][sid][pos] = score
    for arm in arms:
        for sid, values in scores[arm].items():
            if set(values) != set(expected[sid]):
                raise ValueError(f'incomplete prediction coordinates: {arm}/{sid}')
    return scores

def calculate(samples, expected, folds, scores):
    rows = []
    summaries = {}
    arms = sorted(scores)
    for arm in arms:
        pooled_scores = []
        pooled_labels = []
        for sid in sorted(samples):
            pos = sorted(expected[sid])
            y = np.array([expected[sid][p] for p in pos], dtype=np.int8)
            s = np.array([scores[arm][sid][p] for p in pos], dtype=np.float64)
            ap, auc = protein_ap(s, y)
            rows.append({'task_id': TASK, 'sample_id': sid, 'binding_group_id': samples[sid]['binding_group_id'], 'outer_fold': folds[sid], 'arm': arm, 'n_labeled_residues': len(y), 'n_positive': int(y.sum()), 'n_negative': int(len(y) - y.sum()), 'average_precision': ap, 'auroc': auc})
            pooled_scores.append(s)
            pooled_labels.append(y)
        rs = [r for r in rows if r['arm'] == arm]
        summaries[arm] = {'n_proteins': len(rs), 'mean_per_protein_AP': float(np.mean([r['average_precision'] for r in rs])), 'sd_per_protein_AP': float(np.std([r['average_precision'] for r in rs])), 'mean_per_protein_AUROC': float(np.mean([r['auroc'] for r in rs])), 'pooled_residue_AP_secondary': float(sklearn.metrics.average_precision_score(np.concatenate(pooled_labels), np.concatenate(pooled_scores)))}
    if set(arms) == REFERENCE_ARMS:
        summaries['group_bootstrap'] = group_bootstrap(rows, arms)
    else:
        groups = sorted({r['binding_group_id'] for r in rows})
        rng = np.random.default_rng(2026)
        draws = rng.integers(0, len(groups), size=(2000, len(groups)))
        boot = {}
        for arm in arms:
            by_group = {g: [r['average_precision'] for r in rows if r['arm'] == arm and r['binding_group_id'] == g] for g in groups}
            totals = np.array([sum(by_group[g]) for g in groups])
            counts = np.array([len(by_group[g]) for g in groups])
            means = totals[draws].sum(axis=1) / counts[draws].sum(axis=1)
            boot[arm] = {'macro_ap_ci95': np.quantile(means, [0.025, 0.975]).tolist(), 'replicates': 2000, 'bootstrap_unit': 'binding_group_id'}
        summaries['group_bootstrap'] = boot
    return (rows, summaries)

def main():
    ap = argparse.ArgumentParser(description='Score full out-of-fold predictions for the frozen disorder reference task.')
    ap.add_argument('--task-dir', required=True, type=Path)
    ap.add_argument('--predictions', required=True, type=Path)
    ap.add_argument('--model-metadata', type=Path)
    ap.add_argument('--example', action='store_true', help='Score the bundled synthetic software example')
    ap.add_argument('--reference', action='store_true')
    ap.add_argument('--output-dir', required=True, type=Path)
    a = ap.parse_args()
    out = a.output_dir.resolve()
    task = a.task_dir.resolve()
    pred = a.predictions.resolve()
    if out.exists():
        raise ValueError('output directory already exists')
    if a.reference:
        bundle = task.parent.parent / 'results/disorder_reference_esm2'
        check_manifest(bundle)
        if pred != bundle / 'residue_predictions.tsv.gz':
            raise ValueError('reference mode requires the archived reference predictions')
        metadata = {'model_id': 'facebook/esm2_t33_650M_UR50D', 'model_revision': '08e4846e537177426273712802403f7ba8261b6c', 'prediction_method': 'archived audited nested grouped predictions', 'arms': sorted(REFERENCE_ARMS)}
    else:
        if a.model_metadata is None:
            raise ValueError('model metadata is required')
        metadata = json.loads(a.model_metadata.read_text(encoding='utf-8-sig'))
        for k in ['model_id', 'model_revision', 'prediction_method', 'prediction_code_sha256', 'training_protocol']:
            if not isinstance(metadata.get(k), str) or not metadata[k].strip():
                raise ValueError(f'missing model metadata: {k}')
        code_hash = metadata['prediction_code_sha256']
        if len(code_hash) != 64 or any((c not in '0123456789abcdef' for c in code_hash)):
            raise ValueError('invalid prediction code hash')
        if metadata['training_protocol'] != 'frozen_nested_grouped_oof_v1':
            raise ValueError('unsupported training protocol')
        declared = metadata.get('arms', [metadata['model_id']])
        if not isinstance(declared, list) or not declared or any((not isinstance(x, str) or not x for x in declared)) or (len(set(declared)) != len(declared)):
            raise ValueError('invalid declared arms')
        metadata['arms'] = declared
    inputs = {str(task / 'manifest.json'): sha256(task / 'manifest.json'), str(pred): sha256(pred), str(Path(__file__).resolve()): sha256(Path(__file__).resolve()), str(Path(__file__).with_name('residue_scoring_core.py').resolve()): sha256(Path(__file__).with_name('residue_scoring_core.py'))}
    if a.model_metadata is not None:
        inputs[str(a.model_metadata.resolve())] = sha256(a.model_metadata)
    run = {'started_at': stamp(), 'status': 'RUNNING', 'model_metadata': metadata, 'input_sha256': inputs, 'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'scikit_learn': sklearn.__version__}, 'scope': 'software_example' if a.example else 'retrospective_grouped_development', 'warnings': ['Scores verify identities and metrics; training history is supplied provenance.', 'Intervals condition on fixed predictions; no refitting uncertainty.']}
    out.mkdir(parents=True)
    try:
        samples, expected, folds, data_hashes = load_task(task, example=a.example)
        run['dataset_sha256'] = data_hashes
        scores = read_scores(pred, samples, expected, folds, set(metadata['arms']), metadata['model_id'])
        rows, metrics = calculate(samples, expected, folds, scores)
        write_tsv(out / 'protein_metrics.tsv', list(rows[0]), rows)
        (out / 'metrics.json').write_bytes((json.dumps(metrics, indent=2, sort_keys=True) + '\n').encode('utf-8'))
        if load_task(task, example=a.example)[3] != data_hashes or any((sha256(Path(p)) != h for p, h in inputs.items())):
            raise ValueError('input or code changed during scoring')
        run.update(status='COMPLETE', ended_at=stamp(), n_predictions=sum((len(ps) for arm in scores.values() for ps in arm.values())), output_sha256={p.name: sha256(p) for p in out.iterdir() if p.is_file()})
    except Exception as e:
        run.update(status='FAILED', ended_at=stamp(), error=str(e))
        (out / 'run.json').write_bytes((json.dumps(run, indent=2, sort_keys=True) + '\n').encode('utf-8'))
        raise
    (out / 'run.json').write_bytes((json.dumps(run, indent=2, sort_keys=True) + '\n').encode('utf-8'))
    print(json.dumps({'status': run['status'], 'n_predictions': run['n_predictions'], 'output_dir': str(out)}, indent=2))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
