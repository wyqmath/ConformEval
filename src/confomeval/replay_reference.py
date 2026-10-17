from __future__ import annotations
import argparse, csv, gzip, json, math, subprocess, sys
from pathlib import Path
from .evaluate_predictions import evaluate, write_gzip_tsv, sha256, now

def compare(a, b, trail=''):
    if isinstance(a, dict):
        if set(a) != set(b):
            raise ValueError(f'metric keys differ: {trail}')
        for k in a:
            compare(a[k], b[k], trail + '/' + k)
    elif isinstance(a, list):
        if len(a) != len(b):
            raise ValueError(f'metric lengths differ: {trail}')
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, trail + f'/{i}')
    elif isinstance(a, (int, float)):
        if not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError(f'metric differs: {trail}: {a} != {b}')
    elif a != b:
        raise ValueError(f'metric differs: {trail}')

def main():
    ap = argparse.ArgumentParser(description='Replay all archived prediction arms and compare frozen reference metrics; no encoder inference or reader refitting.')
    ap.add_argument('--benchmark-root', required=True, type=Path)
    ap.add_argument('--output-dir', required=True, type=Path)
    a = ap.parse_args()
    root = a.benchmark_root.resolve()
    out = a.output_dir.resolve()
    if out.exists():
        raise ValueError('output directory already exists')
    out.mkdir(parents=True)
    started = now()
    ref = root / 'results/esm2_reference'
    with gzip.open(ref / 'predictions.tsv.gz', 'rt', encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    arms = sorted({r['arm'] for r in rows})
    checks = []
    source_provenance = json.loads((ref / 'provenance.json').read_text(encoding='utf-8-sig'))
    source_configuration = json.loads((ref / 'configuration.json').read_text(encoding='utf-8-sig'))
    metadata = {'model_id': source_provenance['encoder']['id'], 'model_revision': source_provenance['encoder']['revision'], 'prediction_method': 'archived audited nested grouped prediction replay', 'prediction_code_sha256': sha256(Path(__file__)), 'training_protocol': 'frozen_nested_grouped_oof_v1', 'candidate_layers': source_configuration['encoder']['selected_candidates'], 'pooling': 'mean_over_residues_without_special_tokens'}
    meta = out / 'reference_metadata.json'
    meta.write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8', newline='\n')
    for arm in arms:
        preds = out / f'{arm}.tsv.gz'
        write_gzip_tsv(preds, [{k: r[k] for k in ['task_id', 'sample_id', 'outer_fold', 'score']} for r in rows if r['arm'] == arm], ['task_id', 'sample_id', 'outer_fold', 'score'])
        target = out / arm
        evaluate(root / 'data/base', preds, meta, root / 'configs/scoring_protocol.json', target)
        actual = json.loads((target / 'metrics.json').read_text(encoding='utf-8-sig'))['metrics']
        for task, m in actual.items():
            expected = json.loads((ref / f'metrics_{task}.json').read_text(encoding='utf-8-sig'))['arms'][arm]
            if task == 'KNOT-PRESENCE':
                compare(m['auroc'], expected['auroc'])
                compare(m['average_precision_secondary'], expected['average_precision_secondary'])
                compare(m['group_bootstrap']['ci95'], expected['binding_group_bootstrap_auroc_ci95'])
            else:
                for key, value in expected.items():
                    mapped = {'mean_rank_percentile_lower_is_better': 'mean_held_case_rank_percentile', 'binding_group_bootstrap_ci95': 'group_bootstrap_ci95'}.get(key, key)
                    compare(m[mapped], value, task + '/' + arm + '/' + key)
                with (ref / 'case_ranks.tsv').open(newline='', encoding='utf-8-sig') as f:
                    archived = {x['sample_id']: x for x in csv.DictReader(f, delimiter='\t') if x['arm'] == arm}
                if set(archived) != {x['sample_id'] for x in m['case_ranks']}:
                    raise ValueError('FS case-rank identities differ')
                for case in m['case_ranks']:
                    old = archived[case['sample_id']]
                    for field, value in case.items():
                        original = old['group_id' if field == 'binding_group_id' else field]
                        compare(value, float(original) if isinstance(value, (int, float)) else original, 'FS case-rank/' + arm + '/' + field)
            checks.append({'task_id': task, 'arm': arm, 'metrics_and_group_intervals': 'MATCH'})
        print(f'Replayed {arm}: both task metrics and group intervals MATCH', flush=True)
    task = root / 'data/disorder_reference'
    dest = out / 'disorder_reference'
    subprocess.run([sys.executable, '-m', 'confomeval.evaluate_residue_predictions', '--task-dir', str(task), '--predictions', str(root / 'results/disorder_reference_esm2/residue_predictions.tsv.gz'), '--reference', '--output-dir', str(dest)], check=True)
    actual = json.loads((dest / 'metrics.json').read_text(encoding='utf-8-sig'))
    expected = json.loads((root / 'results/disorder_reference_esm2/metrics.json').read_text(encoding='utf-8-sig'))
    compare(actual, expected, 'disorder_reference_all_metrics_and_bootstraps')
    checks.append({'task_id': 'DISORDER-STRUCTURE-REFERENCE-CANDIDATE', 'arms': 7, 'all_metrics_and_bootstraps': 'MATCH'})
    record = {'status': 'PASS', 'started_at': started, 'ended_at': now(), 'scope': 'Archived prediction rescoring; no encoder inference, reader refitting, or full inner-selection trajectory replay. KNOT/FS check point metrics, group intervals and all FS case ranks; disorder checks full metrics and 2000 group-bootstrap intervals.', 'numeric_tolerance': 1e-10, 'checks': checks}
    (out / 'REPLAY.json').write_text(json.dumps(record, indent=2, sort_keys=True) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(record, indent=2))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
