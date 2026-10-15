from __future__ import annotations
import argparse, ast, csv, gzip, io, json, tokenize
from collections import Counter
from pathlib import Path
try:
    from .verify_legacy_snapshot import *
    from .verify_legacy_snapshot import verify as verify_legacy
    from .evaluate_residue_predictions import load_task, check_manifest
except ImportError:
    from verify_legacy_snapshot import *
    from verify_legacy_snapshot import verify as verify_legacy
    from evaluate_residue_predictions import load_task, check_manifest

def verify_stable_controls(root):
    check_manifest(root)
    endpoints = {r['endpoint_sample_id']: r for r in csv.DictReader((root / 'endpoint_provenance.tsv').open(), delimiter='\t')}
    counts = {}
    tables = {}
    for name in ['labels_primary.tsv.gz', 'labels_strict_sensitivity.tsv.gz']:
        seen = set()
        c = Counter()
        values = {}
        with gzip.open(root / name, 'rt') as f:
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

def verify_no_comments(root):
    n = 0
    for p in sorted(root.rglob('*.py')):
        if any((x in p.parts for x in ['build', '.venv', '__pycache__'])):
            continue
        source = p.read_text()
        tree = ast.parse(source)
        if any((t.type == tokenize.COMMENT for t in tokenize.generate_tokens(io.StringIO(source).readline))):
            raise ValueError(f'Python comment found: {p.name}')
        if any((ast.get_docstring(x) for x in ast.walk(tree) if isinstance(x, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)))):
            raise ValueError(f'docstring found: {p.name}')
        n += 1
    return n

def verify_delivery_file_set(root):
    manifest = json.loads((root / 'manifest.json').read_text())
    declared = set(manifest['files'])
    actual = set()
    for p in root.rglob('*'):
        rel = p.relative_to(root)
        if any((part in {'.venv', '.git', 'build', 'dist', 'runs', 'wheelhouse', '__pycache__'} or part.endswith('.egg-info') or part.startswith('.venv_') for part in rel.parts)):
            continue
        if p.is_symlink():
            raise ValueError('symlink in delivery tree')
        if p.is_file() and rel.as_posix() != 'manifest.json':
            actual.add(rel.as_posix())
    if actual != declared:
        raise ValueError(f'delivery file set differs: unregistered={sorted(actual - declared)}; missing={sorted(declared - actual)}')
    return len(actual)

def verify(root):
    check_manifest(root)
    registered_files = verify_delivery_file_set(root)
    result = verify_legacy(root)
    samples, labels, folds, hashes = load_task(root / 'data/disorder_reference_v1')
    if len(samples) != 1279 or sum(map(len, labels.values())) != 391718 or len({s['binding_group_id'] for s in samples.values()}) != 1029:
        raise ValueError('disorder reference denominator changed')
    check_manifest(root / 'results/disorder_reference_esm2_v1')
    result['status'] = 'PASS_LOCAL_DELIVERY'
    result['registered_delivery_files'] = registered_files
    result['disorder_reference'] = {'proteins': 1279, 'unique_sequences': len({s['sequence_id'] for s in samples.values()}), 'binding_groups': 1029, 'labeled_residue_rows': 391718, 'outer_folds': 5, 'inner_folds': 3}
    result['fs_region_controls'] = verify_stable_controls(root / 'data/fs_region_controls_v1')
    result['python_files_without_comments_or_docstrings'] = verify_no_comments(root)
    return result

def main():
    p = argparse.ArgumentParser(description='Verify the complete local first-release bundle.')
    p.add_argument('--benchmark-root', required=True, type=Path)
    a = p.parse_args()
    print(json.dumps(verify(a.benchmark_root.resolve()), indent=2, sort_keys=True))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
