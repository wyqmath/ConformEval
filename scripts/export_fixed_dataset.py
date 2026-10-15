from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import tempfile
from pathlib import Path
COPY_FILES = ['task_samples.tsv.gz', 'sample_sequences.tsv.gz', 'sequences.fasta', 'position_labels_strict.tsv.gz', 'position_labels_structure_reference_candidate.tsv.gz', 'fs_region_full_length_labels.tsv.gz', 'case_modules.tsv', 'case_samples.tsv.gz', 'folds/fold_qc.tsv', 'folds/outer_group_holdout.tsv', 'folds/inner_validation_groups.tsv.gz']
SOURCE_ARTIFACTS = {'results/data_completion/20261012/02_fold_switch/final/cases.tsv': 'fold_switch_cases', 'results/data_completion/20261012/03_fs_regions/final/labels.tsv.gz': 'fold_switch_region_labels', 'results/data_completion/20261012/05_conditions/final/pairs.tsv': 'condition_pairs', 'results/benchmark_fixed/20261014/groups_v03/matched_constraints.tsv': 'matched_control_audit', 'configs/physinfobench_existing_20261012.json': 'accepted_protocol_stop_record'}
PATH_MARKERS = re.compile('(?:/Users/|/tmp/|/private/var/|/home/|[A-Za-z]:\\\\)')
LICENSE_STATUS = {'db_rcsb': 'CC0-1.0; attribution encouraged; PDB data files only', 'db_uniprot': 'CC-BY-4.0; attribution required', 'db_disprot': 'CC-BY-4.0; attribution required', 'db_knotprot': 'UNRESOLVED; source terms not located; do not redistribute dependent labels until cleared or rebuilt from licensed PDB-only inputs', 'model_esm2_t33_650m_p603': 'MIT checkpoint license; weights and cached representations are not bundled'}

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()

def read_tsv(path: Path, compressed: bool=False) -> list[dict[str, str]]:
    opener = gzip.open if compressed else open
    with opener(path, 'rt', newline='', encoding='utf-8') as handle:
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

def portable_case_sources(case_modules_path: Path, project_root: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    case_rows = read_tsv(case_modules_path)
    by_path: dict[str, Path] = {}
    for row in case_rows:
        raw = row.get('source_path', '')
        if raw not in SOURCE_ARTIFACTS:
            raise ValueError(f'unrecognized case source artifact: {raw}')
        path = Path(raw)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError(f'unsafe project source path: {raw}')
        resolved = project_root / path
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        by_path[raw] = resolved
    fold_cases_path = project_root / 'results/data_completion/20261012/02_fold_switch/final/cases.tsv'
    fold_evidence_path = project_root / 'results/data_completion/20261012/02_fold_switch/final/case_experimental_evidence.tsv'
    region_labels_path = project_root / 'results/data_completion/20261012/03_fs_regions/final/labels.tsv.gz'
    condition_pairs_path = project_root / 'results/data_completion/20261012/05_conditions/final/pairs.tsv'
    fold_cases = {row['case_id']: row for row in read_tsv(fold_cases_path)}
    fold_evidence = {row['evidence_id']: row for row in read_tsv(fold_evidence_path)}
    region_evidence_ids: dict[str, set[str]] = {}
    for row in read_tsv(region_labels_path, compressed=True):
        if row.get('evidence_id') not in {'', 'NA'}:
            region_evidence_ids.setdefault(row['target_id'], set()).add(row['evidence_id'])
    condition_pairs = {row['pair_id']: row for row in read_tsv(condition_pairs_path)}
    path_hashes = {raw: sha256(path) for raw, path in by_path.items()}
    portable_rows: list[dict[str, str]] = []
    for index, row in enumerate(sorted(case_rows, key=lambda item: (item['module'], item['case_id'])), start=1):
        raw = row['source_path']
        artifact_id = SOURCE_ARTIFACTS[raw]
        case_id = row['case_id']
        evidence_ids: set[str] = set()
        doi_values: set[str] = set()
        pmcid_values: set[str] = set()
        other_refs: set[str] = set()
        if artifact_id == 'fold_switch_cases':
            case = fold_cases.get(case_id)
            if case is None:
                raise ValueError(f'case missing from fold-switch source table: {case_id}')
            evidence_id = case.get('conversion_evidence_id', '')
            if evidence_id not in {'', 'NA'}:
                evidence_ids.add(evidence_id)
                evidence = fold_evidence.get(evidence_id)
                if evidence:
                    if evidence.get('doi'):
                        doi_values.add(evidence['doi'])
                    if evidence.get('pmcid'):
                        pmcid_values.add(evidence['pmcid'])
        elif artifact_id == 'fold_switch_region_labels':
            evidence_ids.update(region_evidence_ids.get(case_id, set()))
            for evidence_id in evidence_ids:
                evidence = fold_evidence.get(evidence_id)
                if evidence:
                    if evidence.get('doi'):
                        doi_values.add(evidence['doi'])
                    if evidence.get('pmcid'):
                        pmcid_values.add(evidence['pmcid'])
        elif artifact_id == 'condition_pairs':
            pair = condition_pairs.get(case_id)
            if pair is None:
                raise ValueError(f'condition pair missing from source table: {case_id}')
            reference = pair.get('independent_group', '').strip()
            if reference:
                other_refs.add(reference)
        else:
            other_refs.add('internal audit or accepted protocol record; see case_modules detail')
        source_ref = f'CS{index:04d}'
        row['source_ref'] = source_ref
        row.pop('source_path', None)
        portable_rows.append({'source_ref': source_ref, 'module': row['module'], 'case_id': case_id, 'source_artifact_id': artifact_id, 'source_artifact_filename': Path(raw).name, 'source_artifact_sha256': path_hashes[raw], 'evidence_ids': ';'.join(sorted(evidence_ids)), 'doi': ';'.join(sorted(doi_values)), 'pmcid': ';'.join(sorted(pmcid_values)), 'other_source_references': ';'.join(sorted(other_refs))})
    return (case_rows, portable_rows)

def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier = None
    sequence: list[str] = []
    for line in path.read_text(encoding='utf-8').splitlines():
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

def public_representation_coverage(data_dir: Path, overlay_path: Path) -> list[dict[str, str]]:
    base = read_tsv(data_dir / 'cache_coverage.tsv.gz', compressed=True)
    overlay = read_tsv(overlay_path, compressed=True)
    sequence_lengths = {(row['task_id'], row['sample_id'], row['sequence_sha256']): row['length'] for row in read_tsv(data_dir / 'sample_sequences.tsv.gz', compressed=True)}
    keys = lambda r: (r['task_id'], r['sample_id'], r['sequence_sha256'])
    result: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in base:
        key = keys(row)
        result[key] = {'task_id': row['task_id'], 'sample_id': row['sample_id'], 'sequence_sha256': row['sequence_sha256'], 'length': sequence_lengths.get(key, row.get('required_residues', '')), 'layer_set': row.get('layer_set', ''), 'representation_shape': row.get('mean_shape', ''), 'feature_ready': row.get('feature_ready', '0'), 'coverage_source': 'data_v06_cache_index'}
    for row in overlay:
        key = keys(row)
        current = result.get(key)
        if current is None:
            raise ValueError(f'representation overlay identity absent from frozen coverage: {key}')
        current.update({'length': row['length'], 'layer_set': row['layer_set'], 'representation_shape': row['mean_shape'], 'feature_ready': row['feature_ready_after_supplement'], 'coverage_source': row['source']})
    fields = ['task_id', 'sample_id', 'sequence_sha256', 'length', 'layer_set', 'representation_shape', 'feature_ready', 'coverage_source']
    return [result[key] for key in sorted(result)]

def export(data_dir: Path, data_manifest: Path, overlay_path: Path, overlay_manifest: Path, source_registry: Path, project_root: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise ValueError(f'refusing to overwrite existing output: {output_dir}')
    for path in (data_dir, data_manifest, overlay_path, overlay_manifest, source_registry, project_root):
        if not path.exists():
            raise FileNotFoundError(path)
    manifest = json.loads(data_manifest.read_text(encoding='utf-8'))
    overlay_meta = json.loads(overlay_manifest.read_text(encoding='utf-8'))
    for rel, expected in manifest.get('outputs', {}).items():
        path = data_dir / rel
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f'frozen data output missing or changed: {rel}')
    overlay_outputs = overlay_meta.get('outputs', {})
    overlay_rel = 'feature_coverage_overlay.tsv.gz'
    if overlay_rel not in overlay_outputs or sha256(overlay_path) != overlay_outputs[overlay_rel]:
        raise ValueError('independently reviewed representation overlay changed')
    sequence_audit = verify_sequences(data_dir)
    mask_audit = verify_masks(data_dir)
    coordinate_audit = verify_coordinates(data_dir)
    representation_rows = public_representation_coverage(data_dir, overlay_path)
    source_rows = read_tsv(source_registry)
    source_fields = ['source_id', 'name', 'doi_or_url', 'related_tasks', 'version', 'license', 'distribution_status']
    public_sources = []
    for row in source_rows:
        source_id = row['source_id']
        license_value = row.get('license', '') or 'not_recorded'
        status = LICENSE_STATUS.get(source_id)
        if status is None:
            if source_id.startswith('paper_') or source_id in {'legacy_paper_sources', 'legacy_manifests'}:
                status = 'SOURCE_SPECIFIC_REVIEW_REQUIRED; no full text is included'
            elif row.get('category') == 'model_host':
                status = 'weights_not_included; verify checkpoint license before each model download'
            elif row.get('category') == 'exclusion':
                continue
            else:
                status = 'REVIEW_REQUIRED'
        public_sources.append({'source_id': source_id, 'name': row.get('name', ''), 'doi_or_url': row.get('doi_or_url', ''), 'related_tasks': row.get('related_tasks', ''), 'version': row.get('version', ''), 'license': license_value, 'distribution_status': status})
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=f'.{output_dir.name}.', dir=output_dir.parent))
    try:
        for rel in COPY_FILES:
            source = data_dir / rel
            if not source.is_file():
                raise FileNotFoundError(source)
            target = temp / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
        portable_cases, case_sources = portable_case_sources(temp / 'case_modules.tsv', project_root)
        case_fields = ['module', 'case_id', 'status', 'qualification', 'sample_count', 'binding_groups', 'source_ref', 'detail']
        write_tsv(temp / 'case_modules.tsv', case_fields, portable_cases)
        case_source_fields = ['source_ref', 'module', 'case_id', 'source_artifact_id', 'source_artifact_filename', 'source_artifact_sha256', 'evidence_ids', 'doi', 'pmcid', 'other_source_references']
        write_tsv(temp / 'case_sources.tsv', case_source_fields, case_sources)
        protocol = json.loads((data_dir / 'protocol.json').read_text(encoding='utf-8'))
        protocol['protocol_source'] = 'accepted PhysInfoBench data protocol'
        protocol['source_data_lock_sha256'] = sha256(data_manifest)
        protocol['group_map'] = 'binding_group_id and fold assignments are included in this data release; historical identity graph remains in the source project'
        (temp / 'data_protocol.json').write_text(json.dumps(protocol, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        write_gzip_tsv(temp / 'representation_coverage.tsv.gz', ['task_id', 'sample_id', 'sequence_sha256', 'length', 'layer_set', 'representation_shape', 'feature_ready', 'coverage_source'], representation_rows)
        write_tsv(temp / 'sources.tsv', source_fields, public_sources)
        result = {'schema_version': '1.0', 'data_version': '2026-10-14-data-lock-v01-plus-FD03-coverage-overlay-paths-v01', 'distribution_status': 'DRAFT_NOT_CLEARED_FOR_PUBLIC_RELEASE', 'source_data_lock_sha256': sha256(data_manifest), 'representation_overlay_sha256': sha256(overlay_path), 'sequences': sequence_audit, 'residue_label_masks': mask_audit, 'coordinates': coordinate_audit, 'representation_coverage_rows': len(representation_rows), 'files': {}, 'excluded_local_artifacts': ['Original data manifest and QC include machine-specific paths.', 'Project-local case source paths are replaced by source_ref keys resolved in case_sources.tsv; upstream artifact hashes and recorded citation identifiers are retained.', 'Embedding caches and model weights are not included.', 'Raw article full texts and PDB coordinate files are not included.', 'The frozen data_v06 task_status.tsv is not copied because it predates FD.03; the benchmark-level task_catalog.tsv is authoritative for this draft.']}
        (temp / 'README.md').write_text('# PhysInfoBench data snapshot v0.3.1\n\nThis is a portable export of the frozen 2026-10-14 development dataset plus the independently reviewed FD.03 representation-coverage overlay. Labels, masks, sample identities, and folds are copied without modification. Case-module rows point to `case_sources.tsv` through `source_ref`; that table records the upstream artifact fingerprint and available DOI, PMCID, or source notes without exposing project-local paths. `representation_coverage.tsv.gz` contains no local cache paths.\n\nDistribution status: **draft, not cleared for public release**. Consult `sources.tsv` and the benchmark-level distribution notes before copying these files to a public host. No model weights or embedding caches are included.\n', encoding='utf-8')
        for path in temp.rglob('*'):
            if not path.is_file():
                continue
            opener = gzip.open if path.suffix.lower() == '.gz' else open
            try:
                with opener(path, 'rt', encoding='utf-8') as handle:
                    if any((PATH_MARKERS.search(line) for line in handle)):
                        raise ValueError(f'machine-specific path in exported file: {path.relative_to(temp)}')
            except UnicodeDecodeError:
                continue
        files = {}
        for path in sorted((p for p in temp.rglob('*') if p.is_file() and p.name != 'dataset_manifest.json')):
            rel = path.relative_to(temp).as_posix()
            files[rel] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
        result['files'] = files
        (temp / 'dataset_manifest.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        os.rename(temp, output_dir)
        return {'output_dir': str(output_dir), 'files': len(files), 'distribution_status': result['distribution_status'], 'sequences': sequence_audit}
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise

def main() -> int:
    parser = argparse.ArgumentParser(description='Export the frozen PhysInfoBench dataset as a portable, path-free snapshot.')
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--data-manifest', type=Path, required=True)
    parser.add_argument('--overlay', type=Path, required=True)
    parser.add_argument('--overlay-manifest', type=Path, required=True)
    parser.add_argument('--source-registry', type=Path, required=True)
    parser.add_argument('--project-root', type=Path, required=True, help='project root used to fingerprint internal case source artifacts')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    result = export(args.data_dir, args.data_manifest, args.overlay, args.overlay_manifest, args.source_registry, args.project_root, args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
if __name__ == '__main__':
    raise SystemExit(main())
