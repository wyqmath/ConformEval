"""Extract sequence and residue representations without consuming task labels."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from .adapters import HFAdapter, download_snapshot, load_config
from .readers import parse_fasta
from .runtime import RunRecord, json_write


def extract(args, adapter_factory=HFAdapter):
    config = load_config(args.model_config)
    for name in ['device', 'dtype']:
        if getattr(args, name, None):
            config[name] = getattr(args, name)
    if getattr(args, 'sequence_only', False):
        config['store_residues'] = False
    record = RunRecord(args.output_dir.resolve(), config)
    output = record.output
    try:
        record.track(args.model_config)
        record.track(args.input_fasta)
        sequences = parse_fasta(args.input_fasta)
        if not sequences:
            raise ValueError('input FASTA is empty')
        import torch
        torch.set_num_threads(args.threads)
        snapshot = args.snapshot or download_snapshot(config, args.cache_dir, args.local_files_only)
        for path in Path(snapshot).iterdir():
            if path.is_file() and path.suffix in {'.json', '.safetensors', '.bin', '.model', '.txt'}:
                record.track(path)
        adapter = adapter_factory(config, Path(snapshot))
        record.data.update(model_depth=adapter.depth, output_dimension=adapter.dimension,
                           checkpoint_loading=adapter.loading)
        (output / 'archives').mkdir()
        rows = []
        for index, (sequence_id, sequence) in enumerate(sorted(sequences.items()), start=1):
            pooled, residues, segments = adapter.encode(sequence)
            meta = {'sequence_id': sequence_id, 'len_full': len(sequence), 'len_used': len(sequence),
                    'truncated': False, 'model': config['model_id'], 'revision': config['revision'],
                    'layer_set': 'mean_' + '_'.join(map(str, config['layers'])) + '+resid_' + '_'.join(map(str, config['layers'])),
                    'pooling': config['pooling'], 'dimension': adapter.dimension,
                    'extraction_config': config, 'normalization': 'UZOBJ_to_X' if config['adapter']=='prott5' else 'J_to_X',
                    'windows_0based_half_open': segments, 'special_tokens_excluded': True}
            path = output / 'archives' / f'sequence-{index:06d}.npz'
            arrays = {'meta': np.asarray(json.dumps(meta, sort_keys=True)), 'sequence_layers': pooled}
            if config['store_residues']:
                arrays.update(resid_layers=residues, resid_indices=np.arange(1, len(sequence)+1, dtype=np.int64))
            np.savez_compressed(path, **arrays)
            rows.append({'sequence_id': sequence_id, 'length': len(sequence), 'file': path.relative_to(output).as_posix()})
        manifest = {'schema_version': 1, 'model_id': config['model_id'], 'model_revision': config['revision'],
                    'output_dimension': adapter.dimension, 'candidate_layers': config['layers'],
                    'pooling': config['pooling'], 'extraction_config': config,
                    'archives': rows}
        json_write(output / 'representation_manifest.json', manifest)
        record.data.update(unique_sequences=len(rows), total_residues=sum(r['length'] for r in rows),
                           model_id=config['model_id'], model_revision=config['revision'])
        if config['device']=='cuda':
            record.data['peak_cuda_allocated_bytes'] = torch.cuda.max_memory_allocated()
    except Exception as error:
        record.finish(error)
        raise
    error = record.finish()
    if error:
        raise error
    print(json.dumps({'status':'PASS','manifest':str(output / 'representation_manifest.json'),
                      'unique_sequences':len(rows), 'dimension':adapter.dimension}))
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-fasta', required=True, type=Path)
    p.add_argument('--model-config', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--snapshot', type=Path)
    p.add_argument('--cache-dir', type=Path)
    p.add_argument('--local-files-only', action='store_true')
    p.add_argument('--sequence-only', action='store_true')
    p.add_argument('--device', choices=['cpu','cuda','mps'])
    p.add_argument('--dtype', choices=['float32','float16','bfloat16'])
    p.add_argument('--threads', type=int, default=4)
    args = p.parse_args()
    extract(args)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
