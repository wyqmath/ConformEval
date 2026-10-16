from __future__ import annotations
import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
from confomeval.run_nested_readers import fold_indices, fs_rank_percentiles, fit_logistic, nearest_neighbor, run, validate_sample_sequence_mappings
from confomeval.build_feature_manifest import scan_archives
TASK = 'KNOT-PRESENCE'

def write_inner(path: Path, rows: list[dict[str, str]]) -> None:
    fields = ['task_id', 'outer_fold', 'inner_fold', 'group_id']
    with gzip.open(path, 'wt', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def split_fixture(root: Path, contaminate: bool=False):
    fold_rows = []
    samples = []
    group_to_fold = {}
    for outer in range(5):
        for label, suffix in ((1, 'p'), (0, 'n')):
            group = f'G{outer}_{suffix}'
            group_to_fold[TASK, group] = str(outer)
            samples.append({'sample_id': group, 'binding_group_id': group, 'label': label})
        train_outer_folds = [group_outer for group_outer in range(5) if group_outer != outer]
        for inner_index, group_outer in enumerate(train_outer_folds):
            inner = str(inner_index % 3)
            for suffix in ('p', 'n'):
                fold_rows.append({'task_id': TASK, 'outer_fold': str(outer), 'inner_fold': inner, 'group_id': f'G{group_outer}_{suffix}'})
    if contaminate:
        fold_rows.append({'task_id': TASK, 'outer_fold': '0', 'inner_fold': '0', 'group_id': 'G0_p'})
    fold_dir = root / 'folds'
    fold_dir.mkdir(parents=True)
    write_inner(fold_dir / 'inner_validation_groups.tsv.gz', fold_rows)
    return (samples, group_to_fold)

class NestedReaderTests(unittest.TestCase):

    def write_representation(self, root: Path, name: str, digest: str, dimension: int=4, revision: str='r1', value: float=1.0) -> Path:
        path = root / name
        meta = {'seq_sha256': digest, 'len_full': 4, 'len_used': 4, 'truncated': False, 'model': 'test/model', 'revision': revision, 'layer_set': 'mean_17_33'}
        np.savez_compressed(path, meta=np.asarray(json.dumps(meta)), mean_layers=np.full((2, dimension), value, dtype=np.float32))
        return path

    def test_manifest_builder_scans_complete_exact_sequence_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            digest_a, digest_b = ('a' * 64, 'b' * 64)
            self.write_representation(root, 'a.npz', digest_a)
            self.write_representation(root, 'b.npz', digest_b)
            rows, dimension = scan_archives(root, {digest_a: (f'SEQ_{digest_a}', 4), digest_b: (f'SEQ_{digest_b}', 4)}, 'test/model', 'r1')
            self.assertEqual(len(rows), 2)
            self.assertEqual(dimension, 4)
            self.assertTrue(all((len(row['archive_sha256']) == 64 for row in rows)))

    def test_manifest_builder_rejects_missing_conflicting_or_mixed_dimension_cache(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            digest_a, digest_b = ('a' * 64, 'b' * 64)
            self.write_representation(root, 'a.npz', digest_a)
            with self.assertRaisesRegex(ValueError, 'missing 1 eligible sequences'):
                scan_archives(root, {digest_a: (f'SEQ_{digest_a}', 4), digest_b: (f'SEQ_{digest_b}', 4)}, 'test/model', 'r1')
            self.write_representation(root, 'b.npz', digest_b, dimension=5)
            with self.assertRaisesRegex(ValueError, 'dimensions differ'):
                scan_archives(root, {digest_a: (f'SEQ_{digest_a}', 4), digest_b: (f'SEQ_{digest_b}', 4)}, 'test/model', 'r1')
            self.write_representation(root, 'a.npz', digest_a, dimension=5)
            self.write_representation(root, 'a_duplicate.npz', digest_a, dimension=5, value=2.0)
            with self.assertRaisesRegex(ValueError, 'multiple different archives'):
                scan_archives(root, {digest_a: (f'SEQ_{digest_a}', 4), digest_b: (f'SEQ_{digest_b}', 4)}, 'test/model', 'r1')

    def test_inner_partition_is_complete_and_disjoint_from_outer_test(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows, assignment = split_fixture(root)
            labels, plan = fold_indices(TASK, rows, root, assignment)
            self.assertEqual(len(plan), 5)
            self.assertEqual(set(labels), {0, 1})
            for outer, (train, test, inners) in plan.items():
                self.assertFalse(set(train) & set(test))
                self.assertEqual(len(inners), 3)
                for _, inner_train, validation in inners:
                    self.assertFalse(set(inner_train) & set(validation))
                    self.assertFalse(set(test) & (set(inner_train) | set(validation)))

    def test_inner_split_rejects_outer_test_group_contamination(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows, assignment = split_fixture(root, contaminate=True)
            with self.assertRaisesRegex(ValueError, 'inner validation groups do not partition outer training'):
                fold_indices(TASK, rows, root, assignment)

    def test_preprocessing_is_fit_on_training_rows_only(self):
        features = np.asarray([[0.0, 2.0], [2.0, 4.0], [1000.0, -1000.0]])
        labels = np.asarray([0, 1, 0])
        scores, params = fit_logistic(features, labels, np.asarray([0, 1]), np.asarray([2]), 1.0, {'phase': 'test'}, [])
        self.assertEqual(scores.shape, (1,))
        np.testing.assert_allclose(params['scaler_mean'], [1.0, 3.0])

    def test_fs_rank_percentiles_use_average_ties(self):
        observed = fs_rank_percentiles(np.asarray([0.9, 0.9, 0.1, 0.2]), np.asarray([1, 0, 0, 0]))
        np.testing.assert_allclose(observed, [0.375])

    def test_nearest_reference_record_is_json_serializable(self):
        features = np.asarray([[0.0, 1.0], [1.0, 0.0], [0.2, 0.8]])
        labels = np.asarray([0, 1, 1])
        _, reference = nearest_neighbor(features, labels, np.asarray([0, 1]), np.asarray([2]), ['train-negative', 'train-positive', 'test'])
        restored = json.loads(json.dumps(reference, allow_nan=False))
        self.assertEqual(restored['reference_sample_ids'], ['train-negative', 'train-positive'])
        self.assertEqual(len(restored['scaler_mean']), features.shape[1])

    def test_identical_fold_switch_endpoints_keep_two_endpoint_roles(self):
        digest = 'a' * 64
        sample = {'sequence_ids': f'SEQ_{digest};SEQ_{digest}', 'sequence_sha256s': f'{digest};{digest}', 'lengths': '74;74', 'unit_type': 'positive_endpoint_pair'}
        mappings = [{'sequence_id': f'SEQ_{digest}', 'sequence_sha256': digest, 'length': '74', 'endpoint_role': 'endpoint_A'}, {'sequence_id': f'SEQ_{digest}', 'sequence_sha256': digest, 'length': '74', 'endpoint_role': 'endpoint_B'}]
        hashes, lengths = validate_sample_sequence_mappings(('FS-L2-grouped-PU-ranking', 'pair'), sample, mappings)
        self.assertEqual(hashes, [digest, digest])
        self.assertEqual(lengths, [74, 74])
        mappings[1]['endpoint_role'] = 'endpoint_A'
        with self.assertRaisesRegex(ValueError, 'one endpoint_A and one endpoint_B'):
            validate_sample_sequence_mappings(('FS-L2-grouped-PU-ranking', 'pair'), sample, mappings)

    def test_failed_run_leaves_timestamped_record_and_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / 'failed_run'
            args = SimpleNamespace(output_dir=output, data_dir=root / 'missing_data', protocol=root / 'missing_protocol.json', feature_manifest=root / 'missing_features.json')
            with self.assertRaisesRegex(FileNotFoundError, 'missing_data'):
                run(args)
            record = json.loads((output / 'run.json').read_text(encoding='utf-8'))
            self.assertEqual(record['status'], 'FAILED')
            self.assertEqual(record['error_type'], 'FileNotFoundError')
            self.assertIn('started_at', record)
            self.assertIn('completed_at', record)
            with self.assertRaisesRegex(ValueError, 'refusing to overwrite'):
                run(args)
if __name__ == '__main__':
    unittest.main()
