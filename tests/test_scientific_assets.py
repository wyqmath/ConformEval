"""Check the shipped scientific inputs independently of repository organization."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from conformeval import data_validation as validation
from conformeval.evaluate_predictions import verify_dataset
from conformeval.evaluate_residue_predictions import check_manifest, load_task

ROOT = Path(os.environ.get('CONFORMEVAL_REFERENCE_ROOT', Path(__file__).resolve().parents[1]))


@unittest.skipUnless((ROOT / 'data/base/dataset_manifest.json').is_file(),
                     'Set CONFORMEVAL_REFERENCE_ROOT to the extracted reference bundle')
class ScientificAssetTests(unittest.TestCase):

    def test_base_identities_masks_coordinates_and_groups(self):
        data = ROOT / 'data/base'
        verify_dataset(data)
        self.assertEqual(validation.verify_sequences(data)['fasta_sequences'], 11684)
        validation.verify_masks(data)
        validation.verify_coordinates(data)
        validation.verify_case_source_references(data)
        validation.verify_knot_representation_coverage(data)
        with (ROOT / 'data/task_catalog.tsv').open(newline='', encoding='utf-8-sig') as handle:
            tasks = {r['task_id']: r for r in csv.DictReader(handle, delimiter='\t')
                     if r['task_id'] in validation.EXPECTED}
        validation.verify_group_splits(data, tasks)

    def test_reference_prediction_coverage_and_assignment(self):
        predictions = ROOT / 'results/esm2_reference/predictions.tsv.gz'
        self.assertEqual(sum(validation.verify_predictions(predictions).values()), 67291)
        validation.verify_prediction_assignments(predictions, ROOT / 'data/base')

    def test_8m_predictions_use_one_combined_table_per_task(self):
        expected = {'knot': (734, 7), 'fs': (8879, 7)}
        for task, (rows_per_arm, arm_count) in expected.items():
            with self.subTest(task=task):
                path = ROOT / f'results/esm2_8m/{task}/predictions.tsv.gz'
                rows = validation.read_tsv(path, compressed=True)
                arms = {row['arm'] for row in rows}
                self.assertEqual(len(arms), arm_count)
                self.assertEqual(len(rows), rows_per_arm * arm_count)
                self.assertFalse(list(path.parent.glob('PLM-*.tsv.gz')))
                for arm in arms:
                    self.assertEqual(sum(row['arm'] == arm for row in rows), rows_per_arm)

    def test_disorder_reference_identities_and_folds(self):
        samples, labels, folds, hashes = load_task(ROOT / 'data/disorder_reference')
        self.assertEqual(len(samples), 1279)
        self.assertEqual(sum(map(len, labels.values())), 391718)
        self.assertEqual(set(samples), set(folds))

    def test_stable_control_labels(self):
        actual = validation.verify_stable_controls(ROOT / 'data/fs_region_controls')
        self.assertEqual(actual['endpoints'], 20)
        self.assertEqual(actual['primary'], {'1': 124, '0': 848, 'NA': 4746})

    def test_data_integrity_rejects_modified_artifact(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            payload = root / 'labels.tsv'
            payload.write_bytes(('label\tmask\n1\t1\n').encode('utf-8'))
            record = {'sha256': hashlib.sha256(payload.read_bytes()).hexdigest(),
                      'bytes': payload.stat().st_size}
            (root / 'manifest.json').write_bytes((json.dumps({'files': {'labels.tsv': record}})).encode('utf-8'))
            check_manifest(root)
            payload.write_bytes(('label\tmask\n0\t1\n').encode('utf-8'))
            with self.assertRaisesRegex(ValueError, 'hash/size mismatch'):
                check_manifest(root)

    def test_control_identity_mask_coordinate_and_coverage_failures(self):
        source = ROOT / 'data/fs_region_controls'
        for field, value, message in [
            ('sequence_sha256', 'a' * 64, 'endpoint identity'),
            ('mask', 'invalid', 'not masked'),
            ('position_1based', '999999', 'out-of-range'),
            (None, None, 'coverage differs'),
        ]:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temp:
                root = Path(temp) / 'controls'
                shutil.copytree(source, root)
                table = root / 'labels_primary.tsv.gz'
                with gzip.open(table, 'rt', newline='', encoding='utf-8-sig') as handle:
                    reader = csv.DictReader(handle, delimiter='\t')
                    fields, rows = reader.fieldnames, list(reader)
                if field is None:
                    rows.pop(0)
                else:
                    rows[0][field] = value
                validation.write_gzip_tsv(table, fields, rows)
                # Keep the data hash consistent to exercise semantic rejection.
                manifest_path = root / 'manifest.json'
                manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
                manifest['files'][table.name] = {
                    'sha256': hashlib.sha256(table.read_bytes()).hexdigest(),
                    'bytes': table.stat().st_size,
                }
                manifest_path.write_bytes((json.dumps(manifest)).encode('utf-8'))
                with self.assertRaisesRegex(ValueError, message):
                    validation.verify_stable_controls(root)


if __name__ == '__main__':
    unittest.main()
