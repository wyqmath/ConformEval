from __future__ import annotations
import csv
import gzip
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from confomeval import data_validation
ROOT = Path(__file__).resolve().parents[1]


class DataIntegrityTests(unittest.TestCase):

    def write_tsv(self, path: Path, fields: list[str], rows: list[dict], compressed: bool=False):
        path.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if compressed else open
        mode = 'wt'
        with opener(path, mode, newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

    def fixture(self, root: Path):
        sequence = 'ACD'
        digest = hashlib.sha256(sequence.encode('ascii')).hexdigest()
        (root / 'sequences.fasta').write_text(f'>SEQ_{digest}\n{sequence}\n', encoding='utf-8')
        self.write_tsv(root / 'sample_sequences.tsv.gz', ['task_id', 'sample_id', 'endpoint_role', 'sequence_id', 'sequence_sha256', 'length'], [{'task_id': 'FS-REGION', 'sample_id': 'case_A', 'endpoint_role': 'A', 'sequence_id': f'SEQ_{digest}', 'sequence_sha256': digest, 'length': '3'}], compressed=True)
        label_fields = ['sample_id', 'sequence_sha256', 'position', 'label', 'mask']
        row = {'sample_id': 'DP00001', 'sequence_sha256': digest, 'position': '1', 'label': '1', 'mask': '1'}
        self.write_tsv(root / 'position_labels_strict.tsv.gz', label_fields, [row], compressed=True)
        self.write_tsv(root / 'position_labels_structure_reference_candidate.tsv.gz', label_fields, [row], compressed=True)
        self.write_tsv(root / 'fs_region_full_length_labels.tsv.gz', ['endpoint_sample_id', 'sequence_position', 'start', 'end', 'target_value', 'mask'], [{'endpoint_sample_id': 'case_A', 'sequence_position': '1', 'start': '1', 'end': '1', 'target_value': '1', 'mask': '1'}], compressed=True)
        return digest

    def test_unknown_mask_requires_unknown_value(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fields = ['sample_id', 'sequence_sha256', 'position', 'label', 'mask']
            row = {'sample_id': 'DP00001', 'sequence_sha256': 'abc', 'position': '1', 'label': '0', 'mask': '0'}
            self.write_tsv(root / 'position_labels_strict.tsv.gz', fields, [row], compressed=True)
            self.write_tsv(root / 'position_labels_structure_reference_candidate.tsv.gz', fields, [], compressed=True)
            with self.assertRaisesRegex(ValueError, 'unknown residue label'):
                data_validation.verify_masks(root)

    def test_residue_and_fs_region_coordinates_are_bounded(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            digest = self.fixture(root)
            data_validation.verify_coordinates(root)
            fields = ['sample_id', 'sequence_sha256', 'position', 'label', 'mask']
            bad = {'sample_id': 'DP00001', 'sequence_sha256': digest, 'position': '4', 'label': '1', 'mask': '1'}
            self.write_tsv(root / 'position_labels_strict.tsv.gz', fields, [bad], compressed=True)
            with self.assertRaisesRegex(ValueError, 'outside sequence'):
                data_validation.verify_coordinates(root)

    def test_gzip_tsv_is_deterministic(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fields = ['id', 'label']
            rows = [{'id': 'a', 'label': '1'}, {'id': 'b', 'label': 'NA'}]
            data_validation.write_gzip_tsv(root / 'a.tsv.gz', fields, rows)
            data_validation.write_gzip_tsv(root / 'b.tsv.gz', fields, rows)
            self.assertEqual((root / 'a.tsv.gz').read_bytes(), (root / 'b.tsv.gz').read_bytes())

class ReferenceDataTests(unittest.TestCase):

    def write_tsv(self, path: Path, fields: list[str], rows: list[dict], compressed: bool=False):
        path.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if compressed else open
        with opener(path, 'wt', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)

    def test_catalog_declares_denominator_units_and_consistent_case_counts(self):
        result = data_validation.verify_task_catalog(ROOT)
        self.assertEqual(result['registered_tasks'], 13)
        with (ROOT / 'data/task_catalog.tsv').open(newline='', encoding='utf-8') as handle:
            rows = {row['task_id']: row for row in csv.DictReader(handle, delimiter='\t')}
        self.assertEqual(rows['FS-L1-PAIRED']['candidate_units'], '16')
        self.assertEqual(rows['FS-L1-PAIRED']['eligible_units'], '10')
        self.assertEqual(rows['FS-REGION']['candidate_units'], '10')
        self.assertEqual(rows['FS-REGION']['eligible_units'], '3')
        self.assertEqual(rows['FS-REGION']['eligible_count_unit'], 'case_with_supported_region')

    def test_prediction_assignment_audit_joins_samples_folds_and_all_arms(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            data = root / 'data'
            (data / 'folds').mkdir(parents=True)
            self.write_tsv(data / 'task_samples.tsv.gz', ['task_id', 'sample_id', 'data_eligible', 'binding_group_id', 'learning_label', 'learning_mask', 'biological_target', 'biological_mask', 'label', 'mask'], [{'task_id': 'TASK', 'sample_id': 'S0', 'data_eligible': '1', 'binding_group_id': 'G0', 'learning_label': '1', 'learning_mask': '1', 'biological_target': '1', 'biological_mask': '1', 'label': '1', 'mask': '1'}, {'task_id': 'TASK', 'sample_id': 'S1', 'data_eligible': '1', 'binding_group_id': 'G1', 'learning_label': '0', 'learning_mask': '0', 'biological_target': 'null', 'biological_mask': '0', 'label': 'NA', 'mask': '0'}], compressed=True)
            self.write_tsv(data / 'folds/outer_group_holdout.tsv', ['task_id', 'outer_fold', 'group_id', 'role'], [{'task_id': 'TASK', 'outer_fold': '0', 'group_id': 'G0', 'role': 'test'}, {'task_id': 'TASK', 'outer_fold': '1', 'group_id': 'G1', 'role': 'test'}])
            fields = ['task_id', 'outer_fold', 'arm', 'sample_id', 'binding_group_id', 'learning_label', 'learning_mask', 'biological_target', 'biological_mask', 'label', 'mask']
            predictions = []
            for sample, fold, group, learning_label, learning_mask, target, bio_mask, label, mask in (('S0', '0', 'G0', '1', '1', '1', '1', '1', '1'), ('S1', '1', 'G1', '0', '0', 'null', '0', 'NA', '0')):
                for arm in ('A', 'B'):
                    predictions.append({'task_id': 'TASK', 'outer_fold': fold, 'arm': arm, 'sample_id': sample, 'binding_group_id': group, 'learning_label': learning_label, 'learning_mask': learning_mask, 'biological_target': target, 'biological_mask': bio_mask, 'label': label, 'mask': mask})
            path = root / 'predictions.tsv.gz'
            self.write_tsv(path, fields, predictions, compressed=True)
            audit = data_validation.verify_prediction_assignments(path, data, {'TASK'}, {'A', 'B'})
            self.assertEqual(audit['TASK'], {'eligible_samples': 2, 'arms_per_sample': 2})
            wrong_fold = [dict(row) for row in predictions]
            wrong_fold[0]['outer_fold'] = '1'
            self.write_tsv(path, fields, wrong_fold, compressed=True)
            with self.assertRaisesRegex(ValueError, 'frozen outer test fold'):
                data_validation.verify_prediction_assignments(path, data, {'TASK'}, {'A', 'B'})
            missing_arm = predictions[:-1]
            self.write_tsv(path, fields, missing_arm, compressed=True)
            with self.assertRaisesRegex(ValueError, 'arm coverage differs'):
                data_validation.verify_prediction_assignments(path, data, {'TASK'}, {'A', 'B'})

    def test_group_split_audit_accepts_partition_and_rejects_outer_leakage(self):
        eligible = {'TASK': {'a', 'b', 'c', 'd'}}
        outer = [{'task_id': 'TASK', 'outer_fold': '0', 'group_id': 'a', 'role': 'test'}, {'task_id': 'TASK', 'outer_fold': '0', 'group_id': 'b', 'role': 'test'}, {'task_id': 'TASK', 'outer_fold': '1', 'group_id': 'c', 'role': 'test'}, {'task_id': 'TASK', 'outer_fold': '1', 'group_id': 'd', 'role': 'test'}]
        inner = []
        for fold, groups in (('0', ('c', 'd')), ('1', ('a', 'b'))):
            for inner_fold, group in enumerate(groups):
                inner.append({'task_id': 'TASK', 'outer_fold': fold, 'inner_fold': str(inner_fold), 'group_id': group, 'role': 'validation'})
        result = data_validation.audit_group_splits(eligible, outer, inner, {'TASK': 2}, inner_fold_count=2)
        self.assertEqual(result['TASK']['eligible_groups'], 4)
        inner[0]['group_id'] = 'a'
        with self.assertRaisesRegex(ValueError, 'outer test group leaked'):
            data_validation.audit_group_splits(eligible, outer, inner, {'TASK': 2}, inner_fold_count=2)

    def test_fs_l2_unlabeled_background_must_remain_masked(self):
        target_globals = data_validation.verify_predictions.__globals__
        old_expected = target_globals['EXPECTED']
        target_globals['EXPECTED'] = {'FS-L2-grouped-PU-ranking': 1}
        try:
            with tempfile.TemporaryDirectory() as temp:
                path = Path(temp) / 'predictions.tsv.gz'
                fields = ['task_id', 'outer_fold', 'arm', 'sample_id', 'learning_label', 'biological_target', 'biological_mask', 'label', 'mask', 'score']
                valid = {'task_id': 'FS-L2-grouped-PU-ranking', 'outer_fold': '0', 'arm': 'PLM', 'sample_id': 'U1', 'learning_label': '0', 'biological_target': 'null', 'biological_mask': '0', 'label': 'NA', 'mask': '0', 'score': '0.2'}
                with gzip.open(path, 'wt', newline='', encoding='utf-8') as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
                    writer.writeheader()
                    writer.writerow(valid)
                data_validation.verify_predictions(path)
                valid['biological_target'] = '0'
                with gzip.open(path, 'wt', newline='', encoding='utf-8') as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
                    writer.writeheader()
                    writer.writerow(valid)
                with self.assertRaisesRegex(ValueError, 'converted to a biological negative'):
                    data_validation.verify_predictions(path)
        finally:
            target_globals['EXPECTED'] = old_expected



    def test_case_source_references_resolve_and_reject_broken_links(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fields = ['module', 'case_id', 'source_ref', 'detail']
            source_ref = 'CS0001'
            self.write_tsv(root / 'case_modules.tsv', fields, [{'module': 'FS-L1-PAIRED', 'case_id': 'case_A', 'source_ref': source_ref, 'detail': 'audited'}])
            self.write_tsv(root / 'case_sources.tsv', ['source_ref', 'module', 'case_id', 'source_artifact_sha256'], [{'source_ref': source_ref, 'module': 'FS-L1-PAIRED', 'case_id': 'case_A', 'source_artifact_sha256': 'a' * 64}])
            self.assertEqual(data_validation.verify_case_source_references(root), {'case_rows': 1, 'resolvable_source_records': 1})
            self.write_tsv(root / 'case_modules.tsv', fields, [{'module': 'FS-L1-PAIRED', 'case_id': 'case_A', 'source_ref': 'CS9999', 'detail': 'audited'}])
            with self.assertRaisesRegex(ValueError, 'does not resolve'):
                data_validation.verify_case_source_references(root)

    def test_case_source_verifier_rejects_project_relative_source_path(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.write_tsv(root / 'case_modules.tsv', ['module', 'case_id', 'source_path'], [{'module': 'FS-L1-PAIRED', 'case_id': 'case_A', 'source_path': 'results/private.tsv'}])
            self.write_tsv(root / 'case_sources.tsv', ['source_ref', 'module', 'case_id', 'source_artifact_sha256'], [{'source_ref': 'CS0001', 'module': 'FS-L1-PAIRED', 'case_id': 'case_A', 'source_artifact_sha256': 'a' * 64}])
            with self.assertRaisesRegex(ValueError, 'source_path field remains'):
                data_validation.verify_case_source_references(root)
if __name__ == '__main__':
    unittest.main()
