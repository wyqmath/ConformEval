import json
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
from confomeval.model_interface import align_residue_labels, load_mean_representation, load_residue_representation, sha256_file
MODEL = 'facebook/esm2_t33_650M_UR50D'
REVISION = 'test-revision'
SEQ_SHA = 'a' * 64
MEAN_LAYERS = (5, 11, 17, 23, 29, 33)
RESIDUE_LAYERS = (11, 23, 33)

class ConformEvalModelInterfaceTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def save_archive(self, name, **arrays):
        path = self.root / name
        np.savez_compressed(path, **arrays)
        return path

    def meta(self, layer_set, length=4, **overrides):
        values = {'seq_sha256': SEQ_SHA, 'len_full': length, 'len_used': length, 'truncated': False, 'model': MODEL, 'revision': REVISION, 'layer_set': layer_set}
        values.update(overrides)
        return np.asarray(json.dumps(values))

    def test_mean_interface_selects_declared_layers_and_checks_identity(self):
        path = self.save_archive('mean.npz', meta=self.meta('mean_5_11_17_23_29_33'), mean_layers=np.arange(6 * 3, dtype=np.float16).reshape(6, 3))
        selected, dimension = load_mean_representation(path, SEQ_SHA, 4, sha256_file(path), (17, 5), MODEL, REVISION)
        self.assertEqual(selected.shape, (2, 3))
        self.assertEqual(dimension, 3)
        np.testing.assert_array_equal(selected[0], np.asarray([6, 7, 8], dtype=np.float32))
        np.testing.assert_array_equal(selected[1], np.asarray([0, 1, 2], dtype=np.float32))
        with self.assertRaisesRegex(ValueError, 'sequence hash mismatch'):
            load_mean_representation(path, 'b' * 64, 4, sha256_file(path), (17,), MODEL, REVISION)

    def test_mean_interface_rejects_missing_archive_hash_and_nonfinite_values(self):
        path = self.save_archive('mean_bad.npz', meta=self.meta('mean_5_11_17_23_29_33'), mean_layers=np.full((6, 3), np.nan, dtype=np.float32))
        with self.assertRaisesRegex(ValueError, 'archive SHA-256 is required'):
            load_mean_representation(path, SEQ_SHA, 4, '', (17,), MODEL, REVISION)
        with self.assertRaisesRegex(ValueError, 'invalid numeric values'):
            load_mean_representation(path, SEQ_SHA, 4, sha256_file(path), (17,), MODEL, REVISION)

    def test_mean_interface_rejects_float32_overflow_and_nonpositive_layer_ids(self):
        overflow = self.save_archive('mean_overflow.npz', meta=self.meta('mean_5'), mean_layers=np.asarray([[1e+300]], dtype=np.float64))
        with self.assertRaisesRegex(ValueError, 'not finite after float32 conversion'):
            load_mean_representation(overflow, SEQ_SHA, 4, sha256_file(overflow), (5,), MODEL, REVISION)
        for layer in (0, -1):
            path = self.save_archive(f'mean_bad_layer_{layer}.npz', meta=self.meta(f'mean_{layer}'), mean_layers=np.ones((1, 1), dtype=np.float32))
            with self.subTest(layer=layer), self.assertRaisesRegex(ValueError, 'positive 1-based IDs'):
                load_mean_representation(path, SEQ_SHA, 4, sha256_file(path), (layer,), MODEL, REVISION)

    def test_feature_identity_requires_exact_integer_lengths(self):
        for field in ('len_full', 'len_used'):
            path = self.save_archive(f'bad_{field}.npz', meta=self.meta('mean_5', **{field: 3.5}), mean_layers=np.ones((1, 1), dtype=np.float32))
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'length/truncation mismatch'):
                load_mean_representation(path, SEQ_SHA, 3, sha256_file(path), (5,), MODEL, REVISION)
        path = self.save_archive('bad_expected_length.npz', meta=self.meta('mean_5', length=3), mean_layers=np.ones((1, 1), dtype=np.float32))
        with self.assertRaisesRegex(ValueError, 'expected sequence length must be a positive integer'):
            load_mean_representation(path, SEQ_SHA, 3.5, sha256_file(path), (5,), MODEL, REVISION)

    def test_residue_interface_uses_exact_one_based_coordinates_and_masks_unknown(self):
        path = self.save_archive('residue.npz', meta=self.meta('mean_5_11_17_23_29_33+resid_11_23_33|domidx'), resid_layers=np.arange(3 * 3 * 2, dtype=np.float16).reshape(3, 3, 2), resid_indices=np.asarray([1, 3, 4], dtype=np.int32))
        features, positions, dimension = load_residue_representation(path, SEQ_SHA, 4, sha256_file(path), RESIDUE_LAYERS, MODEL, REVISION)
        labels = [{'sequence_sha256': SEQ_SHA, 'position': '1', 'label': '1', 'mask': '1'}, {'sequence_sha256': SEQ_SHA, 'position': '2', 'label': 'NA', 'mask': '0'}, {'sequence_sha256': SEQ_SHA, 'position': '3', 'label': '0', 'mask': '1'}, {'sequence_sha256': SEQ_SHA, 'position': '4', 'label': '1', 'mask': '1'}]
        aligned = align_residue_labels(features, positions, 4, labels, SEQ_SHA)
        self.assertEqual(dimension, 2)
        np.testing.assert_array_equal(aligned['positions_1based'], [1, 3, 4])
        np.testing.assert_array_equal(aligned['labels'], [1, 0, 1])
        np.testing.assert_array_equal(aligned['mask'], [1, 1, 1])
        self.assertEqual(aligned['features'].shape, (3, 3, 2))
        self.assertEqual(aligned['n_unknown'], 1)

    def test_residue_interface_rejects_float32_overflow_and_nonpositive_layer_ids(self):
        overflow = self.save_archive('residue_overflow.npz', meta=self.meta('resid_11'), resid_layers=np.asarray([[[1e+300]]], dtype=np.float64), resid_indices=np.asarray([1], dtype=np.int32))
        with self.assertRaisesRegex(ValueError, 'not finite after float32 conversion'):
            load_residue_representation(overflow, SEQ_SHA, 4, sha256_file(overflow), (11,), MODEL, REVISION)
        for layer in (0, -1):
            path = self.save_archive(f'residue_bad_layer_{layer}.npz', meta=self.meta(f'resid_{layer}'), resid_layers=np.ones((1, 1, 1), dtype=np.float32), resid_indices=np.asarray([1], dtype=np.int32))
            with self.subTest(layer=layer), self.assertRaisesRegex(ValueError, 'positive 1-based IDs'):
                load_residue_representation(path, SEQ_SHA, 4, sha256_file(path), (layer,), MODEL, REVISION)

    def test_residue_interface_fails_closed_on_missing_labeled_coordinates(self):
        with self.assertRaisesRegex(ValueError, 'lack feature coordinates \\(1\\): 2'):
            align_residue_labels(np.zeros((3, 2, 2), dtype=np.float32), np.asarray([1, 3]), 3, [{'sequence_sha256': SEQ_SHA, 'position': '1', 'label': '0', 'mask': '1'}, {'sequence_sha256': SEQ_SHA, 'position': '2', 'label': '1', 'mask': '1'}], SEQ_SHA)

    def test_residue_interface_rejects_invalid_coordinates_and_mask_label_pairs(self):
        with self.assertRaisesRegex(ValueError, 'unique 1-based positions'):
            align_residue_labels(np.zeros((3, 2, 1)), np.asarray([2, 2]), 3, [], SEQ_SHA)
        with self.assertRaisesRegex(ValueError, 'integer 1-based positions'):
            align_residue_labels(np.zeros((3, 1, 1)), np.asarray([1.5]), 3, [], SEQ_SHA)
        with self.assertRaisesRegex(ValueError, 'masked residue has a biological class'):
            align_residue_labels(np.zeros((3, 1, 1)), np.asarray([1]), 1, [{'sequence_sha256': SEQ_SHA, 'position': '1', 'label': '0', 'mask': '0'}], SEQ_SHA)
        with self.assertRaisesRegex(ValueError, 'sequence hash mismatch'):
            align_residue_labels(np.zeros((3, 1, 1)), np.asarray([1]), 1, [{'sequence_sha256': 'b' * 64, 'position': '1', 'label': '1', 'mask': '1'}], SEQ_SHA)
        with self.assertRaisesRegex(ValueError, 'duplicate or out-of-range'):
            align_residue_labels(np.zeros((3, 1, 1)), np.asarray([1]), 1, [{'sequence_sha256': SEQ_SHA, 'position': '1', 'label': 'NA', 'mask': '0'}, {'sequence_sha256': SEQ_SHA, 'position': '1', 'label': '1', 'mask': '1'}], SEQ_SHA)
if __name__ == '__main__':
    unittest.main()
