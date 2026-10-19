from __future__ import annotations
import csv
import tempfile
import unittest
from pathlib import Path
from conformeval.evaluate_predictions import load_predictions, rank_positive

def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

class PredictionEvaluatorTests(unittest.TestCase):

    def test_rank_uses_fold_pool_and_average_ties(self):
        result = rank_positive(scores=__import__('numpy').asarray([0.9, 0.9, 0.1, 0.2]), labels=__import__('numpy').asarray([1, 0, 0, 0]))
        self.assertEqual(result['rank'], 1.5)
        self.assertEqual(result['universe'], 4)
        self.assertEqual(result['percentile'], 0.375)
        self.assertEqual(result['recall_at_1'], 0.5)

    def test_wrong_fold_and_missing_sample_are_rejected(self):
        eligible = {('KNOT-PRESENCE', 'a'): {'binding_group_id': 'g1'}, ('KNOT-PRESENCE', 'b'): {'binding_group_id': 'g2'}}
        group_to_fold = {('KNOT-PRESENCE', 'g1'): '0', ('KNOT-PRESENCE', 'g2'): '1'}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'predictions.tsv'
            write_tsv(path, [{'task_id': 'KNOT-PRESENCE', 'sample_id': 'a', 'outer_fold': '1', 'score': '0.8'}, {'task_id': 'KNOT-PRESENCE', 'sample_id': 'b', 'outer_fold': '1', 'score': '0.2'}], ['task_id', 'sample_id', 'outer_fold', 'score'])
            with self.assertRaisesRegex(ValueError, 'frozen outer assignment'):
                load_predictions(path, eligible, group_to_fold)
            write_tsv(path, [{'task_id': 'KNOT-PRESENCE', 'sample_id': 'a', 'outer_fold': '0', 'score': '0.8'}], ['task_id', 'sample_id', 'outer_fold', 'score'])
            with self.assertRaisesRegex(ValueError, 'coverage differs'):
                load_predictions(path, eligible, group_to_fold)

    def test_prediction_table_rejects_labels_and_extra_columns(self):
        eligible = {('KNOT-PRESENCE', 'a'): {'binding_group_id': 'g1'}}
        group_to_fold = {('KNOT-PRESENCE', 'g1'): '0'}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'predictions.tsv'
            write_tsv(path, [{'task_id': 'KNOT-PRESENCE', 'sample_id': 'a', 'outer_fold': '0', 'score': '0.8', 'label': '1'}], ['task_id', 'sample_id', 'outer_fold', 'score', 'label'])
            with self.assertRaisesRegex(ValueError, 'exactly these fields'):
                load_predictions(path, eligible, group_to_fold)

if __name__ == '__main__':
    unittest.main()
