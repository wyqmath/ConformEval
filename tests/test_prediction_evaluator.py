from __future__ import annotations
import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/base'
REFERENCE = ROOT / 'results/esm2_reference'
PROTOCOL = ROOT / 'configs/scoring_protocol.json'
from confomeval.evaluate_predictions import evaluate, load_predictions, rank_positive

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

    def test_evaluator_reproduces_existing_esm2_primary_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prediction_path = root / 'nested_predictions.tsv.gz'
            source_path = REFERENCE / 'predictions.tsv.gz'
            rows = []
            with gzip.open(source_path, 'rt', newline='', encoding='utf-8-sig') as source:
                for row in csv.DictReader(source, delimiter='\t'):
                    if row['arm'] == 'PLM-nested':
                        rows.append({key: row[key] for key in ('task_id', 'sample_id', 'outer_fold', 'score')})
            with gzip.open(prediction_path, 'wt', newline='', encoding='utf-8') as out:
                writer = csv.DictWriter(out, fieldnames=['task_id', 'sample_id', 'outer_fold', 'score'], delimiter='\t', lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
            provenance = json.loads((REFERENCE / 'provenance.json').read_text(encoding='utf-8-sig'))
            configuration = json.loads((REFERENCE / 'configuration.json').read_text(encoding='utf-8-sig'))
            metadata = {'model_id': provenance['encoder']['id'], 'model_revision': provenance['encoder']['revision'], 'prediction_method': 'nested balanced logistic regression', 'prediction_code_sha256': provenance['evaluator_sha256'], 'training_protocol': 'frozen_nested_grouped_oof_v1', 'candidate_layers': configuration['encoder']['selected_candidates'], 'pooling': 'mean_over_residues_without_special_tokens'}
            metadata_path = root / 'model_metadata.json'
            metadata_path.write_text(json.dumps(metadata), encoding='utf-8', newline='\n')
            result = evaluate(DATA, prediction_path, metadata_path, PROTOCOL, root / 'score')
            observed = json.loads((root / 'score/metrics.json').read_text(encoding='utf-8-sig'))['metrics']
            expected_knot = json.loads((REFERENCE / 'metrics_KNOT-PRESENCE.json').read_text(encoding='utf-8-sig'))
            expected_fs = json.loads((REFERENCE / 'metrics_FS-L2-grouped-PU-ranking.json').read_text(encoding='utf-8-sig'))
            self.assertAlmostEqual(observed['KNOT-PRESENCE']['auroc'], expected_knot['arms']['PLM-nested']['auroc'], places=12)
            self.assertAlmostEqual(observed['FS-L2-grouped-PU-ranking']['mean_held_case_rank_percentile'], expected_fs['arms']['PLM-nested']['mean_rank_percentile_lower_is_better'], places=12)
            self.assertEqual(result['model_training_protocol_self_attested'], True)
            self.assertEqual(result['scorer_uses_outer_test_labels_for_model_selection'], False)
if __name__ == '__main__':
    unittest.main()
