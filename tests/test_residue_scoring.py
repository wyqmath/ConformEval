import csv, json, math, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from evaluate_residue_predictions import TASK, read_scores, calculate
from residue_scoring_core import validate_fold_roles
from verify_release import verify_delivery_file_set

class ResidueScorerTests(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'pred.tsv'
        self.samples = {'a': {'binding_group_id': 'g'}}
        self.expected = {'a': {1: 0, 3: 1}}
        self.folds = {'a': 0}
        self.rows = [dict(task_id=TASK, sample_id='a', position_1based=p, outer_fold=0, score=s) for p, s in [(1, 0.2), (3, 0.8)]]

    def read(self, rows):
        with self.path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(self.rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        return read_scores(self.path, self.samples, self.expected, self.folds, {'m'}, 'm')

    def test_complete_scores_and_metric(self):
        s = self.read(self.rows)
        rows, metrics = calculate(self.samples, self.expected, self.folds, s)
        self.assertEqual(metrics['m']['mean_per_protein_AP'], 1.0)
        self.assertEqual(metrics['group_bootstrap']['m']['macro_ap_ci95'], [1.0, 1.0])

    def test_unknown_missing_duplicate_nonfinite_wrong_fold(self):
        for changes, error in [({'position_1based': 2}, 'coordinate'), ({'outer_fold': 1}, 'frozen outer'), ({'score': float('nan')}, 'non-finite')]:
            bad = [dict(r) for r in self.rows]
            bad[0].update(changes)
            with self.subTest(changes=changes), self.assertRaisesRegex(ValueError, error):
                self.read(bad)
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            self.read(self.rows[:1])
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.read(self.rows + [self.rows[0]])

    def test_unregistered_delivery_file_is_rejected(self):
        root = Path(self.tmp.name)
        (root / 'manifest.json').write_text(json.dumps({'files': {}}))
        self.assertEqual(verify_delivery_file_set(root), 0)
        (root / 'unexpected.tsv').write_text('new data')
        with self.assertRaisesRegex(ValueError, 'unregistered'):
            verify_delivery_file_set(root)

    def test_exact_groups_cannot_cross_outer_roles(self):
        samples = {'a': {'binding_group_id': 'g'}, 'b': {'binding_group_id': 'g'}, 'c': {'binding_group_id': 'h'}}
        outer = [{'outer_fold': 0, 'sample_id': sid, 'binding_group_id': s['binding_group_id'], 'role': 'test' if sid == 'a' else 'train'} for sid, s in samples.items()]
        with self.assertRaisesRegex(ValueError, 'crosses outer'):
            validate_fold_roles(samples, outer, [], expected_outer=1, expected_inner=1)

    def test_outer_test_cannot_enter_inner_selection(self):
        samples = {'a': {'binding_group_id': 'g'}, 'b': {'binding_group_id': 'h'}}
        outer = [{'outer_fold': fold, 'sample_id': sid, 'binding_group_id': s['binding_group_id'], 'role': 'test' if (fold == 0) == (sid == 'a') else 'train'} for fold in [0, 1] for sid, s in samples.items()]
        inner = [{'outer_fold': 0, 'inner_fold': 0, 'sample_id': 'a', 'binding_group_id': 'g', 'role': 'validation'}]
        with self.assertRaisesRegex(ValueError, 'outer test sample entered'):
            validate_fold_roles(samples, outer, inner, expected_outer=2, expected_inner=1)
if __name__ == '__main__':
    unittest.main()
