"""Exercise the short example command and prevent accidental mixing of inputs."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ExampleCLITests(unittest.TestCase):
    def run_cli(self, *args, cwd=ROOT):
        env = dict(os.environ)
        env['PYTHONPATH'] = os.pathsep.join(filter(None, [str(ROOT / 'src'), env.get('PYTHONPATH')]))
        return subprocess.run(
            [sys.executable, '-m', 'conformeval.evaluate_predictions', *map(str, args)],
            cwd=cwd, env=env, capture_output=True, text=True, encoding='utf-8',
        )

    def test_example_matches_explicit_scoring_and_rejects_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            short = Path(temp) / 'short'
            explicit = Path(temp) / 'explicit'
            run = self.run_cli('--example', '--output-dir', short)
            self.assertEqual(run.returncode, 0, run.stderr)
            run = self.run_cli(
                '--task', 'KNOT-PRESENCE', '--data-dir', 'examples/base',
                '--predictions', 'examples/predictions/knot.tsv',
                '--model-metadata', 'examples/predictions/metadata.json',
                '--protocol', 'configs/scoring.json', '--output-dir', explicit,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            for name in ('metrics.json', 'scored_predictions.tsv.gz'):
                self.assertEqual((short / name).read_bytes(), (explicit / name).read_bytes())
            metrics = json.loads((short / 'metrics.json').read_text(encoding='utf-8'))['metrics']
            self.assertEqual(set(metrics), {'KNOT-PRESENCE'})
            self.assertEqual(metrics['KNOT-PRESENCE']['auroc'], 1.0)
            self.assertEqual(metrics['KNOT-PRESENCE']['n_samples'], 20)
            self.assertEqual(json.loads((short / 'run.json').read_text(encoding='utf-8'))['scope'], 'software_example')
            before = (short / 'run.json').read_bytes()
            run = self.run_cli('--example', '--output-dir', short)
            self.assertEqual(run.returncode, 2)
            self.assertIn('refusing to overwrite', run.stderr)
            self.assertEqual((short / 'run.json').read_bytes(), before)

    def test_example_does_not_mix_with_custom_inputs(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'output'
            for option, value in (
                ('--task', 'KNOT-PRESENCE'), ('--data-dir', 'data/base'),
                ('--predictions', 'custom.tsv'), ('--model-metadata', 'custom.json'),
                ('--protocol', 'custom.json'),
            ):
                with self.subTest(option=option):
                    run = self.run_cli('--example', option, value, '--output-dir', output)
                    self.assertEqual(run.returncode, 2)
                    self.assertIn('--example cannot be combined with ' + option, run.stderr)
                    self.assertFalse(output.exists())

    def test_missing_real_inputs_and_example_files_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp) / 'output'
            run = self.run_cli('--output-dir', output)
            self.assertEqual(run.returncode, 2)
            self.assertIn('required without --example', run.stderr)
            run = self.run_cli('--example', '--output-dir', output, cwd=temp)
            self.assertEqual(run.returncode, 2)
            self.assertIn('run it from the ConformEval repository root', run.stderr)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
