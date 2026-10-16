# Operating-system validation

ConformEval 0.5.1 was installed from the same pure-Python wheel on the systems below on 2026-10-17. Each installed package's 14 Python modules matched the final source bytes. The wheel SHA-256 is `1a5622c8b8908d5d9a85f4f23e038aac3493c7ca9c4a01dcb5d3d67fede3efd5`.

| Actual system | Python | Installation | Latest local suite | Complete archived reference scoring |
|---|---|---|---|---|
| macOS 26.4, arm64 | 3.14.5 | PASS, base-only venv | 55 passed, 0 skipped; 28.113 s | PASS; all three tasks, seven arms per task; 42.65 s |
| macOS 26.4, arm64 | 3.11.15 | PASS, base-only venv | 55 passed, 0 skipped; 25.510 s | Not repeated for this Python version; scoring is exercised by the test suite |
| Linux kernel 4.18.0-553.el8_10, x86_64, glibc 2.28 | 3.12.3 | PASS, fresh package in venv using existing Conda system-site libraries | 55 passed, 0 skipped; 40.343 s | PASS; all three tasks, seven arms per task; 78.36 s |
| Windows | 未测试 / Not tested | 未测试 | 未测试 | 未测试 |

The Linux checks ran on an actual remote Linux host. These rows describe base-package installation, training/scoring interface tests and archived prediction scoring. Model encoder inference was not rerun for this patch; its actual device/precision coverage is recorded separately in [model validation](model_validation.md). No Windows installation or execution is claimed.

[Machine-readable validation and logs](../results/platform_validation_v1/validation.json) record the actual Python/dependency versions, code hashes, test counts, replay summaries and source-log hashes. Local absolute paths are replaced by portable labels in exported evidence; original record/log hashes remain available. The latest full suites use the final text tests. Earlier macOS scoring records use the same final source modules; their original full-suite check remains in the source record.

## Actual dependencies

| Environment | NumPy | SciPy | scikit-learn |
|---|---|---|---|
| macOS Python 3.14 | 2.5.4 | 1.18.1 | 1.9.1 |
| macOS Python 3.11 | 2.4.6 | 1.17.1 | 1.9.1 |
| Linux Python 3.12 | 1.26.4 | 1.12.0 | 1.5.2 |

Linux used the preinstalled numerical libraries through a new venv's `--system-site-packages`; installing the wheel did not replace those libraries. The macOS base environments do not contain Torch, Transformers or SentencePiece. Third-party EncodingWarnings from the older Linux numerical stack are outside the project code; text tests reject warnings originating in ConformEval while also requiring both actual scoring CLIs to succeed with UTF-8 mode disabled.

## Text and checkout behavior

- JSON, FASTA and TSV inputs are read as UTF-8, accepting an optional initial BOM. Invalid UTF-8 raises an error.
- TSV readers accept LF/CRLF records and quoted multiline cells. FASTA readers accept LF/CRLF. New JSON/FASTA text and generated TSV record terminators use LF.
- Relative archive paths and output-hash keys use `/`, including when generated on Windows. The Windows runtime remains untested.
- `.gitattributes` preserves exact bytes for versioned data, results, configs, mini data, `.gz` and `.npz`. Source/document text checks out with LF. Tests exercise Git with `core.autocrlf=true` and `core.eol=crlf` without changing scientific files.
- Windows has an automatic `tzdata` dependency for the named evaluation timezone. The resource-usage module remains optional where unavailable.

These changes preserve the existing labels, masks, coordinates, group folds, predictions, parameters and numerical metrics. They do not make an untested model/task or operating system a completed evaluation.

## Local reproduction

Install the package with base dependencies and run from the repository root:

```bash
python -m pip install .
python -m unittest discover -s tests -v
confomeval-replay-reference --benchmark-root . --output-dir runs/platform-reference-replay
```

The output directory must be new. `REPLAY.json` should report `PASS`, with archived point metrics and group-bootstrap intervals matching at tolerance `1e-10`. [Installation](installation.md) supplies model-specific extras and Windows activation commands; [testing](testing.md) describes the scientific failure checks.
