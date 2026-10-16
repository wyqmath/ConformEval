# Installation

Use Python 3.11 or newer. The base package installs only numerical scoring and reader dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install .
```

For Windows commands and actual system coverage, see the platform section below. For editable development use `python -m pip install -e .`. The wheel contains Python code; the cloned repository supplies versioned data and reference predictions through explicit paths.

## Model-specific optional dependencies

```bash
python -m pip install ".[esm2]"
# Or, for ProtT5:
python -m pip install ".[prott5]"
```

ESM-2 uses PyTorch, Transformers and huggingface-hub. ProtT5 additionally uses SentencePiece and protobuf. The base package does not install any of these inference dependencies. Both sequence and residue linear/ridge reader training use NumPy/scikit-learn.

Dependency ranges are declared in [`pyproject.toml`](../pyproject.toml). [`requirements.txt`](../requirements.txt) contains base ranges. The inference validation environment and archived reference environment are listed separately in [resources](resources.md); actual versions are saved with each run.

## Download pinned model files

```bash
confomeval-download --model-config configs/models/esm2_8m.json
confomeval-download --model-config configs/models/esm2_650m.json
confomeval-download --model-config configs/models/prott5.json
```

Run only the download(s) you need. The command downloads the declared immutable model revision, encoder weights, configuration and tokenizer into the Hugging Face cache. It prints the resolved snapshot path. Use `--cache-dir /path/to/model-cache` to choose the cache location. Extraction uses the same `--cache-dir`, or accepts `--snapshot /path/to/downloaded/snapshot`. Add `--local-files-only` to avoid network access during extraction. Model weights are external inputs and are not committed to this repository.

For inference choose `--device cpu`, `--device cuda` or `--device mps`; the requested device must exist. `--dtype float32`, `float16` and `bfloat16` override the configured precision. CPU float32 is the default. Select one precision and long-sequence policy before comparing models, and keep the run configuration.

## Verify the installation

```bash
python -m unittest discover -s tests -v
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

The tests run with base dependencies, without downloading models. Real-checkpoint inference and evaluation evidence is listed in [model validation](model_validation.md).

## Text files and operating systems

Use UTF-8 for JSON, FASTA and TSV inputs; an optional UTF-8 BOM is accepted. TSV and FASTA readers accept LF or CRLF line endings. New JSON/FASTA text and generated TSV records use LF. Invalid UTF-8 is rejected. Git checkout rules preserve frozen scientific/configuration files byte-for-byte, including any original CRLF, so their recorded hashes remain valid.

On Windows, the base dependency set includes `tzdata` for the named evaluation timezone. These PowerShell commands correspond to the base installation and scoring route:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install .
confomeval-replay-reference --benchmark-root . --output-dir runs/reference-replay
```

For `cmd.exe`, activate with `.venv\Scripts\activate.bat`. Windows execution is **未测试 / Not tested**. Actual macOS/Linux package-installation and scoring evidence, with Python versions, is recorded in the [platform matrix](platforms.md). The model inference device/precision coverage is separately recorded in [model validation](model_validation.md).
