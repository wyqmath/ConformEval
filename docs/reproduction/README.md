# Reproduction environments

Install ConformEval using the [compatible dependency ranges](../installation.md).
For a specific archived experiment, consult its saved environment below.

| Record | Contents |
|---|---|
| [Historical 650M requirements](requirements-reference.txt) | Exact NumPy, SciPy, scikit-learn and PyTorch versions recorded for the archived reference environment; the file is preserved byte-for-byte from its former root location |
| [Measured resources](../resources.md) | Exact Python and dependency versions for encoder extraction and reader training, plus measured hardware, memory and timings |
| [Platform validation](../platforms.md) | Exact Python and numerical-library versions for tested macOS/Linux installations and their scoring results |
| [Platform machine-readable record](../../results/platform_validation/validation.json) | Actual environment versions, source hashes and test/replay evidence |
| [Adapter machine-readable record](../../results/adapter_validation/validation.json) | Checkpoint revisions, extraction environments and completed model runs |

Each result's `run.json` or `run_manifest.json` records its actual dependencies
alongside configuration and input/code hashes. Historical requirements describe
the original environment; the older PyTorch version in that file predates the
current model adapters. Install current inference extras using the installation
guide. Reference score replay reads saved predictions and needs only the base
package.
