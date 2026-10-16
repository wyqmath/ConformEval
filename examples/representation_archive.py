"""Write and validate a tiny synthetic mean-representation archive.

Run after installing the package. The vectors demonstrate the archive schema;
they are not model predictions or benchmark results.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from confomeval.model_interface import load_mean_representation, sha256_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    sequence = "ACDEFG"
    digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
    metadata = {
        "seq_sha256": digest, "len_full": len(sequence), "len_used": len(sequence),
        "truncated": False, "model": "synthetic-schema-example", "revision": "v1",
        "layer_set": "mean_17_33",
    }
    path = output / "example.npz"
    # Save one vector per declared layer, without special-token positions.
    np.savez_compressed(path, meta=np.array(json.dumps(metadata)),
                        mean_layers=np.arange(8, dtype=np.float32).reshape(2, 4))
    values, dimension = load_mean_representation(
        path, digest, len(sequence), sha256_file(path), (17, 33),
        metadata["model"], metadata["revision"],
    )
    print(json.dumps({"archive": str(path), "shape": list(values.shape),
                      "dimension": dimension, "sequence_sha256": digest}))


if __name__ == "__main__":
    main()
