from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


def validate_config(config: dict) -> dict:
    if not isinstance(config, dict):
        raise ValueError("model configuration must be a JSON object")
    if type(config.get("schema_version")) is not int or config["schema_version"] != 1:
        raise ValueError("unsupported model configuration schema")
    if not isinstance(config.get("adapter"), str) or config["adapter"] not in {
        "esm2",
        "prott5",
    }:
        raise ValueError("unsupported model adapter")
    model_id = config.get("model_id")
    if (
        not isinstance(model_id, str)
        or not re.fullmatch(
            r"[A-Za-z0-9][A-Za-z0-9._-]*(?:/[A-Za-z0-9][A-Za-z0-9._-]*)?", model_id
        )
        or ".." in model_id
        or "--" in model_id
        or model_id.endswith((".git", ".", "-"))
        or len(model_id.rsplit("/", 1)[-1]) > 96
    ):
        raise ValueError("model_id must be a nonempty model repository name")
    revision = config.get("revision")
    if not isinstance(revision, str) or not re.fullmatch("[0-9a-f]{40}", revision):
        raise ValueError("model revision must be an immutable 40-character commit")
    layers = config.get("layers")
    if (
        not isinstance(layers, list)
        or not layers
        or any(type(x) is not int or x < 1 for x in layers)
        or len(set(layers)) != len(layers)
    ):
        raise ValueError("layers must be distinct positive transformer layer IDs")
    if not isinstance(config.get("pooling"), str) or config["pooling"] not in {
        "mean",
        "max",
    }:
        raise ValueError("pooling must be mean or max")
    long = config.get("long_sequence")
    if not isinstance(long, dict):
        raise ValueError("long_sequence must be a configuration object")
    limit, overlap = long.get("max_residues"), long.get("overlap")
    if (
        not isinstance(long.get("policy"), str)
        or long["policy"] not in {"window", "reject"}
        or type(limit) is not int
        or limit < 1
        or type(overlap) is not int
        or not 0 <= overlap < limit
    ):
        raise ValueError("invalid long-sequence window configuration")
    if not isinstance(config.get("dtype"), str) or config["dtype"] not in {
        "float32",
        "float16",
        "bfloat16",
    }:
        raise ValueError("unsupported inference dtype")
    if not isinstance(config.get("device"), str) or config["device"] not in {
        "cpu",
        "cuda",
        "mps",
    }:
        raise ValueError("device must be cpu, cuda or mps")
    if type(config.get("store_residues")) is not bool:
        raise ValueError("store_residues must be a boolean")
    return config


def load_config(path: Path) -> dict:
    return validate_config(json.loads(path.read_text(encoding="utf-8-sig")))


def download_snapshot(config: dict, cache_dir=None, local_only=False) -> Path:
    config = validate_config(config)
    try:
        from huggingface_hub import HfApi, snapshot_download
    except ImportError as error:
        raise ImportError(
            "Install conformeval[esm2] or conformeval[prott5] for model inference"
        ) from error
    common = [
        "config.json",
        "tokenizer_config.json",
        "special_tokens_map.json",
        "vocab.txt",
        "spiece.model",
    ]
    formats = ["model.safetensors", "pytorch_model.bin"]
    if not local_only:
        info = HfApi().model_info(config["model_id"], revision=config["revision"])
        names = {s.rfilename for s in info.siblings}
        formats = [name for name in formats if name in names][:1]
        if not formats:
            raise ValueError(
                "checkpoint must provide model.safetensors or pytorch_model.bin"
            )
    last_error = None
    for filename in formats:
        try:
            snapshot = Path(
                snapshot_download(
                    config["model_id"],
                    revision=config["revision"],
                    cache_dir=cache_dir,
                    local_files_only=local_only,
                    allow_patterns=common + [filename],
                    max_workers=2,
                )
            )
            if not (snapshot / filename).is_file():
                raise FileNotFoundError(f"checkpoint file is absent: {filename}")
            return snapshot
        except FileNotFoundError as error:
            if not local_only:
                raise
            last_error = error
    raise FileNotFoundError(
        "no complete cached checkpoint with model.safetensors or pytorch_model.bin"
    ) from last_error


def token_residue_map(ids, attention, special_mask, tokenizer, normalized_sequence):
    """Exclude every special/padded token and validate the residue reconstruction."""
    pieces = tokenizer.convert_ids_to_tokens(ids)
    mapping, reconstructed = [], ""
    for index, (piece, active, special) in enumerate(
        zip(pieces, attention, special_mask)
    ):
        if not active or special:
            continue
        residue = piece.replace("▁", "").replace("Ġ", "")
        if not re.fullmatch("[ACDEFGHIKLMNPQRSTVWYXUZOBJ]+", residue):
            raise ValueError(f"non-residue token in model output: {piece!r}")
        reconstructed += residue
        mapping.extend([index] * len(residue))
    if reconstructed != normalized_sequence or len(mapping) != len(normalized_sequence):
        raise ValueError("token/residue identity mismatch or tokenizer truncation")
    return np.asarray(mapping, dtype=np.int64)


def windows(length, limit, overlap, policy):
    if length <= limit:
        return [(0, length)]
    if policy == "reject":
        raise ValueError(f"sequence length {length} exceeds configured maximum {limit}")
    starts = list(range(0, length - limit + 1, limit - overlap))
    if starts[-1] != length - limit:
        starts.append(length - limit)
    return [(s, s + limit) for s in starts]


class HFAdapter:
    def __init__(self, config: dict, snapshot: Path):
        config = validate_config(config)
        if Path(snapshot).name != config["revision"]:
            raise ValueError(
                "snapshot directory must identify the pinned model revision"
            )
        try:
            import torch
            from transformers import EsmModel, EsmTokenizer, T5EncoderModel, T5Tokenizer
        except ImportError as error:
            raise ImportError(
                f"Install conformeval[{config['adapter']}] for this adapter"
            ) from error
        self.torch, self.config = torch, config
        self.device = config["device"]
        if self.device == "cuda" and not torch.cuda.is_available():
            raise ValueError("CUDA was requested but is unavailable")
        if self.device == "mps" and not torch.backends.mps.is_available():
            raise ValueError("MPS was requested but is unavailable")
        dtype = getattr(torch, config["dtype"])
        if config["adapter"] == "esm2":
            self.tokenizer = EsmTokenizer.from_pretrained(
                snapshot, local_files_only=True
            )
            self.model, loading = EsmModel.from_pretrained(
                snapshot,
                local_files_only=True,
                add_pooling_layer=False,
                torch_dtype=dtype,
                output_loading_info=True,
            )
        else:
            self.tokenizer = T5Tokenizer.from_pretrained(
                snapshot, do_lower_case=False, local_files_only=True
            )
            self.model, loading = T5EncoderModel.from_pretrained(
                snapshot,
                local_files_only=True,
                torch_dtype=dtype,
                output_loading_info=True,
            )
        if any(
            loading.get(k) for k in ["missing_keys", "mismatched_keys", "error_msgs"]
        ):
            raise ValueError(f"incomplete model checkpoint loading: {loading}")
        self.loading = {
            k: list(v) if isinstance(v, set) else v for k, v in loading.items()
        }
        self.model.to(self.device).eval()
        self.depth = int(
            getattr(self.model.config, "num_hidden_layers", None)
            or self.model.config.num_layers
        )
        self.dimension = int(
            getattr(self.model.config, "hidden_size", None) or self.model.config.d_model
        )
        if max(config["layers"]) > self.depth:
            raise ValueError(f"requested layer exceeds model depth {self.depth}")
        if config["adapter"] == "esm2":
            # ESM position IDs start above the padding ID; CLS/EOS also occupy positions.
            position_offset = int(self.model.config.pad_token_id) + 1
            special_tokens = self.tokenizer.num_special_tokens_to_add(pair=False)
            limit = (
                int(self.model.config.max_position_embeddings)
                - position_offset
                - special_tokens
            )
            if config["long_sequence"]["max_residues"] > limit:
                raise ValueError(
                    f"window exceeds the configured ESM residue limit {limit}"
                )

    def encode_window(self, sequence: str):
        normalized = (
            re.sub("[UZOBJ]", "X", sequence)
            if self.config["adapter"] == "prott5"
            else sequence.replace("J", "X")
        )
        text = (
            " ".join(normalized) if self.config["adapter"] == "prott5" else normalized
        )
        encoded = self.tokenizer(
            text,
            add_special_tokens=True,
            truncation=False,
            return_special_tokens_mask=True,
            return_tensors="pt",
        )
        mapping = token_residue_map(
            encoded["input_ids"][0].tolist(),
            encoded["attention_mask"][0].tolist(),
            encoded["special_tokens_mask"][0].tolist(),
            self.tokenizer,
            normalized,
        )
        inputs = {
            k: v.to(self.device)
            for k, v in encoded.items()
            if k in {"input_ids", "attention_mask"}
        }
        with self.torch.inference_mode():
            result = self.model(**inputs, output_hidden_states=True, return_dict=True)
        hidden = result.hidden_states
        if hidden is None or len(hidden) != self.depth + 1:
            raise ValueError("model returned an unexpected transformer layer layout")
        values = np.stack(
            [
                hidden[layer][0, mapping].detach().float().cpu().numpy()
                for layer in self.config["layers"]
            ]
        )
        if (
            values.shape != (len(self.config["layers"]), len(sequence), self.dimension)
            or not np.isfinite(values).all()
        ):
            raise ValueError("invalid residue representation shape or values")
        return values

    def encode(self, sequence: str):
        if not sequence or not re.fullmatch("[ACDEFGHIKLMNPQRSTVWYXUZOBJ]+", sequence):
            raise ValueError("sequence must contain uppercase protein residue letters")
        long = self.config["long_sequence"]
        segments = windows(
            len(sequence), long["max_residues"], long["overlap"], long["policy"]
        )
        values = np.zeros(
            (len(self.config["layers"]), len(sequence), self.dimension),
            dtype=np.float32,
        )
        coverage = np.zeros(len(sequence), dtype=np.int32)
        for start, end in segments:
            values[:, start:end] += self.encode_window(sequence[start:end])
            coverage[start:end] += 1
        if np.any(coverage < 1):
            raise ValueError("windowing left an unrepresented residue")
        values /= coverage[None, :, None]
        pooled = (
            values.mean(axis=1)
            if self.config["pooling"] == "mean"
            else values.max(axis=1)
        )
        return pooled, values, segments


def download_main():
    parser = argparse.ArgumentParser(
        description="Download pinned model weights and tokenizer for one adapter"
    )
    parser.add_argument("--model-config", required=True, type=Path)
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args()
    config = load_config(args.model_config)
    print(
        json.dumps(
            {
                "model_id": config["model_id"],
                "revision": config["revision"],
                "snapshot": str(download_snapshot(config, args.cache_dir)),
            }
        )
    )
    return 0
