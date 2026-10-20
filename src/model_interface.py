from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def _layer_ids(meta: dict, kind: str) -> tuple[int, ...]:
    if kind not in {"mean", "resid"}:
        raise ValueError(f"unsupported representation kind: {kind}")
    layer_set = meta.get("layer_set")
    if not isinstance(layer_set, str):
        raise ValueError("feature metadata has no layer_set")
    prefix = f"{kind}_"
    matches = []
    for part in layer_set.split("+"):
        part = part.split("|", 1)[0]
        if part.startswith(prefix):
            raw = part[len(prefix) :]
            try:
                ids = tuple((int(value) for value in raw.split("_")))
            except ValueError as exc:
                raise ValueError(f"invalid {kind} layer metadata: {layer_set}") from exc
            if any((layer < 1 for layer in ids)):
                raise ValueError(
                    f"{kind} layers must use positive 1-based IDs: {layer_set}"
                )
            matches.append(ids)
    if len(matches) != 1 or not matches[0] or len(set(matches[0])) != len(matches[0]):
        raise ValueError(f"ambiguous or invalid {kind} layer metadata: {layer_set}")
    return matches[0]


def _load_meta(archive: np.lib.npyio.NpzFile) -> dict:
    if "meta" not in archive.files:
        raise ValueError("feature archive has no metadata")
    try:
        meta = json.loads(str(archive["meta"].item()))
    except (ValueError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError("feature archive metadata is not valid JSON") from exc
    if not isinstance(meta, dict):
        raise ValueError("feature archive metadata must be a JSON object")
    return meta


def _validate_identity(
    meta: dict,
    expected_sequence_id: str,
    expected_length: int,
    expected_model: str,
    expected_revision: str,
) -> None:
    if type(expected_length) is not int or expected_length < 1:
        raise ValueError("expected sequence length must be a positive integer")
    if not expected_sequence_id:
        raise ValueError("expected sequence ID is required")
    if meta.get("sequence_id") != expected_sequence_id:
        raise ValueError("feature sequence ID mismatch")
    full_length = meta.get("len_full")
    used_length = meta.get("len_used")
    if (
        type(full_length) is not int
        or type(used_length) is not int
        or full_length != expected_length
        or (used_length != expected_length)
    ):
        raise ValueError("feature length/truncation mismatch")
    if meta.get("truncated") not in (False, 0, "false"):
        raise ValueError("truncated feature is not eligible")
    if meta.get("model") != expected_model or meta.get("revision") != expected_revision:
        raise ValueError("feature model revision mismatch")


def load_mean_representation(
    path: Path,
    expected_sequence_id: str,
    expected_length: int,
    selected_layers: tuple[int, ...],
    model: str,
    revision: str,
) -> tuple[np.ndarray, int]:
    with np.load(path, allow_pickle=False) as archive:
        key = "sequence_layers" if "sequence_layers" in archive.files else "mean_layers"
        if key not in archive.files:
            raise ValueError(f"feature archive has no mean_layers: {path.name}")
        meta = _load_meta(archive)
        values = np.asarray(archive[key])
    _validate_identity(meta, expected_sequence_id, expected_length, model, revision)
    layer_ids = _layer_ids(meta, "mean")
    if values.ndim != 2 or values.shape[0] != len(layer_ids) or values.shape[1] < 1:
        raise ValueError(
            f"invalid mean representation shape: {path.name} {values.shape}"
        )
    if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
        raise ValueError(f"mean representation has invalid numeric values: {path.name}")
    if not selected_layers or len(set(selected_layers)) != len(selected_layers):
        raise ValueError("selected mean layers must be nonempty and unique")
    missing = [layer for layer in selected_layers if layer not in layer_ids]
    if missing:
        raise ValueError(f"requested mean layers are absent: {missing}")
    with np.errstate(over="ignore", invalid="ignore"):
        selected = np.asarray(
            values[[layer_ids.index(layer) for layer in selected_layers]],
            dtype=np.float32,
        )
    if not np.isfinite(selected).all():
        raise ValueError(
            f"mean representation is not finite after float32 conversion: {path.name}"
        )
    return (selected, int(values.shape[1]))


def load_residue_representation(
    path: Path,
    expected_sequence_id: str,
    expected_length: int,
    selected_layers: tuple[int, ...],
    model: str,
    revision: str,
) -> tuple[np.ndarray, np.ndarray, int]:
    with np.load(path, allow_pickle=False) as archive:
        required = {"resid_layers", "resid_indices"}
        if not required.issubset(archive.files):
            raise ValueError(f"feature archive lacks residue arrays: {path.name}")
        meta = _load_meta(archive)
        values = np.asarray(archive["resid_layers"])
        raw_indices = np.asarray(archive["resid_indices"])
    _validate_identity(meta, expected_sequence_id, expected_length, model, revision)
    layer_ids = _layer_ids(meta, "resid")
    if (
        values.ndim != 3
        or values.shape[0] != len(layer_ids)
        or values.shape[1] < 1
        or (values.shape[2] < 1)
    ):
        raise ValueError(
            f"invalid residue representation shape: {path.name} {values.shape}"
        )
    if not np.issubdtype(values.dtype, np.number) or not np.isfinite(values).all():
        raise ValueError(
            f"residue representation has invalid numeric values: {path.name}"
        )
    if (
        raw_indices.ndim != 1
        or len(raw_indices) != values.shape[1]
        or (not np.issubdtype(raw_indices.dtype, np.integer))
    ):
        raise ValueError(f"residue coordinate array is invalid: {path.name}")
    indices = raw_indices.astype(np.int64, copy=False)
    if (
        len(np.unique(indices)) != len(indices)
        or np.any(indices < 1)
        or np.any(indices > expected_length)
    ):
        raise ValueError(
            f"residue coordinates must be unique 1-based positions within sequence: {path.name}"
        )
    if np.any(indices[1:] <= indices[:-1]):
        raise ValueError(
            f"residue coordinates must be strictly increasing: {path.name}"
        )
    if not selected_layers or len(set(selected_layers)) != len(selected_layers):
        raise ValueError("selected residue layers must be nonempty and unique")
    missing = [layer for layer in selected_layers if layer not in layer_ids]
    if missing:
        raise ValueError(f"requested residue layers are absent: {missing}")
    with np.errstate(over="ignore", invalid="ignore"):
        selected = np.asarray(
            values[[layer_ids.index(layer) for layer in selected_layers]],
            dtype=np.float32,
        )
    if not np.isfinite(selected).all():
        raise ValueError(
            f"residue representation is not finite after float32 conversion: {path.name}"
        )
    return (selected, indices, int(values.shape[2]))
