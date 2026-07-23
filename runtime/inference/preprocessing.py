from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PreparedPlane:
    """Pure NumPy representation immediately before SparseTensor construction."""

    coordinates: np.ndarray
    features: np.ndarray
    available: float
    nnz: int


def _validate_feature_dimensions(payload_feature_dim: int, model_in_ch: int) -> None:
    if payload_feature_dim not in (2, 3):
        raise RuntimeError(
            f"Unsupported sparse payload feature_dim={payload_feature_dim}. Expected 2 "
            "([adc, nu_slice]) or 3 ([adc, nu_slice, dead])."
        )
    if model_in_ch not in (2, 3, 4):
        raise RuntimeError(
            f"Unsupported model in_ch={model_in_ch}; expected one of 2, 3, or 4"
        )


def validate_feature_schema(payload_feature_dim: int, model_in_ch: int) -> None:
    """Enforce the established production payload-to-model compatibility matrix."""

    _validate_feature_dimensions(payload_feature_dim, model_in_ch)
    compatible = payload_feature_dim == model_in_ch or (
        payload_feature_dim == 2 and model_in_ch in (3, 4)
    )
    if not compatible:
        raise RuntimeError(
            f"Sparse payload feature_dim={payload_feature_dim} is incompatible "
            f"with model in_ch={model_in_ch}. Only the 2->3/4 compatibility "
            "path is supported."
        )


def _coordinates_as_int32(coordinates: np.ndarray) -> np.ndarray:
    raw = np.asarray(coordinates)
    if raw.size == 0:
        return np.empty(0, dtype=np.int32)
    if np.issubdtype(raw.dtype, np.integer):
        as_int64 = raw.astype(np.int64, copy=False)
    else:
        try:
            numeric = raw.astype(np.float64)
        except (TypeError, ValueError, OverflowError) as error:
            raise RuntimeError("Sparse coordinates must be integers") from error
        if not np.isfinite(numeric).all() or not np.equal(numeric, np.floor(numeric)).all():
            raise RuntimeError("Sparse coordinates must be finite integers")
        as_int64 = numeric.astype(np.int64)
    info = np.iinfo(np.int32)
    if np.any(as_int64 < info.min) or np.any(as_int64 > info.max):
        raise RuntimeError("Sparse coordinates exceed the int32 range")
    return as_int64.astype(np.int32, copy=False).reshape(-1)


def prepare_sparse_plane(
    coordinates_flat: np.ndarray,
    features_flat: np.ndarray,
    *,
    payload_feature_dim: int,
    model_in_ch: int,
    height: int,
    width: int,
    threshold: float,
    batch_index: int = 0,
) -> PreparedPlane:
    """
    Apply the training-time sparse feature transform without PyTorch or ME.

    Model channel schemas:
      - C2: ``[occupancy, log1p(adc)]``
      - C3: ``[log1p(adc), nu_slice, dead]``
      - C4: ``[occupancy, log1p(adc), nu_slice, dead]``
    """

    _validate_feature_dimensions(payload_feature_dim, model_in_ch)
    if height <= 0 or width <= 0:
        raise RuntimeError(f"Invalid image dimensions {width}x{height}")
    if not math.isfinite(threshold):
        raise RuntimeError("Sparse ADC threshold must be finite")
    if not np.iinfo(np.int32).min <= batch_index <= np.iinfo(np.int32).max:
        raise RuntimeError("batch_index exceeds the int32 range")

    coordinates = _coordinates_as_int32(coordinates_flat)
    if coordinates.size % 2:
        raise RuntimeError("Sparse coords payload must contain row/col pairs")
    source_nnz = coordinates.size // 2

    try:
        raw_features = np.asarray(features_flat)
        if raw_features.size and not np.isfinite(raw_features).all():
            raise RuntimeError("Sparse features must be finite")
        features = np.asarray(raw_features, dtype=np.float32).reshape(-1)
    except (TypeError, ValueError, OverflowError) as error:
        raise RuntimeError("Sparse features must be finite numeric values") from error
    if features.size and not np.isfinite(features).all():
        raise RuntimeError("Sparse features exceed the finite float32 range")
    if features.size != source_nnz * payload_feature_dim:
        raise RuntimeError("Sparse coords and features have mismatched lengths")

    coordinate_pairs = coordinates.reshape(-1, 2)
    if source_nnz:
        rows = coordinate_pairs[:, 0]
        columns = coordinate_pairs[:, 1]
        if np.any(rows < 0) or np.any(rows >= height):
            raise RuntimeError("Sparse plane contains out-of-range rows")
        if np.any(columns < 0) or np.any(columns >= width):
            raise RuntimeError("Sparse plane contains out-of-range columns")

    features = features.reshape(-1, payload_feature_dim)
    adc = features[:, 0]
    nu_slice = features[:, 1]
    dead = (
        features[:, 2]
        if payload_feature_dim == 3
        else np.zeros_like(adc, dtype=np.float32)
    )

    keep = adc > threshold
    coordinate_pairs = coordinate_pairs[keep]
    adc = adc[keep]
    nu_slice = nu_slice[keep]
    dead = dead[keep]
    nnz = int(coordinate_pairs.shape[0])

    if nnz == 0:
        return PreparedPlane(
            coordinates=np.array([[batch_index, 0, 0]], dtype=np.int32),
            features=np.zeros((1, model_in_ch), dtype=np.float32),
            available=0.0,
            nnz=0,
        )

    sparse_coordinates = np.empty((nnz, 3), dtype=np.int32)
    sparse_coordinates[:, 0] = np.int32(batch_index)
    sparse_coordinates[:, 1:] = coordinate_pairs

    logq = adc.astype(np.float32, copy=True)
    np.maximum(logq, 0.0, out=logq)
    np.log1p(logq, out=logq)

    model_features = np.zeros((nnz, model_in_ch), dtype=np.float32)
    if model_in_ch == 2:
        model_features[:, 0] = 1.0
        model_features[:, 1] = logq
    elif model_in_ch == 3:
        model_features[:, 0] = logq
        model_features[:, 1] = nu_slice
        model_features[:, 2] = dead
    else:
        model_features[:, 0] = 1.0
        model_features[:, 1] = logq
        model_features[:, 2] = nu_slice
        model_features[:, 3] = dead

    return PreparedPlane(
        coordinates=sparse_coordinates,
        features=model_features,
        available=1.0,
        nnz=nnz,
    )
