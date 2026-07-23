from __future__ import annotations

import io
import os
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Optional, Sequence, Tuple, Union

import numpy as np


REQUEST_MAGIC = b"IASP"
REQUEST_VERSION = 2
RESULT_MAGIC = b"IAOK"
RESULT_VERSION = 1

REQUEST_HEADER = struct.Struct("<4sIIIII")
PLANE_HEADER = struct.Struct("<Q")
RESULT_HEADER = struct.Struct("<4sIIIQ")

INT32_LE = np.dtype("<i4")
FLOAT32_LE = np.dtype("<f4")

PathLike = Union[str, os.PathLike]


class ProtocolError(RuntimeError):
    """A malformed, unsupported, or unreasonably large inference message."""


@dataclass(frozen=True)
class ProtocolLimits:
    """Resource ceilings checked before any payload-sized allocation."""

    max_file_bytes: int = 1 << 30
    max_planes: int = 16
    max_height: int = 16384
    max_width: int = 16384
    max_feature_dim: int = 16
    max_nnz_per_plane: int = 50_000_000
    max_result_values: int = 1_000_000


DEFAULT_LIMITS = ProtocolLimits()


@dataclass(frozen=True)
class SparsePlane:
    """One plane in wire order, with flat ``row,col`` coordinates and features."""

    coordinates: np.ndarray
    features: np.ndarray

    @property
    def nnz(self) -> int:
        return int(self.coordinates.size // 2)


@dataclass(frozen=True)
class SparseRequest:
    height: int
    width: int
    feature_dim: int
    planes: Tuple[SparsePlane, ...]


class _ExactReader:
    def __init__(self, stream: BinaryIO, total_bytes: int, source: str):
        self._stream = stream
        self._remaining = total_bytes
        self._source = source

    @property
    def remaining(self) -> int:
        return self._remaining

    def read(self, size: int, description: str) -> bytes:
        if size < 0 or size > self._remaining:
            raise ProtocolError(
                f"Short read for {description} in {self._source}: "
                f"need {size} bytes, have {self._remaining}"
            )
        data = self._stream.read(size)
        if len(data) != size:
            raise ProtocolError(
                f"Short read for {description} in {self._source}: "
                f"expected {size} bytes, got {len(data)}"
            )
        self._remaining -= size
        return data


def _checked_file_size(size: int, limits: ProtocolLimits, source: str) -> None:
    if size < REQUEST_HEADER.size:
        raise ProtocolError(f"Incomplete sparse request header in {source}")
    if size > limits.max_file_bytes:
        raise ProtocolError(
            f"Sparse request {source} is {size} bytes; limit is {limits.max_file_bytes}"
        )


def _validate_header(
    nplanes: int,
    height: int,
    width: int,
    feature_dim: int,
    *,
    expected_planes: Optional[int],
    expected_height: Optional[int],
    expected_width: Optional[int],
    limits: ProtocolLimits,
    source: str,
) -> None:
    if not 1 <= nplanes <= limits.max_planes:
        raise ProtocolError(
            f"Invalid sparse plane count {nplanes} in {source}; limit is {limits.max_planes}"
        )
    if not 1 <= height <= limits.max_height:
        raise ProtocolError(
            f"Invalid sparse height {height} in {source}; limit is {limits.max_height}"
        )
    if not 1 <= width <= limits.max_width:
        raise ProtocolError(
            f"Invalid sparse width {width} in {source}; limit is {limits.max_width}"
        )
    if not 1 <= feature_dim <= limits.max_feature_dim:
        raise ProtocolError(
            f"Invalid sparse feature dimension {feature_dim} in {source}; "
            f"limit is {limits.max_feature_dim}"
        )
    if expected_planes is not None and nplanes != expected_planes:
        raise ProtocolError(f"Expected {expected_planes} planes in {source}, got {nplanes}")
    if expected_height is not None and height != expected_height:
        raise ProtocolError(
            f"Request height {height} does not match CLI argument {expected_height}"
        )
    if expected_width is not None and width != expected_width:
        raise ProtocolError(
            f"Request width {width} does not match CLI argument {expected_width}"
        )


def _read_sparse_request(
    stream: BinaryIO,
    total_bytes: int,
    *,
    source: str,
    expected_planes: Optional[int],
    expected_height: Optional[int],
    expected_width: Optional[int],
    limits: ProtocolLimits,
) -> SparseRequest:
    _checked_file_size(total_bytes, limits, source)
    reader = _ExactReader(stream, total_bytes, source)

    magic, version, nplanes, height, width, feature_dim = REQUEST_HEADER.unpack(
        reader.read(REQUEST_HEADER.size, "sparse request header")
    )
    if magic != REQUEST_MAGIC:
        raise ProtocolError(f"Bad sparse request magic {magic!r} in {source}")
    if version != REQUEST_VERSION:
        raise ProtocolError(
            f"Unsupported sparse request version {version} in {source}; "
            f"expected {REQUEST_VERSION}"
        )
    _validate_header(
        nplanes,
        height,
        width,
        feature_dim,
        expected_planes=expected_planes,
        expected_height=expected_height,
        expected_width=expected_width,
        limits=limits,
        source=source,
    )

    max_sites = height * width
    planes = []
    for plane_index in range(nplanes):
        (nnz,) = PLANE_HEADER.unpack(
            reader.read(PLANE_HEADER.size, f"plane {plane_index} header")
        )
        if nnz > limits.max_nnz_per_plane:
            raise ProtocolError(
                f"Plane {plane_index} nnz={nnz} exceeds limit "
                f"{limits.max_nnz_per_plane} in {source}"
            )
        if nnz > max_sites:
            raise ProtocolError(
                f"Plane {plane_index} nnz={nnz} exceeds image capacity "
                f"{max_sites} in {source}"
            )

        coordinate_bytes = nnz * 2 * INT32_LE.itemsize
        feature_bytes = nnz * feature_dim * FLOAT32_LE.itemsize
        plane_bytes = coordinate_bytes + feature_bytes
        if plane_bytes > reader.remaining:
            raise ProtocolError(
                f"Plane {plane_index} declares {plane_bytes} payload bytes, "
                f"but only {reader.remaining} remain in {source}"
            )

        coordinates = np.frombuffer(
            reader.read(coordinate_bytes, f"plane {plane_index} coordinates"),
            dtype=INT32_LE,
        ).copy()
        features = np.frombuffer(
            reader.read(feature_bytes, f"plane {plane_index} features"),
            dtype=FLOAT32_LE,
        ).copy()

        if nnz:
            coordinate_pairs = coordinates.reshape(-1, 2)
            rows = coordinate_pairs[:, 0]
            columns = coordinate_pairs[:, 1]
            if np.any(rows < 0) or np.any(rows >= height):
                raise ProtocolError(
                    f"Plane {plane_index} contains out-of-range sparse rows in {source}"
                )
            if np.any(columns < 0) or np.any(columns >= width):
                raise ProtocolError(
                    f"Plane {plane_index} contains out-of-range sparse columns in {source}"
                )
            if not np.isfinite(features).all():
                raise ProtocolError(
                    f"Plane {plane_index} contains non-finite sparse features in {source}"
                )

        planes.append(SparsePlane(coordinates=coordinates, features=features))

    if reader.remaining:
        raise ProtocolError(
            f"Unexpected trailing data ({reader.remaining} bytes) in sparse request {source}"
        )

    return SparseRequest(
        height=height,
        width=width,
        feature_dim=feature_dim,
        planes=tuple(planes),
    )


def read_sparse_request(
    path: PathLike,
    *,
    expected_planes: Optional[int] = None,
    expected_height: Optional[int] = None,
    expected_width: Optional[int] = None,
    limits: ProtocolLimits = DEFAULT_LIMITS,
) -> SparseRequest:
    """Read and strictly validate an IASP v2 request from ``path``."""

    request_path = Path(path)
    with request_path.open("rb") as stream:
        total_bytes = os.fstat(stream.fileno()).st_size
        return _read_sparse_request(
            stream,
            total_bytes,
            source=str(request_path),
            expected_planes=expected_planes,
            expected_height=expected_height,
            expected_width=expected_width,
            limits=limits,
        )


def decode_sparse_request(
    payload: bytes,
    *,
    expected_planes: Optional[int] = None,
    expected_height: Optional[int] = None,
    expected_width: Optional[int] = None,
    limits: ProtocolLimits = DEFAULT_LIMITS,
) -> SparseRequest:
    """In-memory form used by golden tests and protocol tooling."""

    return _read_sparse_request(
        io.BytesIO(payload),
        len(payload),
        source="<bytes>",
        expected_planes=expected_planes,
        expected_height=expected_height,
        expected_width=expected_width,
        limits=limits,
    )


def encode_results(scores: Sequence[float], *, limits: ProtocolLimits = DEFAULT_LIMITS) -> bytes:
    """Encode raw class logits as an exact little-endian IAOK v1 response."""

    raw = np.asarray(scores).reshape(-1)
    if raw.size > limits.max_result_values:
        raise ProtocolError(
            f"Result contains {raw.size} values; limit is {limits.max_result_values}"
        )
    try:
        if not np.isfinite(raw).all():
            raise ProtocolError("Result contains non-finite class scores")
        with np.errstate(over="ignore", invalid="ignore"):
            values = np.asarray(raw, dtype=FLOAT32_LE)
    except (TypeError, ValueError, OverflowError) as error:
        raise ProtocolError("Result class scores must be finite numeric values") from error
    if not np.isfinite(values).all():
        raise ProtocolError("Result contains values outside the finite float32 range")

    count = int(values.size)
    header = RESULT_HEADER.pack(
        RESULT_MAGIC,
        RESULT_VERSION,
        count,
        0,
        count * FLOAT32_LE.itemsize,
    )
    return header + values.tobytes(order="C")


def write_results(path: PathLike, scores: Sequence[float]) -> None:
    """Write an IAOK response, truncating any prior file."""

    with Path(path).open("wb") as stream:
        stream.write(encode_results(scores))
