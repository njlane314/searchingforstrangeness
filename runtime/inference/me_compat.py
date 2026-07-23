from __future__ import annotations

from typing import TYPE_CHECKING

import torch
import torch.nn as nn

import MinkowskiEngine as ME

if TYPE_CHECKING:
    from preprocessing import PreparedPlane


def replace_feature(x: ME.SparseTensor, new_features: torch.Tensor) -> ME.SparseTensor:
    """Replace sparse features across both current and legacy ME APIs."""

    if hasattr(x, "replace_feature"):
        return x.replace_feature(new_features)

    coordinate_map_key = getattr(x, "coordinate_map_key", None)
    coordinate_manager = getattr(x, "coordinate_manager", None)
    legacy_key = getattr(x, "coords_key", None)
    legacy_manager = getattr(x, "coords_man", None)
    if legacy_manager is None:
        legacy_manager = getattr(x, "coords_manager", None)

    if coordinate_map_key is not None and coordinate_manager is not None:
        try:
            return ME.SparseTensor(
                features=new_features,
                coordinate_map_key=coordinate_map_key,
                coordinate_manager=coordinate_manager,
                device=new_features.device,
            )
        except TypeError:
            pass

    if legacy_key is not None and legacy_manager is not None:
        try:
            return ME.SparseTensor(
                features=new_features,
                coords_key=legacy_key,
                coords_manager=legacy_manager,
                device=new_features.device,
            )
        except TypeError:
            pass

    return ME.SparseTensor(
        features=new_features,
        coordinates=x.C,
        device=new_features.device,
    )


def submanifold_convolution(
    in_ch: int,
    out_ch: int,
    *,
    kernel_size: int,
    dimension: int,
    bias: bool = False,
) -> nn.Module:
    """Construct a submanifold convolution across supported ME releases."""

    if hasattr(ME, "MinkowskiSubmanifoldConvolution"):
        return ME.MinkowskiSubmanifoldConvolution(
            int(in_ch),
            int(out_ch),
            kernel_size=kernel_size,
            dimension=dimension,
            bias=bias,
        )
    try:
        return ME.MinkowskiConvolution(
            int(in_ch),
            int(out_ch),
            kernel_size=kernel_size,
            stride=1,
            dimension=dimension,
            bias=bias,
            expand_coordinates=False,
        )
    except TypeError:
        return ME.MinkowskiConvolution(
            int(in_ch),
            int(out_ch),
            kernel_size=kernel_size,
            stride=1,
            dimension=dimension,
            bias=bias,
        )


def make_sparse_tensor(plane: "PreparedPlane", *, device: str) -> ME.SparseTensor:
    """Cross the pure-NumPy preprocessing boundary into ME."""

    coordinates = torch.from_numpy(plane.coordinates).to(dtype=torch.int32, device="cpu")
    features = torch.from_numpy(plane.features).to(
        dtype=torch.float32,
        device=device,
        non_blocking=True,
    )
    return ME.SparseTensor(features=features, coordinates=coordinates, device=device)
