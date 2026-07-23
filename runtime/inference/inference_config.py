from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Optional, Tuple


_CUDA_DEVICE = re.compile(r"cuda(?::[0-9]+)?\Z")


@dataclass(frozen=True)
class InferenceConfig:
    """Configuration encoded in the legacy comma-separated ``--arch`` value."""

    seed: int = 12345
    device: str = "cpu"
    threshold: float = 0.0
    backbone: Optional[str] = None
    embed_dim: Optional[int] = None
    base: Optional[int] = None
    blocks: Optional[Tuple[int, ...]] = None
    in_ch: int = 2


def _strip_arch_prefix(value: str) -> str:
    """Remove an optional ``name:`` prefix without breaking ``device=cuda:0``."""

    first_field = value.split(",", 1)[0]
    if ":" in first_field and "=" not in first_field.split(":", 1)[0]:
        return value.split(":", 1)[1]
    return value


def parse_arch(value: str) -> InferenceConfig:
    """
    Parse the backward-compatible ``--arch`` mini-language.

    Unknown fields and fields without ``=`` remain ignored, as they were in the
    original production entry point.
    """

    options = {
        "seed": 12345,
        "device": "cpu",
        "threshold": 0.0,
        "backbone": None,
        "embed_dim": None,
        "base": None,
        "blocks": None,
        "in_ch": 2,
    }
    if value:
        for field in _strip_arch_prefix(value).split(","):
            if not field or "=" not in field:
                continue
            key, raw = field.split("=", 1)
            key = key.strip().lower()
            raw = raw.strip()
            if key in ("seed", "s"):
                options["seed"] = int(raw)
            elif key in ("device", "dev"):
                options["device"] = raw.lower()
            elif key in ("thr", "threshold"):
                options["threshold"] = float(raw)
            elif key in ("backbone", "bb"):
                options["backbone"] = raw
            elif key in ("embed_dim", "embed", "d"):
                options["embed_dim"] = int(raw)
            elif key == "base":
                options["base"] = int(raw)
            elif key == "blocks":
                options["blocks"] = tuple(int(part) for part in raw.split("|") if part)
            elif key in ("in_ch", "inch", "channels"):
                options["in_ch"] = int(raw)

    config = InferenceConfig(**options)
    _validate_config(config)
    return config


def _validate_config(config: InferenceConfig) -> None:
    if not math.isfinite(config.threshold):
        raise ValueError("Inference threshold must be finite")
    if config.in_ch not in (2, 3, 4):
        raise ValueError(f"Unsupported model in_ch={config.in_ch}; expected 2, 3, or 4")
    if config.embed_dim is not None and config.embed_dim <= 0:
        raise ValueError("embed_dim must be positive")
    if config.base is not None and config.base <= 0:
        raise ValueError("base must be positive")
    if config.blocks is not None:
        if not config.blocks or any(count <= 0 for count in config.blocks):
            raise ValueError("blocks must contain positive counts")
    if not config.device:
        raise ValueError("device must not be empty")


def resolve_device(requested: str, *, cuda_available: bool) -> str:
    """
    Resolve a configured device without importing PyTorch.

    The established ``cuda``-to-CPU fallback is retained when CUDA is absent,
    while indexed devices such as ``cuda:0`` are now recognised correctly.
    """

    device = requested.strip().lower()
    if device == "cpu":
        return "cpu"
    if _CUDA_DEVICE.fullmatch(device):
        return device if cuda_available else "cpu"
    raise ValueError(
        f"Unsupported inference device {requested!r}; expected cpu, cuda, or cuda:<index>"
    )
