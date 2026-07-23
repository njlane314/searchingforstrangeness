from __future__ import annotations

import inspect
from collections.abc import Mapping
from typing import Any, Dict, Optional, Sequence, Set, Tuple

from model import SparseResNet2D, make_backbone
from fusion import MultiViewSetClassifier
from inference_config import InferenceConfig


def extract_state_dict(checkpoint: Any) -> Dict[str, Any]:
    """Extract a raw model state dict from supported training checkpoints."""

    if not isinstance(checkpoint, Mapping):
        raise RuntimeError(
            f"Unsupported weights object type: {type(checkpoint).__name__}"
        )
    if isinstance(checkpoint.get("model"), Mapping):
        return dict(checkpoint["model"])
    if isinstance(checkpoint.get("state_dict"), Mapping):
        return dict(checkpoint["state_dict"])

    import torch

    if any(torch.is_tensor(value) for value in checkpoint.values()):
        return dict(checkpoint)
    raise RuntimeError(
        "Weights file is a dict but does not contain a usable state_dict "
        "(expected keys: 'model' or 'state_dict', or a raw parameter dict)."
    )


def strip_module_prefix(state_dict: Mapping[str, Any]) -> Dict[str, Any]:
    """
    Remove DataParallel's ``module.`` prefix per key.

    Per-key handling avoids corrupting the unprefixed entries in a mixed state
    dict. A collision is rejected rather than silently choosing one tensor.
    """

    stripped: Dict[str, Any] = {}
    for key, value in state_dict.items():
        if not isinstance(key, str):
            raise RuntimeError("State-dict keys must be strings")
        output_key = key[len("module.") :] if key.startswith("module.") else key
        if output_key in stripped:
            raise RuntimeError(
                f"State-dict key collision while stripping module prefix: {output_key!r}"
            )
        stripped[output_key] = value
    return stripped


def _find_key_ending(state_dict: Mapping[str, Any], suffix: str) -> Optional[str]:
    return next((key for key in state_dict if key.endswith(suffix)), None)


def infer_hparams_from_state_dict(
    state_dict: Mapping[str, Any],
) -> Tuple[int, Tuple[int, ...], int, int]:
    """Infer ``base, blocks, embed_dim, num_views`` from checkpoint shapes."""

    plane_embedding_key = _find_key_ending(state_dict, "plane_emb.weight")
    if plane_embedding_key is None:
        raise RuntimeError(
            "Cannot infer embed_dim: missing 'plane_emb.weight' in state_dict."
        )
    plane_embedding = state_dict[plane_embedding_key]
    if len(plane_embedding.shape) != 2:
        raise RuntimeError(
            f"plane_emb.weight has unexpected shape {tuple(plane_embedding.shape)}"
        )
    num_views, embed_dim = map(int, plane_embedding.shape)

    base_key = _find_key_ending(state_dict, "backbone.stem.1.ln.weight")
    if base_key is None:
        raise RuntimeError(
            "Cannot infer base: missing 'backbone.stem.1.ln.weight' in state_dict."
        )
    base_shape = tuple(int(value) for value in state_dict[base_key].shape)
    if len(base_shape) != 1:
        raise RuntimeError(
            f"backbone.stem.1.ln.weight has unexpected shape {base_shape}"
        )
    base = base_shape[0]

    levels: Dict[int, Set[int]] = {}
    for key in state_dict:
        parts = key.split(".")
        for index, token in enumerate(parts):
            if token != "blocks" or index + 2 >= len(parts):
                continue
            try:
                level = int(parts[index + 1])
                block = int(parts[index + 2])
            except ValueError:
                continue
            levels.setdefault(level, set()).add(block)

    if not levels:
        raise RuntimeError(
            "Cannot infer blocks: no 'backbone.blocks.<li>.<bi>.*' keys found "
            "in state_dict."
        )
    ordered_levels = sorted(levels)
    if ordered_levels != list(range(len(ordered_levels))):
        raise RuntimeError(f"Checkpoint has non-contiguous block levels {ordered_levels}")
    for level in ordered_levels:
        indices = sorted(levels[level])
        if indices != list(range(len(indices))):
            raise RuntimeError(
                f"Checkpoint level {level} has non-contiguous block indices {indices}"
            )

    blocks = tuple(len(levels[level]) for level in ordered_levels)
    return base, blocks, embed_dim, num_views


def _torch_load_options(torch_load: Any) -> Dict[str, Any]:
    """Keep the trusted production checkpoint load stable across PyTorch APIs."""

    options: Dict[str, Any] = {"map_location": "cpu"}
    if "weights_only" in inspect.signature(torch_load).parameters:
        # The bundled training checkpoint also stores optimizer and NumPy RNG
        # state. It is a repository-controlled artifact, not a weights-only
        # archive, so modern PyTorch must be told to use its legacy loader.
        options["weights_only"] = False
    return options


def _load_checkpoint(path: str) -> Any:
    """Load the repository-controlled training checkpoint on CPU."""

    import torch

    return torch.load(path, **_torch_load_options(torch.load))


def build_llr_model(
    *,
    weights_path: str,
    device: str,
    config: InferenceConfig,
    plane_names: Sequence[str] = ("u", "v", "w"),
) -> Any:
    """Build the registered production architecture and strictly load its weights."""

    names = tuple(plane_names)
    if weights_path.startswith("random://"):
        if config.backbone is None or config.embed_dim is None:
            raise RuntimeError(
                "weights='random://' requires --arch to specify backbone and embed_dim "
                "(e.g. --arch "
                "'llr:backbone=base,embed_dim=128,thr=0.1,device=cpu')"
            )
        backbone = make_backbone(
            config.backbone,
            in_ch=config.in_ch,
            embed_dim=config.embed_dim,
        )
        model = MultiViewSetClassifier(
            backbone=backbone,
            embed_dim=config.embed_dim,
            plane_names=names,
        )
        model.to(device)
        model.eval()
        return model

    state_dict = strip_module_prefix(
        extract_state_dict(_load_checkpoint(weights_path))
    )
    if config.backbone is not None and config.embed_dim is not None:
        backbone = make_backbone(
            config.backbone,
            in_ch=config.in_ch,
            embed_dim=config.embed_dim,
        )
        model = MultiViewSetClassifier(
            backbone=backbone,
            embed_dim=config.embed_dim,
            plane_names=names,
        )
    else:
        base, blocks, embed_dim, num_views = infer_hparams_from_state_dict(state_dict)
        if num_views != len(names):
            raise RuntimeError(
                f"weights expect num_views={num_views}, but plane_names={names}"
            )

        if config.base is not None:
            base = config.base
        if config.blocks is not None:
            blocks = config.blocks
        if config.embed_dim is not None:
            embed_dim = config.embed_dim

        backbone = SparseResNet2D(
            in_ch=config.in_ch,
            base=base,
            blocks=blocks,
            embed_dim=embed_dim,
        )
        model = MultiViewSetClassifier(
            backbone=backbone,
            embed_dim=embed_dim,
            plane_names=names,
        )

    model.load_state_dict(state_dict, strict=True)
    model.to(device)
    model.eval()
    return model
