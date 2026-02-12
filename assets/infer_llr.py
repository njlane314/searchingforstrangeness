#!/usr/bin/env python3

from __future__ import annotations

import argparse
import struct
import time
from typing import Dict, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import MinkowskiEngine as ME

from model import SparseResNet2D, make_backbone
from fusion import MultiViewSetClassifier


def parse_arch(s: str) -> Dict:
    """
    Backward-compatible parsing of --arch (comma-separated after optional 'prefix:').

    Existing keys (legacy inference):
      seed=12345
      device=cpu|cuda
      thr=0.0

    Optional keys for this model:
      backbone=tiny|small|base|wide
      embed_dim=128
      base=32
      blocks=2|2|2|2     (only if you want to override; usually inferred from weights)
    """
    cfg = dict(
        seed=12345,
        device="cpu",
        thr=0.0,
        backbone=None,
        embed_dim=None,
        base=None,
        blocks=None,
    )
    if not s:
        return cfg
    if ":" in s:
        s = s.split(":", 1)[1]
    for kv in s.split(","):
        if not kv or "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        k = k.strip().lower()
        v = v.strip()
        if k in ("seed", "s"):
            cfg["seed"] = int(v)
        elif k in ("device", "dev"):
            cfg["device"] = v
        elif k in ("thr", "threshold"):
            cfg["thr"] = float(v)
        elif k in ("backbone", "bb"):
            cfg["backbone"] = v
        elif k in ("embed_dim", "embed", "d"):
            cfg["embed_dim"] = int(v)
        elif k == "base":
            cfg["base"] = int(v)
        elif k == "blocks":
            # blocks=2|2|2|2
            parts = [p for p in v.split("|") if p]
            cfg["blocks"] = tuple(int(p) for p in parts)
    return cfg


def read_chw_f32(path: str, C: int, H: int, W: int) -> np.ndarray:
    arr = np.fromfile(path, dtype=np.float32)
    exp = C * H * W
    if arr.size != exp:
        raise RuntimeError(f"Expected {exp} floats in {path}, got {arr.size}")
    return arr.reshape(C, H, W)


def write_results_bin(path: str, cls_f32: np.ndarray) -> None:
    K = int(cls_f32.size)
    header = struct.pack("<4sIIIQ", b"IAOK", 1, K, 0, K * 4)
    with open(path, "wb") as f:
        f.write(header)
        if K:
            f.write(cls_f32.astype(np.float32, copy=False).tobytes(order="C"))


def _strip_module_prefix(sd: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
    if not any(k.startswith("module.") for k in sd.keys()):
        return sd
    return {k[len("module.") :]: v for k, v in sd.items()}


def _extract_state_dict(obj) -> Dict[str, torch.Tensor]:
    """
    Accept:
      - raw state_dict
      - checkpoint dict with 'model' or 'state_dict'
    """
    if not isinstance(obj, dict):
        raise RuntimeError(f"Unsupported weights object type: {type(obj).__name__}")
    if "model" in obj and isinstance(obj["model"], dict):
        return obj["model"]
    if "state_dict" in obj and isinstance(obj["state_dict"], dict):
        return obj["state_dict"]
    # If it looks like a raw state_dict already
    if any(isinstance(v, torch.Tensor) for v in obj.values()):
        return obj
    raise RuntimeError(
        "Weights file is a dict but does not contain a usable state_dict "
        "(expected keys: 'model' or 'state_dict', or a raw parameter dict)."
    )


def _find_key_ending(sd: Dict[str, torch.Tensor], suffix: str) -> Optional[str]:
    for k in sd.keys():
        if k.endswith(suffix):
            return k
    return None


def infer_hparams_from_state_dict(sd: Dict[str, torch.Tensor]) -> Tuple[int, Tuple[int, ...], int, int]:
    """
    Infer (base, blocks_tuple, embed_dim, num_views) from MultiViewSetClassifier state_dict.
    """
    k_pe = _find_key_ending(sd, "plane_emb.weight")
    if k_pe is None:
        raise RuntimeError("Cannot infer embed_dim: missing 'plane_emb.weight' in state_dict.")
    pe = sd[k_pe]
    if pe.ndim != 2:
        raise RuntimeError(f"plane_emb.weight has unexpected shape {tuple(pe.shape)}")
    num_views, embed_dim = int(pe.shape[0]), int(pe.shape[1])

    k_base = _find_key_ending(sd, "backbone.stem.1.ln.weight")
    if k_base is None:
        raise RuntimeError("Cannot infer base: missing 'backbone.stem.1.ln.weight' in state_dict.")
    base = int(sd[k_base].numel())

    levels: Dict[int, set] = {}
    for k in sd.keys():
        parts = k.split(".")
        for i, tok in enumerate(parts):
            if tok != "blocks":
                continue
            if i + 2 >= len(parts):
                continue
            try:
                li = int(parts[i + 1])
                bi = int(parts[i + 2])
            except Exception:
                continue
            levels.setdefault(li, set()).add(bi)

    if not levels:
        raise RuntimeError("Cannot infer blocks: no 'backbone.blocks.<li>.<bi>.*' keys found in state_dict.")

    blocks = tuple(len(levels[li]) for li in sorted(levels.keys()))
    return base, blocks, embed_dim, num_views


def build_llr_model(
    *,
    weights_path: str,
    device: str,
    arch_cfg: Dict,
    plane_names: Tuple[str, ...] = ("u", "v", "w"),
) -> nn.Module:
    """
    Create the model (matching training code) and load weights.

    Preference order:
      1) If --arch supplies backbone + embed_dim: use make_backbone() presets.
      2) Else: infer base/blocks/embed_dim from weights state_dict.
    """
    if weights_path.startswith("random://"):
        if arch_cfg.get("backbone") is None or arch_cfg.get("embed_dim") is None:
            raise RuntimeError(
                "weights='random://' requires --arch to specify backbone and embed_dim "
                "(e.g. --arch 'llr:backbone=base,embed_dim=128,thr=0.1,device=cpu')"
            )
        backbone = make_backbone(arch_cfg["backbone"], in_ch=2, embed_dim=int(arch_cfg["embed_dim"]))
        model = MultiViewSetClassifier(backbone=backbone, embed_dim=int(arch_cfg["embed_dim"]), plane_names=plane_names)
        model.to(device)
        model.eval()
        return model

    raw = torch.load(weights_path, map_location="cpu")
    sd = _strip_module_prefix(_extract_state_dict(raw))

    if arch_cfg.get("backbone") is not None and arch_cfg.get("embed_dim") is not None:
        backbone = make_backbone(arch_cfg["backbone"], in_ch=2, embed_dim=int(arch_cfg["embed_dim"]))
        model = MultiViewSetClassifier(backbone=backbone, embed_dim=int(arch_cfg["embed_dim"]), plane_names=plane_names)
    else:
        base, blocks, embed_dim, num_views = infer_hparams_from_state_dict(sd)
        if num_views != len(plane_names):
            raise RuntimeError(f"weights expect num_views={num_views}, but plane_names={plane_names}")

        # Optional overrides via --arch (for debugging)
        if arch_cfg.get("base") is not None:
            base = int(arch_cfg["base"])
        if arch_cfg.get("blocks") is not None:
            blocks = tuple(int(x) for x in arch_cfg["blocks"])
        if arch_cfg.get("embed_dim") is not None:
            embed_dim = int(arch_cfg["embed_dim"])

        backbone = SparseResNet2D(in_ch=2, base=int(base), blocks=tuple(blocks), embed_dim=int(embed_dim))
        model = MultiViewSetClassifier(backbone=backbone, embed_dim=int(embed_dim), plane_names=plane_names)

    model.load_state_dict(sd, strict=True)
    model.to(device)
    model.eval()
    return model


def plane_dense_to_sparse2d(
    plane_hw: np.ndarray,
    *,
    H: int,
    W: int,
    thr: float,
    device: str,
    batch_idx: int = 0,
) -> Tuple[ME.SparseTensor, float, int]:
    """
    Training-consistent sparsification:
      idx = plane > thr
      coords: [batch, y, x]
      feats : [occ=1, logq=log1p(max(adc,0))]
    Returns: (SparseTensor, available_flag, nnz)
    """
    flat = np.asarray(plane_hw, dtype=np.float32).reshape(-1)
    if flat.size != H * W:
        raise RuntimeError(f"Plane size {flat.size} != H*W ({H}*{W})")

    idx = np.flatnonzero(flat > thr)
    nnz = int(idx.size)

    if nnz == 0:
        coords = np.array([[batch_idx, 0, 0]], dtype=np.int32)
        feats = np.zeros((1, 2), dtype=np.float32)
        available = 0.0
    else:
        y = (idx // W).astype(np.int32, copy=False)
        x = (idx % W).astype(np.int32, copy=False)

        coords = np.empty((nnz, 3), dtype=np.int32)
        coords[:, 0] = np.int32(batch_idx)
        coords[:, 1] = y
        coords[:, 2] = x

        vals = flat[idx].astype(np.float32, copy=True)
        np.maximum(vals, 0.0, out=vals)
        np.log1p(vals, out=vals)

        feats = np.empty((nnz, 2), dtype=np.float32)
        feats[:, 0] = 1.0
        feats[:, 1] = vals

        available = 1.0

    coords_t = torch.from_numpy(coords).to(dtype=torch.int32, device="cpu")
    feats_t = torch.from_numpy(feats).to(dtype=torch.float32, device=device, non_blocking=True)
    st = ME.SparseTensor(features=feats_t, coordinates=coords_t, device=device)
    return st, available, nnz


def _sync_if_cuda(device: str) -> None:
    if str(device).startswith("cuda"):
        torch.cuda.synchronize()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", required=True)
    ap.add_argument("--metrics", dest="metrics_path", required=True)
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--H", type=int, required=True)
    ap.add_argument("--arch", default="")
    ap.add_argument("--weights", required=True)
    args = ap.parse_args()

    t0 = time.time()
    cfg = parse_arch(args.arch)

    torch.manual_seed(int(cfg["seed"]))
    np.random.seed(int(cfg["seed"]))

    req_dev = str(cfg["device"]).lower()
    if req_dev == "cuda" and torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"

    C, H, W = 3, int(args.H), int(args.W)
    chw = read_chw_f32(args.in_path, C, H, W)

    model = build_llr_model(weights_path=args.weights, device=device, arch_cfg=cfg, plane_names=("u", "v", "w"))

    thr = float(cfg["thr"])
    t_setup0 = time.time()

    inputs: Dict[str, ME.SparseTensor] = {}
    avail = []
    nnz = []
    for v, name in enumerate(("u", "v", "w")):
        st, a, n = plane_dense_to_sparse2d(chw[v], H=H, W=W, thr=thr, device=device)
        inputs[name] = st
        avail.append(a)
        nnz.append(n)

    available_mask = torch.tensor([avail], dtype=torch.float32, device=device)
    t_setup = time.time()

    print(
        f"[infer_llr] device={device} thr={thr} nnz(u,v,w)={tuple(nnz)} avail={tuple(int(x) for x in avail)}",
        flush=True,
    )

    _sync_if_cuda(device)
    with torch.no_grad():
        logits = model(inputs, available_mask=available_mask)  # [1,1]
    _sync_if_cuda(device)
    t_infer = time.time()

    logit = float(logits.reshape(-1)[0].detach().cpu().item())
    cls = np.array([logit], dtype=np.float32)
    write_results_bin(args.out_path, cls)

    t_write = time.time()

    # Best-effort RSS
    max_rss_mb = 0.0
    try:
        import resource

        # Linux: ru_maxrss is KB
        max_rss_mb = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0
    except Exception:
        max_rss_mb = 0.0

    with open(args.metrics_path, "w") as m:
        m.write(f"t_total_ms={(t_write - t0) * 1000.0:.3f}\n")
        m.write(f"t_setup_ms={(t_setup - t0) * 1000.0:.3f}\n")
        m.write(f"t_infer_ms={(t_infer - t_setup) * 1000.0:.3f}\n")
        m.write(f"t_post_ms={(t_write - t_infer) * 1000.0:.3f}\n")
        m.write(f"max_rss_mb={max_rss_mb:.3f}\n")
        m.write(f"seed={int(cfg['seed'])}\n")
        m.write(f"thr={thr}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
