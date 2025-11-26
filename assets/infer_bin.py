#!/usr/bin/env python3

import argparse
import struct
import time

import numpy as np
import torch
import MinkowskiEngine as ME

from me_binary_model import build_model


def parse_arch(s: str):
    cfg = dict(seed=12345, device="cpu")
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
    return cfg


def read_chw_f32(path, C, H, W):
    arr = np.fromfile(path, dtype=np.float32)
    exp = C * H * W
    if arr.size != exp:
        raise RuntimeError(f"Expected {exp} floats in {path}, got {arr.size}")
    return arr.reshape(C, H, W)


def write_results_bin(path, cls_f32: np.ndarray):
    K = int(cls_f32.size)
    header = struct.pack("<4sIIIQ", b"IAOK", 1, K, 0, K * 4)
    with open(path, "wb") as f:
        f.write(header)
        if K:
            f.write(cls_f32.astype(np.float32, copy=False).tobytes(order="C"))


def chw_to_sparse(chw: np.ndarray, device: str) -> ME.SparseTensor:
    C, H, W = chw.shape
    if C != 3:
        raise RuntimeError(f"Expected 3 channels (U,V,W), got {C}")

    ys, xs = np.meshgrid(
        np.arange(H, dtype=np.int32),
        np.arange(W, dtype=np.int32),
        indexing="ij",
    )
    batch = np.zeros_like(ys, dtype=np.int32)
    coords = np.stack([batch, ys, xs], axis=-1).reshape(-1, 3)

    hw = H * W
    feats = chw.reshape(C, hw).T.astype(np.float32)
    if feats.shape[1] == 3:
        pad = np.zeros((feats.shape[0], 1), dtype=np.float32)
        feats = np.concatenate([feats, pad], axis=1)

    coords_t = torch.from_numpy(coords).int()
    feats_t = torch.from_numpy(feats).float()

    st = ME.SparseTensor(
        features=feats_t,
        coordinates=coords_t,
        device=device,
    )
    return st


def load_model(weights_path: str, device: str):
    model = build_model()

    if weights_path.startswith("random://"):
        print(
            f"[infer_bin] Using randomly initialised weights "
            f"(weights='{weights_path}')",
            flush=True,
        )
    else:
        state = torch.load(weights_path, map_location=device)
        if isinstance(state, dict) and "state_dict" in state:
            state = state["state_dict"]
        model.load_state_dict(state)

    model.to(device)
    model.eval()
    return model


def main():
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

    torch.manual_seed(cfg["seed"])
    np.random.seed(cfg["seed"])

    req_dev = cfg["device"].lower()
    if req_dev == "cuda" and torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"

    C, H, W = 3, args.H, args.W
    chw = read_chw_f32(args.in_path, C, H, W)

    model = load_model(args.weights, device=device)
    sparse = chw_to_sparse(chw, device=device)

    t_setup = time.time()

    with torch.no_grad():
        out = model(sparse)

    logits = out.F
    prob_pos = torch.sigmoid(logits)[0, 0].item()
    t_infer = time.time()

    cls = np.array([float(prob_pos)], dtype=np.float32)
    write_results_bin(args.out_path, cls)

    t_write = time.time()

    with open(args.metrics_path, "w") as m:
        m.write(f"t_total_ms={(t_write - t0) * 1000.0:.3f}\n")
        m.write(f"t_setup_ms={(t_setup - t0) * 1000.0:.3f}\n")
        m.write(f"t_infer_ms={(t_infer - t_setup) * 1000.0:.3f}\n")
        m.write(f"t_post_ms={(t_write - t_infer) * 1000.0:.3f}\n")
        m.write("max_rss_mb=0.0\n")
        m.write("feat_dim=1\n")
        m.write(f"seed={cfg['seed']}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
