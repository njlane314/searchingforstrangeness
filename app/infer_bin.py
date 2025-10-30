#!/usr/bin/env python3
# Wrapper invoked by InferenceProduction::runInferenceDetailed
# Reads CHW float32 planes (U,V,W), extracts random sparse features with ME,
# writes a SIMPLE sidecar "<out>.feat.f32" (flat float32[D]) and a tiny
# result header with K=0 (or K=1 if you want a scalar score).

import argparse, struct, time
import numpy as np
import torch
from rand_sparse_model import extract_features_from_chw

def parse_arch(s: str):
    cfg = dict(w=256, d=2, thr=0.0, seed=12345, device='cpu')
    if not s:
        return cfg
    if ':' in s:
        s = s.split(':', 1)[1]
    for kv in s.split(','):
        if not kv or '=' not in kv: continue
        k, v = kv.split('=', 1)
        k = k.strip().lower(); v = v.strip()
        if k in ('w','width','c'): cfg['w'] = int(v)
        elif k in ('d','depth','l'): cfg['d'] = int(v)
        elif k in ('thr','th','t'): cfg['thr'] = float(v)
        elif k in ('seed','s'): cfg['seed'] = int(v)
        elif k in ('device','dev'): cfg['device'] = v
    return cfg

def read_chw_f32(path, C, H, W):
    arr = np.fromfile(path, dtype=np.float32)
    exp = C*H*W
    if arr.size != exp:
        raise RuntimeError(f"Expected {exp} floats in {path}, got {arr.size}")
    return arr.reshape(C, H, W)

def write_results_bin(path, cls_f32: np.ndarray):
    K = int(cls_f32.size)
    header = struct.pack('<4s5I3Q', b'IAOK', 1, K, 0, 0, 0, K*4, 0, 0)
    with open(path, 'wb') as f:
        f.write(header)
        if K:
            f.write(cls_f32.astype(np.float32, copy=False).tobytes(order='C'))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='in_path', required=True)
    ap.add_argument('--out', dest='out_path', required=True)
    ap.add_argument('--metrics', dest='metrics_path', required=True)
    ap.add_argument('--W', type=int, required=True)
    ap.add_argument('--H', type=int, required=True)
    ap.add_argument('--arch', default='')
    ap.add_argument('--weights', default='')
    args = ap.parse_args()

    t0 = time.time()
    cfg = parse_arch(args.arch)
    torch.manual_seed(cfg['seed']); np.random.seed(cfg['seed'])

    C,H,W = 3, args.H, args.W
    chw = read_chw_f32(args.in_path, C, H, W)
    device = 'cuda' if (cfg['device'].lower()=='cuda' and torch.cuda.is_available()) else 'cpu'

    t_setup = time.time()
    feat = extract_features_from_chw(chw, width=cfg['w'], depth=cfg['d'],
                                     thr=cfg['thr'], seed=cfg['seed'], device=device)
    t_infer = time.time()

    # Simple sidecar: flat float32
    features_path = args.out_path + ".feat.f32"
    with open(features_path, "wb") as f:
        f.write(feat.tobytes(order='C'))

    # Minimal cls: always emit a single scalar score (e.g., L2 norm)
    cls = np.array([float(np.linalg.norm(feat))], dtype=np.float32)
    write_results_bin(args.out_path, cls)

    t_write = time.time()
    # Metrics text parsed by C++ (plus pointers to sidecar)
    with open(args.metrics_path, "w") as m:
        m.write(f"t_total_ms={(t_write-t0)*1000.0:.3f}\n")
        m.write(f"t_setup_ms={(t_setup-t0)*1000.0:.3f}\n")
        m.write(f"t_infer_ms={(t_infer-t_setup)*1000.0:.3f}\n")
        m.write(f"t_post_ms={(t_write-t_infer)*1000.0:.3f}\n")
        m.write("max_rss_mb=0.0\n")
        m.write("cuda_mem_mb=0.0\n")
        m.write(f"features_path={features_path}\n")
        m.write(f"feat_dim={feat.size}\n")
        m.write(f"seed={cfg['seed']}\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
