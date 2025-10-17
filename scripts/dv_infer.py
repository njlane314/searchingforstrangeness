#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analytical displaced-vertex (DV) feature extractor (training-free, CPU-friendly).

Interface matches your InferenceEngine::runInferenceDetailed:
  --in <CHW float32>   : path to input planes, shape [3, H, W], row-major
  --out <results.bin>  : IAOK header file (K=0 unless --want-cls 1)
  --metrics <txt>      : metrics file; we write features_path & feat_dim here
  --W <int>, --H <int> : width and height
  --arch <str>         : config, e.g.
     "dv:thr=0.0,rs=2,rl=5,nms=3,topk=1,roi=63,
         angles=0|15|30|45|60|75|90|105|120|135|150|165,
         lengths=9|21,rel_thr=0.30,min_sep=20,
         rings=6|12|24|32,alpha=1.0,beta=0.0,globalM=0,
         device=cpu,seed=0"
  --weights <ignored>
  --want-cls 0/1       : if 1, we output ||feat||_2 as a single class score
  --want-seg  0        : ignored here
  --ref  <optional>    : text file with 3 lines "y x" (PV per plane). If absent,
                          PV displacement features fall back to 0 and are gated
                          by a cross-plane "pv_coverage" stat.

All operations are analytic: box averages (DoN), fixed thin-line banks for orientation,
non-maximum suppression, radial ring fractions, prong peak picking, entropy/anisotropy.
"""

import argparse
import math
import os
import struct
import time
from typing import List, Optional, Tuple, Union

import numpy as np
import torch
import torch.nn.functional as F


# ---------------------------
# Parsing, I/O, utilities
# ---------------------------

def parse_arch(s: str):
    # Sensible defaults for 512x512 and ~0.1% occupancy
    cfg = dict(
        thr=0.0,       # |ADC| threshold for occupancy
        rs=2,          # DoN small radius
        rl=5,          # DoN large radius (> rs)
        nms=3,         # NMS window
        topk=1,        # # candidates per plane (we assemble features for top-1)
        roi=63,        # ROI side (odd)
        angles="0|15|30|45|60|75|90|105|120|135|150|165",
        lengths="9|21",
        rel_thr=0.30,  # prong peak threshold (relative to max bin)
        min_sep=20.0,  # min prong angular separation (deg)
        rings="6|12|24|32",  # ring edges (px)
        alpha=1.0,     # weight for DoN term in vertexness
        beta=0.0,      # weight for multi-orientation surplus term in vertexness
        globalM=0,     # if 1, compute surplus M over whole image (slower); else compute at ROI only
        device="cpu",
        seed=0,
    )
    if s and ":" in s:
        s = s.split(":", 1)[1]
    if s:
        for kv in s.split(","):
            if "=" in kv:
                k, v = kv.split("=", 1)
                cfg[k.strip()] = v.strip()
    # casts
    cfg["thr"] = float(cfg["thr"])
    cfg["rs"] = int(cfg["rs"])
    cfg["rl"] = int(cfg["rl"])
    cfg["nms"] = int(cfg["nms"])
    cfg["topk"] = int(cfg["topk"])
    cfg["roi"] = int(cfg["roi"])
    cfg["angles"] = [float(x) for x in cfg["angles"].split("|") if x]
    cfg["lengths"] = [int(x) for x in cfg["lengths"].split("|") if x]
    cfg["rel_thr"] = float(cfg["rel_thr"])
    cfg["min_sep"] = float(cfg["min_sep"])
    cfg["rings"] = [int(x) for x in cfg["rings"].split("|") if x]
    cfg["alpha"] = float(cfg["alpha"])
    cfg["beta"] = float(cfg["beta"])
    cfg["globalM"] = int(cfg["globalM"]) != 0
    dev = cfg["device"].lower()
    cfg["device"] = "cuda" if (dev == "cuda" and torch.cuda.is_available()) else "cpu"
    cfg["seed"] = int(cfg["seed"])
    return cfg


def read_chw_f32(path, C, H, W):
    arr = np.fromfile(path, dtype=np.float32)
    exp = C * H * W
    if arr.size != exp:
        raise RuntimeError(f"Expected {exp} floats in {path}, got {arr.size}")
    return arr.reshape(C, H, W)


def read_ref_file(path: str) -> Optional[List[Tuple[float, float]]]:
    if not path or not os.path.exists(path):
        return None
    out = []
    with open(path, "r", encoding="utf-8") as f:
        for _, line in zip(range(3), f):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            y, x = float(parts[0]), float(parts[1])
            out.append((y, x))
    if not out:
        return None
    while len(out) < 3:
        out.append((float("nan"), float("nan")))
    return out


def write_results_bin(path, cls_vec: np.ndarray):
    K = int(cls_vec.size)
    header = struct.pack('<4s5I3Q', b'IAOK', 1, K, 0, 0, 0, K * 4, 0, 0)
    with open(path, "wb") as f:
        f.write(header)
        if K:
            f.write(cls_vec.astype(np.float32, copy=False).tobytes(order="C"))


# ---------------------------
# Core analytic primitives
# ---------------------------

def to_timg(x: np.ndarray, device: Union[str, torch.device]) -> torch.Tensor:
    """np.float32 [H,W] -> torch [1,1,H,W]"""
    t = torch.from_numpy(x).to(device)
    return t.view(1, 1, x.shape[0], x.shape[1])


def box_avg(timg: torch.Tensor, r: int) -> torch.Tensor:
    """Average in (2r+1)^2 window with same-padding. timg: [1,1,H,W]"""
    k = 2 * r + 1
    w = torch.ones((1, 1, k, k), dtype=timg.dtype, device=timg.device) / float(k * k)
    return F.conv2d(timg, w, padding=r)


def nms_mask(score: torch.Tensor, k: int) -> torch.Tensor:
    """score [1,1,H,W] -> bool mask [1,1,H,W] of local maxima in kxk windows."""
    pooled = F.max_pool2d(score, kernel_size=k, stride=1, padding=k // 2)
    return (score >= pooled) & (score > 0)


def candidate_from_map(score: torch.Tensor, nms: int, topk: int):
    mask = nms_mask(score, nms)
    vals = score[mask]
    if vals.numel() == 0:
        return []
    idx = torch.argsort(vals, descending=True)[:topk]
    coords = mask.nonzero(as_tuple=False)  # [N,4] b,c,y,x
    sel = coords[idx][:, 2:].cpu().numpy()
    v = vals[idx].cpu().numpy()
    return [(float(v[i]), int(sel[i, 0]), int(sel[i, 1])) for i in range(sel.shape[0])]


def extract_roi(img: np.ndarray, cy: int, cx: int, S: int) -> np.ndarray:
    """Crop SxS around (cy,cx); zero-pad outside."""
    H, W = img.shape
    h = S // 2
    y0, y1 = cy - h, cy + h + 1
    x0, x1 = cx - h, cx + h + 1
    out = np.zeros((S, S), np.float32)
    ys0, ys1 = max(0, y0), min(H, y1)
    xs0, xs1 = max(0, x0), min(W, x1)
    out[(ys0 - y0):(ys1 - y0), (xs0 - x0):(xs1 - x0)] = img[ys0:ys1, xs0:xs1]
    return out


def thin_line_bank(ks: int, angles_deg: List[float]) -> np.ndarray:
    """Return [n_ang, ks, ks] thin-line, zero-mean, L2-normalized kernels."""
    k = (ks - 1) // 2
    yy, xx = np.mgrid[-k:k + 1, -k:k + 1]
    bank = []
    for th in angles_deg:
        t = math.radians(th)
        dist = np.abs(-math.sin(t) * xx + math.cos(t) * yy)
        K = (dist <= 0.5).astype(np.float32)
        if K.sum() == 0:
            K[k, k] = 1.0
        K -= K.mean()
        K /= (np.linalg.norm(K) + 1e-6)
        bank.append(K)
    return np.stack(bank, 0)


def orientation_surplus_map(
    occ: np.ndarray,
    angles: List[float],
    lengths: List[int],
    device: str,
    compute_full: bool,
):
    """
    If compute_full=True: compute multi-orientation surplus M(x) over the whole image:
      M(x) = max_L [ sum_theta R_{theta,L}(x) - max_theta R_{theta,L}(x) ], with R >= 0.
    Else: return None and let caller compute surplus only inside ROI(s).
    """
    if not compute_full:
        return None
    timg = to_timg(occ, device)
    M = None
    for L in lengths:
        ks = L if (L % 2) else (L + 1)
        bank = thin_line_bank(ks, angles)    # [n_ang,ks,ks]
        Wt = torch.from_numpy(bank).to(device=device, dtype=timg.dtype).unsqueeze(1)
        pad = ks // 2
        resp = F.conv2d(timg, Wt, padding=pad)
        resp = resp.abs()
        sum_theta = resp.sum(dim=1, keepdim=True)
        max_theta, _ = resp.max(dim=1, keepdim=True)
        M_L = sum_theta - max_theta
        M = M_L if M is None else torch.maximum(M, M_L)
    Mmax = float(M.max().item()) if M is not None else 0.0
    if Mmax > 0:
        M = M / Mmax
    return M


def orientation_histogram_in_roi(
    roi_occ: np.ndarray,
    angles: List[float],
    lengths: List[int],
) -> np.ndarray:
    """
    Compute orientation histogram inside ROI using thin-line banks at given lengths.
    Sum magnitudes across lengths; apply mild radial weight so prongs away from center are emphasized.
    Returns h[angles] (float32).
    """
    S = roi_occ.shape[0]
    device = "cpu"
    img = to_timg(roi_occ, device)
    yy, xx = np.mgrid[:S, :S]
    cy, cx = S // 2, S // 2
    r = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2).astype(np.float32)
    rw = torch.from_numpy(r / (r.max() + 1e-6)).to(device)[None, None]
    hist = None
    for L in lengths:
        ks = L if (L % 2) else (L + 1)
        bank = thin_line_bank(ks, angles)
        Wt = torch.from_numpy(bank).to(device=device, dtype=img.dtype).unsqueeze(1)
        pad = ks // 2
        resp = F.conv2d(img, Wt, padding=pad).abs()
        val = (resp * rw).sum(dim=(2, 3))
        hist = val if hist is None else (hist + val)
    return hist[0].detach().cpu().numpy().astype(np.float32)


def prong_peaks(
    hist: np.ndarray,
    angles: List[float],
    rel_thr: float,
    min_sep_deg: float,
    topm: int = 3,
):
    """Greedy peak picking on a circular histogram (180° periodicity)."""
    if hist.size == 0:
        return 0, []
    ang = np.array(angles, np.float32)
    h = hist.copy()
    Hmax = float(h.max())
    if Hmax <= 0:
        return 0, []
    thr = rel_thr * Hmax
    picks = []
    while len(picks) < topm and h.max() >= thr:
        i = int(np.argmax(h))
        a = float(ang[i])
        picks.append(a)
        wrap = (np.abs(((ang - a + 90.0) % 180.0) - 90.0) <= (min_sep_deg / 2.0))
        h[wrap] = -np.inf
    return len(picks), picks


def ring_fractions(roi_occ: np.ndarray, edges: List[int]) -> List[float]:
    S = roi_occ.shape[0]
    cy = cx = S // 2
    yy, xx = np.mgrid[:S, :S]
    rr = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    tot = float(roi_occ.sum()) + 1e-6
    out = []
    prev = 0.0
    for e in edges:
        mask = (rr > prev) & (rr <= e)
        out.append(float(roi_occ[mask].sum()) / tot)
        prev = e
    return out


def norm_coords(y: int, x: int, H: int, W: int) -> Tuple[float, float]:
    return y / float(H), x / float(W)


def disp_to(y: float, x: float, y0: float, x0: float, H: int, W: int) -> float:
    """Normalized Euclidean distance between (y,x) and (y0,x0) by image diagonal."""
    if any([math.isnan(y), math.isnan(x), math.isnan(y0), math.isnan(x0)]):
        return 0.0
    dy = y - y0
    dx = x - x0
    return float(math.sqrt(dy * dy + dx * dx) / math.sqrt(H * H + W * W))


# ---------------------------
# DV features per plane
# ---------------------------

def dv_features_for_plane(
    ch: np.ndarray,
    H: int,
    W: int,
    cfg,
    pv_xy: Optional[Tuple[float, float]] = None,
    fullM: Optional[torch.Tensor] = None,
) -> np.ndarray:
    """
    Returns 24-dim feature vector for one plane.
    Layout:
      0   : V_peak          (vertexness at candidate)
      1   : DoN_peak        (mu_s - mu_l at candidate, >=0)
      2   : M_peak          (multi-orientation surplus at candidate; 0 if not computed)
      3-4 : y_norm, x_norm
      5   : disp_center_norm
      6   : disp_pv_norm    (0 if PV absent)
      7   : roi_occ_frac
      8-11: ring fractions (4)
      12  : orient_anisotropy (max - mean of orientation histogram)
      13  : orient_entropy     (Shannon / log(n_angles))
      14  : nprongs (0..3)
      15-20: cos/sin of top-3 prong angles (6)
      21-23: opening angles (3) normalized by 180°
    """
    occ = (np.abs(ch) > cfg["thr"]).astype(np.float32)

    dev = fullM.device if fullM is not None else torch.device(cfg["device"])
    t = to_timg(occ, dev)
    mu_s = box_avg(t, cfg["rs"])
    mu_l = box_avg(t, cfg["rl"])
    don = (mu_s - mu_l).clamp_min(0.0)
    mu_hat = don / (float(don.max().item()) + 1e-12)

    if fullM is not None and cfg["beta"] != 0.0:
        if fullM.device != don.device or fullM.dtype != don.dtype:
            M_hat = fullM.to(device=don.device, dtype=don.dtype)
        else:
            M_hat = fullM
    else:
        M_hat = torch.zeros_like(mu_hat)
    V = cfg["alpha"] * mu_hat + cfg["beta"] * M_hat
    peaks = candidate_from_map(V, cfg["nms"], topk=max(1, cfg["topk"]))
    if not peaks:
        return np.zeros(24, np.float32)

    Vpeak, cy, cx = peaks[0]
    don_peak = float(don[0, 0, cy, cx].item())
    Mpeak = float(M_hat[0, 0, cy, cx].item()) if M_hat is not None else 0.0

    S = cfg["roi"]
    roi_occ = extract_roi(occ, cy, cx, S)
    occ_frac = float(roi_occ.mean())

    hist = orientation_histogram_in_roi(roi_occ, cfg["angles"], cfg["lengths"])
    hsum = float(hist.sum()) + 1e-6
    hnorm = hist / hsum
    anis = float(hnorm.max() - hnorm.mean())
    ent = float(-(hnorm * np.log(hnorm + 1e-8)).sum() / math.log(hnorm.size + 1e-8))

    npr, pang = prong_peaks(hist, cfg["angles"], cfg["rel_thr"], cfg["min_sep"], topm=3)
    cs = []
    for i in range(3):
        if i < len(pang):
            a = math.radians(pang[i])
            cs += [math.cos(a), math.sin(a)]
        else:
            cs += [0.0, 0.0]
    opens = []
    if len(pang) >= 2:
        for i in range(len(pang)):
            for j in range(i + 1, len(pang)):
                d = abs(pang[i] - pang[j]) % 180.0
                d = d if d <= 90.0 else (180.0 - d)
                opens.append(d / 180.0)
    while len(opens) < 3:
        opens.append(0.0)

    yn, xn = norm_coords(cy, cx, H, W)
    disp_center = disp_to(cy, cx, H / 2.0, W / 2.0, H, W)
    disp_pv = 0.0
    if pv_xy is not None and all([not math.isnan(pv_xy[0]), not math.isnan(pv_xy[1])]):
        disp_pv = disp_to(cy, cx, pv_xy[0], pv_xy[1], H, W)

    vec = np.array(
        [
            float(Vpeak),
            float(don_peak),
            float(Mpeak),
            yn,
            xn,
            disp_center,
            float(disp_pv),
            occ_frac,
        ]
        + ring_fractions(roi_occ, cfg["rings"])[:4]
        + [anis, ent, float(npr)]
        + cs[:6]
        + opens[:3],
        dtype=np.float32,
    )
    return vec


# ---------------------------
# Per-event assembly
# ---------------------------

def dv_features_event(
    chw: np.ndarray,
    cfg,
    pv_list: Optional[List[Tuple[float, float]]] = None,
) -> np.ndarray:
    """
    Build final feature vector:
      per-plane (24 dims) * 3 = 72
      cross-plane (5 dims)
      Total = 77 dims.
    """
    C, H, W = chw.shape
    assert C == 3, "Expect 3 planes (U,V,W)."

    fullM = None
    if cfg["globalM"]:
        fullM = orientation_surplus_map(
            (np.abs(chw).max(axis=0) > cfg["thr"]).astype(np.float32),
            cfg["angles"],
            cfg["lengths"],
            cfg["device"],
            compute_full=True,
        )

    feats = []
    yxs = []
    disp_c_list = []
    disp_pv_list = []
    planes_with_candidate = 0
    planes_with_pv = 0

    for c in range(3):
        pv_xy = None if pv_list is None else pv_list[c]
        vec = dv_features_for_plane(chw[c], H, W, cfg, pv_xy=pv_xy, fullM=fullM)
        feats.append(vec)
        Vpeak = vec[0]
        if Vpeak > 0:
            planes_with_candidate += 1
            y_norm, x_norm = vec[3], vec[4]
            yxs.append((y_norm, x_norm))
            disp_c_list.append(vec[5])
            if pv_xy is not None and not any(map(math.isnan, pv_xy)):
                planes_with_pv += 1
                disp_pv_list.append(vec[6])

    if yxs:
        yn = np.array([p[0] for p in yxs], np.float32)
        xn = np.array([p[1] for p in yxs], np.float32)
        ystd = float(yn.std()) if yn.size > 1 else 0.0
        xstd = float(xn.std()) if xn.size > 1 else 0.0
        disp_c_mean = float(np.mean(disp_c_list)) if disp_c_list else 0.0
    else:
        ystd = xstd = disp_c_mean = 0.0

    disp_pv_mean = float(np.mean(disp_pv_list)) if disp_pv_list else 0.0
    coverage = planes_with_candidate / 3.0

    cross = np.array([disp_c_mean, disp_pv_mean, ystd, xstd, coverage], np.float32)

    return np.concatenate([feats[0], feats[1], feats[2], cross], axis=0)


# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--out", dest="out_path", required=True)
    ap.add_argument("--metrics", dest="metrics_path", required=True)
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--H", type=int, required=True)
    ap.add_argument(
        "--arch",
        default=(
            "dv:thr=0.0,rs=2,rl=5,nms=3,topk=1,roi=63,"
            "angles=0|15|30|45|60|75|90|105|120|135|150|165,"
            "lengths=9|21,rel_thr=0.30,min_sep=20,"
            "rings=6|12|24|32,alpha=1.0,beta=0.0,globalM=0,"
            "device=cpu,seed=0"
        ),
    )
    ap.add_argument("--weights", default="")
    ap.add_argument("--want-cls", type=int, default=0)
    ap.add_argument("--want-seg", type=int, default=0)
    ap.add_argument("--ref", default="", help='Optional path to PV pixels file (3 lines: "y x")')
    args = ap.parse_args()

    t0 = time.time()
    cfg = parse_arch(args.arch)
    torch.manual_seed(cfg["seed"])
    np.random.seed(cfg["seed"])

    C, H, W = 3, args.H, args.W
    chw = read_chw_f32(args.in_path, C, H, W)
    pv_list = read_ref_file(args.ref) if args.ref else None

    t_setup = time.time()
    feat = dv_features_event(chw, cfg, pv_list=pv_list)
    t_infer = time.time()

    sidecar = args.out_path + ".feat.f32"
    with open(sidecar, "wb") as f:
        f.write(feat.astype(np.float32).tobytes(order="C"))

    cls = (
        np.array([float(np.linalg.norm(feat))], np.float32)
        if args.want_cls
        else np.array([], np.float32)
    )
    write_results_bin(args.out_path, cls)
    t_write = time.time()

    with open(args.metrics_path, "w", encoding="utf-8") as m:
        m.write(f"t_total_ms={(t_write - t0) * 1000.0:.3f}\n")
        m.write(f"t_setup_ms={(t_setup - t0) * 1000.0:.3f}\n")
        m.write(f"t_infer_ms={(t_infer - t_setup) * 1000.0:.3f}\n")
        m.write(f"t_post_ms={(t_write - t_infer) * 1000.0:.3f}\n")
        m.write("max_rss_mb=0.0\ncuda_mem_mb=0.0\n")
        m.write(f"features_path={sidecar}\n")
        m.write(f"feat_dim={feat.size}\n")
        m.write(f"seed={cfg['seed']}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
