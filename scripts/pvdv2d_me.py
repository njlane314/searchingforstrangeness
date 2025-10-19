#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pvdv2d_me.py — minimal PV→VD in 2D (U,V,W) with MinkowskiEngine.

CLI:
  --in <planes.bin>   float32 blob: [U(HxW), V(HxW), W(HxW)]
  --out <result.bin>  IAOK header + 1 float (global VD response)
  --metrics <meta.txt>
  --W <int> --H <int>
  --arch pvdv2d_me
  --beam2d-u bx by
  --beam2d-v bx by
  --beam2d-w bx by
"""

import argparse, struct, time
import numpy as np
import torch
import torch.nn as nn
import MinkowskiEngine as ME

# ---------- IO ----------
def read_chw_f32(path, H, W):
    n = 3*H*W
    a = np.fromfile(path, dtype=np.float32, count=n)
    if a.size != n: raise RuntimeError(f"{path}: expected {n} floats, got {a.size}")
    a = a.reshape(3, H, W)
    return a[0], a[1], a[2]

def write_result_IAOK_cls1(path, value):
    # IAOK v1 (little-endian), K=1, cls_bytes=4, no seg
    header = struct.pack("<4sIIIIIQQQ",
                         b"IAOK", 1, 1, 0, 0, 0,
                         np.uint64(4), np.uint64(0), np.uint64(0))
    with open(path, "wb") as f:
        f.write(header)
        f.write(np.float32(value).tobytes())

def write_metrics(meta_path, t_total_ms, t_setup_ms, t_infer_ms, t_post_ms,
                  features_path, feat_dim, max_rss_mb=0.0, cuda_mem_mb=0.0):
    with open(meta_path, "w") as m:
        m.write(f"t_total_ms={t_total_ms:.3f}\n")
        m.write(f"t_setup_ms={t_setup_ms:.3f}\n")
        m.write(f"t_infer_ms={t_infer_ms:.3f}\n")
        m.write(f"t_post_ms={t_post_ms:.3f}\n")
        m.write(f"max_rss_mb={max_rss_mb:.3f}\n")
        m.write(f"cuda_mem_mb={cuda_mem_mb:.3f}\n")
        m.write(f"features_path={features_path}\n")
        m.write(f"feat_dim={feat_dim}\n")

# ---------- fixed geometry ----------
DIRS8 = [( 1,0),(-1,0),(0,1),(0,-1),( 1,1),( 1,-1),(-1,1),(-1,-1)]
UNIT8 = torch.tensor([(dx,dy)/np.hypot(dx,dy) for (dx,dy) in DIRS8], dtype=torch.float32)

def lin_index_2d(dx, dy, k=3):
    r = k//2; return (dx + r) + k*(dy + r)

def set_fixed_kernel(conv, weight):
    w = weight.contiguous().float()
    if hasattr(conv, "kernel"): conv.kernel = torch.nn.Parameter(w, requires_grad=False)
    elif hasattr(conv, "weight"): conv.weight = torch.nn.Parameter(w, requires_grad=False)
    for p in conv.parameters(): p.requires_grad_(False)

# ---------- ME layers ----------
class OrientedBank2D(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv = ME.MinkowskiSubmanifoldConvolution(1,8,kernel_size=3,dimension=2,bias=False)
        w = torch.zeros(8,1,9)
        for k,(dx,dy) in enumerate(DIRS8): w[k,0,lin_index_2d(dx,dy)] = 1.0
        set_fixed_kernel(self.conv, w)
    def forward(self, st): return self.conv(st)

class MomentMix2D(nn.Module):
    def __init__(self):
        super().__init__()
        self.mix = ME.MinkowskiConvolution(8,3,kernel_size=1,dimension=2,bias=False)
        U = UNIT8; M = torch.stack([U[:,0]*U[:,0], U[:,1]*U[:,1], U[:,0]*U[:,1]], dim=0)
        w = torch.zeros(3,8,1); w[:,:,0] = M
        set_fixed_kernel(self.mix, w)
    def forward(self, r8): return self.mix(r8)

class DipoleMix2D(nn.Module):
    def __init__(self):
        super().__init__()
        self.mix = ME.MinkowskiConvolution(8,2,kernel_size=1,dimension=2,bias=False)
        w = torch.zeros(2,8,1); w[:,:,0] = UNIT8.T
        set_fixed_kernel(self.mix, w)
    def forward(self, r8): return self.mix(r8)

class CentralDiff2D(nn.Module):
    def __init__(self, axis: str):
        super().__init__()
        self.conv = ME.MinkowskiSubmanifoldConvolution(1,1,kernel_size=3,dimension=2,bias=False)
        w = torch.zeros(1,1,9)
        if axis == "x":
            w[0,0,lin_index_2d(+1,0)] = +0.5; w[0,0,lin_index_2d(-1,0)] = -0.5
        elif axis == "y":
            w[0,0,lin_index_2d(0,+1)] = +0.5; w[0,0,lin_index_2d(0,-1)] = -0.5
        else:
            raise ValueError("axis must be 'x' or 'y'")
        set_fixed_kernel(self.conv, w)
    def forward(self, st1): return self.conv(st1)

# ---------- PV→VD core ----------
class PV_DV_2D_ME(nn.Module):
    def __init__(self):
        super().__init__()
        self.bank = OrientedBank2D()
        self.moment = MomentMix2D()
        self.dipole = DipoleMix2D()
        self.dx = CentralDiff2D("x"); self.dy = CentralDiff2D("y")
        self.nms = ME.MinkowskiMaxPooling(kernel_size=3, stride=1, dimension=2)

    @staticmethod
    def _eig_invariants_2d(MF):
        Mxx, Myy, Mxy = MF[:,0], MF[:,1], MF[:,2]
        T = Mxx + Myy
        rad = torch.sqrt((Mxx-Myy)**2 + 4.0*(Mxy**2) + 1e-20)
        lam2 = 0.5 * (T - rad)
        C = (lam2 / (T + 1e-9)).clamp(min=0.0)
        return C, T.clamp_min(1e-9)

    def forward(self, st: ME.SparseTensor, beam2d_xy):
        device = st.F.device
        x = st.replace_feature((st.F > 0).float())
        r8 = self.bank(x); M3 = self.moment(r8); m2 = self.dipole(r8)
        dmx_dx = self.dx(x.replace_feature(m2.F[:,0:1])).F[:,0]
        dmy_dy = self.dy(x.replace_feature(m2.F[:,1:2])).F[:,0]
        div = dmx_dx + dmy_dy

        C, Tr = self._eig_invariants_2d(M3.F)
        S = torch.clamp(div, min=0.0) / (Tr + 1e-9)

        b = torch.tensor(beam2d_xy, dtype=torch.float32, device=device)
        b = b / (torch.linalg.norm(b) + 1e-9)
        m_norm = torch.linalg.norm(m2.F, dim=1) + 1e-9
        Fwd = torch.clamp((m2.F @ b), min=0.0) / m_norm

        proj = (UNIT8.to(device) @ b).view(8)
        fmask = (proj > 0).to(r8.F.dtype); bmask = (proj < 0).to(r8.F.dtype)
        Rplus  = (r8.F * fmask.unsqueeze(0)).sum(dim=1)
        Rminus = (r8.F * bmask.unsqueeze(0)).sum(dim=1)
        Bpv = Rplus / (Rplus + Rminus + 1e-9)

        lex = S + 1e-3*Bpv + 1e-6*C + 1e-9*Fwd
        pv_idx = torch.argmax(lex)

        coords = st.C[:,1:].float()             # [y,x]
        pv_yx = coords[pv_idx,:]
        dyx = coords - pv_yx[None,:]
        R_pix = torch.sqrt(torch.mean(torch.sum(dyx**2, dim=1))).item() + 1e-9

        disp = torch.clamp((dyx @ b) / R_pix, min=0.0)
        rnorm = torch.linalg.norm(dyx, dim=1) + 1e-9
        uhat = dyx / rnorm[:,None]
        Hf = torch.clamp((m2.F * uhat).sum(dim=1), min=0.0) / m_norm
        Q = 1.0 - (C / (C.max().clamp_min(1e-9)))
        DV = C * disp * Fwd * Hf * Q

        DV_st = x.replace_feature(DV.unsqueeze(1))
        DV_pooled = self.nms(DV_st).F[:,0]
        is_peak = (DV >= DV_pooled) & (DV > 0)
        peak_idx = torch.nonzero(is_peak, as_tuple=False).squeeze(1)

        peaks = []; top_score = 0.0
        if peak_idx.numel() > 0:
            order = torch.argsort(DV[peak_idx], descending=True)
            peak_idx = peak_idx[order]
            for j in range(min(10, peak_idx.numel())):
                i = peak_idx[j].item()
                y = int(coords[i,0].item()); xpix = int(coords[i,1].item())
                s = float(DV[i].item())
                peaks.append((xpix, y, s))
            top_score = float(DV[peak_idx[0]].item())

        return {
            "pv_rc": (int(pv_yx[0].item()), int(pv_yx[1].item())),
            "R_pix": float(R_pix),
            "peaks": peaks,
            "top": top_score,
            "coords": coords.cpu().numpy().astype(np.int32),
            "dv": DV.cpu().numpy().astype(np.float32),
            "pv": lex.cpu().numpy().astype(np.float32)
        }

def plane_to_sparse(img_hw, batch_id: int):
    H, W = img_hw.shape
    ys, xs = np.where(img_hw > 0.0)
    if ys.size == 0:
        C = np.array([[batch_id, 0, 0]], dtype=np.int32)
        F = np.zeros((1,1), dtype=np.float32)
        return C, F, True
    C = np.stack([np.full_like(ys, batch_id), ys.astype(np.int32), xs.astype(np.int32)], axis=1)
    F = np.ones((ys.size, 1), dtype=np.float32)
    return C.astype(np.int32), F, False

def run_plane(model, img_hw, beam2d_xy, batch_id):
    C_np, F_np, empty = plane_to_sparse(img_hw, batch_id)
    st = ME.SparseTensor(features=torch.from_numpy(F_np),
                         coordinates=torch.from_numpy(C_np),
                         device=torch.device("cpu"))
    if empty:
        return {"pv_rc": (0,0), "R_pix": 1.0, "peaks": [], "top": 0.0,
                "coords": np.zeros((0,2), np.int32), "dv": np.zeros((0,), np.float32),
                "pv": np.zeros((0,), np.float32)}
    with torch.no_grad():
        return model(st, beam2d_xy)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--H", type=int, required=True)
    ap.add_argument("--arch", required=True)
    ap.add_argument("--beam2d-u", type=float, nargs=2, required=True)
    ap.add_argument("--beam2d-v", type=float, nargs=2, required=True)
    ap.add_argument("--beam2d-w", type=float, nargs=2, required=True)
    args = ap.parse_args()

    if args.arch != "pvdv2d_me":
        raise SystemExit(f"Expected --arch pvdv2d_me (got {args.arch})")

    t0 = time.time()
    U, V, W = read_chw_f32(args.inp, args.H, args.W)
    t1 = time.time()

    model = PV_DV_2D_ME().eval()

    outU = run_plane(model, np.maximum(U,0).astype(np.float32), tuple(args.beam2d_u), 0)
    outV = run_plane(model, np.maximum(V,0).astype(np.float32), tuple(args.beam2d_v), 1)
    outW = run_plane(model, np.maximum(W,0).astype(np.float32), tuple(args.beam2d_w), 2)
    t2 = time.time()

    # Global response
    global_score = max(outU["top"], outV["top"], outW["top"])

    # Features sidecar (U,V,W; 34 floats per plane)
    def pack(o):
        pv_r, pv_c = o["pv_rc"]
        K = min(10, len(o["peaks"]))
        flat = [float(o["pv_rc"][1]), float(o["pv_rc"][0]), float(o["R_pix"]), float(K)]
        for i in range(K):
            c,r,s = o["peaks"][i]
            flat.extend([float(c), float(r), float(s)])
        # pad to 34 floats if fewer peaks
        while len(flat) < 34:
            flat.append(0.0)
        return flat
    feats = pack(outU) + pack(outV) + pack(outW)
    feat_path = args.out + ".feat.f32"
    np.asarray(feats, dtype=np.float32).tofile(feat_path)

    # IAOK (K=1, cls_bytes=4, no seg)
    write_result_IAOK_cls1(args.out, global_score)
    t3 = time.time()

    # Metrics include timings + sidecar pointer
    write_metrics(args.metrics,
                  t_total_ms=1000*(t3 - t0),
                  t_setup_ms=1000*(t1 - t0),
                  t_infer_ms=1000*(t2 - t1),
                  t_post_ms =1000*(t3 - t2),
                  features_path=feat_path,
                  feat_dim=len(feats))

if __name__ == "__main__":
    torch.set_grad_enabled(False)
    main()
