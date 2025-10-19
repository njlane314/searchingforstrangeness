#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pvdv2d_me.py — Parameter-free PV→DV with MinkowskiEngine in pure 2D (U,V,W planes).

Inputs (wrapper contract):
  --in <planes.bin>     float32 blob: [U(HxW), V(HxW), W(HxW)] in this order (no header)
  --W <int>             width  (columns)
  --H <int>             height (rows)
  --arch pvdv2d_me
  --weights <ignored>   unused (no learning)
  --want-cls / --want-seg  supported flags; only cls produced here

Optional per-plane 2D beam directions (pixel coords; normalized internally):
  --beam2d-u bx by    default 0 1
  --beam2d-v bx by
  --beam2d-w bx by

Outputs:
  --out <result.bin>   IAOK header + K=1 classification float (max DV score across planes)
  --metrics <meta.txt> timings + (features_path, feat_dim)
  <out>.feat.f32       flat float32 features per plane:
                       [ pv_col, pv_row, R_pix, K,
                         dv0_col, dv0_row, dv0_s, ..., up to 10 triplets ] × 3 planes
"""

import argparse, struct, time
import numpy as np
import torch
import torch.nn as nn
import MinkowskiEngine as ME


# ---------------- I/O helpers ----------------
def read_chw_f32(path, H, W):
    n = 3 * H * W
    a = np.fromfile(path, dtype=np.float32, count=n)
    if a.size != n:
        raise RuntimeError(f"{path}: expected {n} floats, got {a.size}")
    a = a.reshape(3, H, W)
    return a[0], a[1], a[2]  # U, V, W

def write_result_IAOK(path, cls_floats, segW=0, segH=0, seg_bytes=b"", conf_bytes=b""):
    magic = b"IAOK"; version = 1; K = len(cls_floats)
    has_conf = 1 if conf_bytes else 0
    header = struct.pack("<4sIIIIQQQ", magic, version, K, segW, segH, has_conf,
                         np.uint64(4*K), np.uint64(len(seg_bytes)), np.uint64(len(conf_bytes)))
    with open(path, "wb") as f:
        f.write(header)
        if K:
            f.write(np.asarray(cls_floats, dtype=np.float32).tobytes())
        if seg_bytes:
            f.write(seg_bytes)
        if conf_bytes:
            f.write(conf_bytes)

def write_metrics(meta_path, t_total_ms, t_setup_ms, t_infer_ms, t_post_ms,
                  features_path="", feat_dim=0, seed=0, max_rss_mb=0.0, cuda_mem_mb=0.0):
    with open(meta_path, "w") as m:
        m.write(f"t_total_ms={t_total_ms:.3f}\n")
        m.write(f"t_setup_ms={t_setup_ms:.3f}\n")
        m.write(f"t_infer_ms={t_infer_ms:.3f}\n")
        m.write(f"t_post_ms={t_post_ms:.3f}\n")
        m.write(f"max_rss_mb={max_rss_mb:.3f}\n")
        m.write(f"cuda_mem_mb={cuda_mem_mb:.3f}\n")
        if features_path:
            m.write(f"features_path={features_path}\n")
            m.write(f"feat_dim={feat_dim}\n")
            m.write(f"seed={seed}\n")


# ---------------- 2D fixed geometry ----------------
# 8-neighbor unit directions in (dx,dy) with x=cols, y=rows
DIRS8 = [
    ( 1, 0), (-1, 0), (0, 1), (0,-1),
    ( 1, 1), ( 1,-1), (-1, 1), (-1,-1)
]
UNIT8 = torch.tensor(
    [(dx,dy)/np.hypot(dx,dy) for (dx,dy) in DIRS8],
    dtype=torch.float32
)  # (8,2)

def lin_index_2d(dx, dy, k=3):
    """Flatten (dx,dy) in {-1,0,1}^2 to kernel index with x fastest, then y."""
    r = k // 2
    ix = dx + r
    iy = dy + r
    return ix + k * iy

def set_fixed_kernel(conv, weight):
    """Assign fixed ME kernel/weight and freeze."""
    w = weight.contiguous().float()
    if hasattr(conv, "kernel"):
        conv.kernel = torch.nn.Parameter(w, requires_grad=False)
    elif hasattr(conv, "weight"):
        conv.weight = torch.nn.Parameter(w, requires_grad=False)
    for p in conv.parameters():
        p.requires_grad_(False)


# ---------------- ME layers (2D) ----------------
class OrientedBank2D(nn.Module):
    """Submanifold conv 1->8, each output channel = one neighbor offset (3x3)."""
    def __init__(self):
        super().__init__()
        self.conv = ME.MinkowskiSubmanifoldConvolution(
            in_channels=1, out_channels=8, kernel_size=3, dimension=2, bias=False
        )
        w = torch.zeros(8, 1, 9, dtype=torch.float32)
        for k,(dx,dy) in enumerate(DIRS8):
            w[k, 0, lin_index_2d(dx,dy,3)] = 1.0
        set_fixed_kernel(self.conv, w)
    def forward(self, st): return self.conv(st)

class MomentMix2D(nn.Module):
    """1x1 conv 8->3 for [Mxx, Myy, Mxy] = sum r_i u_i u_i^T"""
    def __init__(self):
        super().__init__()
        self.mix = ME.MinkowskiConvolution(8, 3, kernel_size=1, dimension=2, bias=False)
        U = UNIT8
        ux, uy = U[:,0], U[:,1]
        M = torch.stack([ux*ux, uy*uy, ux*uy], dim=0)  # (3,8)
        w = torch.zeros(3, 8, 1); w[:,:,0] = M
        set_fixed_kernel(self.mix, w)
    def forward(self, r8): return self.mix(r8)

class DipoleMix2D(nn.Module):
    """1x1 conv 8->2 for m = sum r_i u_i"""
    def __init__(self):
        super().__init__()
        self.mix = ME.MinkowskiConvolution(8, 2, kernel_size=1, dimension=2, bias=False)
        w = torch.zeros(2, 8, 1); w[:,:,0] = UNIT8.T
        set_fixed_kernel(self.mix, w)
    def forward(self, r8): return self.mix(r8)

class CentralDiff2D(nn.Module):
    """Central difference along one axis for 1-channel input (1->1) with 3x3 kernel."""
    def __init__(self, axis: str):
        super().__init__()
        self.conv = ME.MinkowskiSubmanifoldConvolution(1,1,kernel_size=3,dimension=2,bias=False)
        w = torch.zeros(1,1,9, dtype=torch.float32)
        if axis == "x":
            w[0,0, lin_index_2d(+1,0)] = +0.5
            w[0,0, lin_index_2d(-1,0)] = -0.5
        elif axis == "y":
            w[0,0, lin_index_2d(0,+1)] = +0.5
            w[0,0, lin_index_2d(0,-1)] = -0.5
        else:
            raise ValueError("axis must be 'x' or 'y'")
        set_fixed_kernel(self.conv, w)
    def forward(self, st1): return self.conv(st1)


# ---------------- PV→DV model (2D, parameter-free, per plane) ----------------
class PV_DV_2D_ME(nn.Module):
    def __init__(self):
        super().__init__()
        self.bank = OrientedBank2D()
        self.moment = MomentMix2D()
        self.dipole = DipoleMix2D()
        self.dx = CentralDiff2D("x")
        self.dy = CentralDiff2D("y")
        self.nms = ME.MinkowskiMaxPooling(kernel_size=3, stride=1, dimension=2)

    @staticmethod
    def _eig_invariants_2d(MF):
        """
        MF: features Nx3 = [Mxx, Myy, Mxy]
        Returns:
          C  = cornerness = lam2 / trace (∈[0,0.5])
          Tr = trace
        """
        Mxx, Myy, Mxy = MF[:,0], MF[:,1], MF[:,2]
        T = Mxx + Myy
        rad = torch.sqrt((Mxx - Myy)**2 + 4.0*(Mxy**2) + 1e-20)
        lam2 = 0.5 * (T - rad)
        C = (lam2 / (T + 1e-9)).clamp(min=0.0)
        return C, T.clamp_min(1e-9)

    def forward(self, st: ME.SparseTensor, beam2d_xy):
        """
        st: SparseTensor with (N,1) features (occupancy/charge) and coords [b,y,x]
        beam2d_xy: (bx,by) pixel-space vector (normalized inside)
        Returns dict:
          pv_rc, R_pix, peaks (list of (c,r,score)), top_score
        """
        device = st.F.device
        # binary occupancy for robustness
        x = st.replace_feature((st.F > 0).float())

        r8 = self.bank(x)       # (N,8)
        M3 = self.moment(r8)    # (N,3) -> Mxx,Myy,Mxy
        m2 = self.dipole(r8)    # (N,2) -> mx,my

        # divergence: central diffs on mx,my channels
        dmx_dx = self.dx(x.replace_feature(m2.F[:,0:1])).F[:,0]
        dmy_dy = self.dy(x.replace_feature(m2.F[:,1:2])).F[:,0]
        div = dmx_dx + dmy_dy                      # (N,)

        # invariants
        C, Tr = self._eig_invariants_2d(M3.F)      # (N,), (N,)
        S = torch.clamp(div, min=0.0) / (Tr + 1e-9)

        # beam forwardness (per site): (m·b)+ / ||m||
        b = torch.tensor(beam2d_xy, dtype=torch.float32, device=device)
        b = b / (torch.linalg.norm(b) + 1e-9)
        m_norm = torch.linalg.norm(m2.F, dim=1) + 1e-9
        Fwd = torch.clamp((m2.F @ b), min=0.0) / m_norm   # (N,)

        # PV backward-quietness from oriented responses
        proj = (UNIT8.to(device) @ b).view(8)
        fmask = (proj > 0).to(r8.F.dtype)
        bmask = (proj < 0).to(r8.F.dtype)
        Rplus  = (r8.F * fmask.unsqueeze(0)).sum(dim=1)
        Rminus = (r8.F * bmask.unsqueeze(0)).sum(dim=1)
        Bpv = Rplus / (Rplus + Rminus + 1e-9)             # (N,)

        # PV index (lexicographic via tiny weights)
        lex = S + 1e-3*Bpv + 1e-6*C + 1e-9*Fwd
        pv_idx = torch.argmax(lex)

        # coordinates (pixels)
        coords = st.C[:,1:].float() # (N,2) [y,x]
        pv_yx = coords[pv_idx, :]   # (y,x)

        # Event RMS radius in pixels over all active sites
        dyx = coords - pv_yx[None,:]
        R_pix = torch.sqrt(torch.mean(torch.sum(dyx**2, dim=1))).item() + 1e-9

        # Compose DV score on active sites
        # Displacement (downstream along b), normalized by R_pix
        disp = torch.clamp((dyx @ b) / R_pix, min=0.0)     # (N,)

        # Hb = forwardness to beam (already Fwd); Hf = projection onto unit(PV->x)
        rnorm = torch.linalg.norm(dyx, dim=1) + 1e-9
        uhat = dyx / rnorm[:,None]
        Hf = torch.clamp((m2.F * uhat).sum(dim=1), min=0.0) / m_norm

        # PV suppression via cornerness
        Q = 1.0 - (C / (C.max().clamp_min(1e-9)))

        DV = C * disp * Fwd * Hf * Q                       # (N,)

        # NMS: local maxima on the sparse 2D graph
        DV_st = x.replace_feature(DV.unsqueeze(1))         # (N,1)
        DV_pooled = self.nms(DV_st).F[:,0]
        is_peak = (DV >= DV_pooled) & (DV > 0)
        peak_idx = torch.nonzero(is_peak, as_tuple=False).squeeze(1)

        peaks = []
        top_score = 0.0
        if peak_idx.numel() > 0:
            order = torch.argsort(DV[peak_idx], descending=True)
            peak_idx = peak_idx[order]
            # keep up to 10 peaks
            for j in range(min(10, peak_idx.numel())):
                i = peak_idx[j].item()
                y = int(coords[i,0].item()); xpix = int(coords[i,1].item())
                s = float(DV[i].item())
                peaks.append((xpix, y, s))
            top_score = float(DV[peak_idx[0]].item())

        return {
            "pv_rc": (int(pv_yx[0].item()), int(pv_yx[1].item())),  # (row, col)
            "R_pix": float(R_pix),
            "peaks": peaks,
            "top": top_score
        }


# ---------------- per-plane runner ----------------
def plane_from_dense_to_sparse(img_hw, batch_id: int):
    """Convert HxW float32 to ME sparse coords [b,y,x] and occupancy features (1.0)."""
    H, W = img_hw.shape
    mask = img_hw > 0.0
    ys, xs = np.where(mask)
    if ys.size == 0:
        # create a single dummy inactive site to keep ME happy (no peaks will be found)
        C = np.array([[batch_id, 0, 0]], dtype=np.int32)
        F = np.zeros((1,1), dtype=np.float32)
        return C, F, True
    C = np.stack([np.full_like(ys, batch_id), ys.astype(np.int32), xs.astype(np.int32)], axis=1)
    F = np.ones((ys.size, 1), dtype=np.float32)  # occupancy; robust and scale-free
    return C.astype(np.int32), F, False


def run_plane(model, img_hw, beam2d_xy, batch_id):
    C_np, F_np, empty = plane_from_dense_to_sparse(img_hw, batch_id)
    st = ME.SparseTensor(
        features=torch.from_numpy(F_np),
        coordinates=torch.from_numpy(C_np),
        device=torch.device("cpu")
    )
    if empty:
        return {"pv_rc": (0,0), "R_pix": 1.0, "peaks": [], "top": 0.0}
    with torch.no_grad():
        out = model(st, beam2d_xy)
    return out


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--W", type=int, required=True)
    ap.add_argument("--H", type=int, required=True)
    ap.add_argument("--arch", required=True)
    ap.add_argument("--weights", required=False)  # ignored
    ap.add_argument("--want-cls", type=int, default=1)
    ap.add_argument("--want-seg", type=int, default=0)
    # optional per-plane 2D beam directions (pixels)
    ap.add_argument("--beam2d-u", type=float, nargs=2, default=[0.0, 1.0])
    ap.add_argument("--beam2d-v", type=float, nargs=2, default=[0.0, 1.0])
    ap.add_argument("--beam2d-w", type=float, nargs=2, default=[0.0, 1.0])
    args = ap.parse_args()

    if args.arch != "pvdv2d_me":
        raise SystemExit(f"Expected --arch pvdv2d_me (got {args.arch})")

    t0 = time.time()
    U, V, W = read_chw_f32(args.inp, args.H, args.W)
    t1 = time.time()

    model = PV_DV_2D_ME().eval()

    outU = run_plane(model, np.maximum(U, 0.0).astype(np.float32), tuple(args.beam2d_u), batch_id=0)
    outV = run_plane(model, np.maximum(V, 0.0).astype(np.float32), tuple(args.beam2d_v), batch_id=1)
    outW = run_plane(model, np.maximum(W, 0.0).astype(np.float32), tuple(args.beam2d_w), batch_id=2)
    t2 = time.time()

    # Classification: top DV score across planes
    top_cls = max(outU["top"], outV["top"], outW["top"]) if args.want_cls else 0.0

    # Features sidecar
    def pack(o):
        pv_r, pv_c = o["pv_rc"]
        K = min(10, len(o["peaks"]))
        flat = [float(pv_c), float(pv_r), float(o["R_pix"]), float(K)]
        for i in range(K):
            c,r,s = o["peaks"][i]
            flat.extend([float(c), float(r), float(s)])
        return flat
    feats = pack(outU) + pack(outV) + pack(outW)

    feat_path = args.out + ".feat.f32"
    np.asarray(feats, dtype=np.float32).tofile(feat_path)

    write_result_IAOK(args.out, [float(top_cls)])
    t3 = time.time()
    write_metrics(
        args.metrics,
        t_total_ms=1000*(t3 - t0),
        t_setup_ms=1000*(t1 - t0),
        t_infer_ms=1000*(t2 - t1),
        t_post_ms=1000*(t3 - t2),
        features_path=feat_path,
        feat_dim=len(feats),
        seed=0
    )

if __name__ == "__main__":
    torch.set_grad_enabled(False)
    main()
