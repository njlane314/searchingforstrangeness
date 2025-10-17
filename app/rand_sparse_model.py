import numpy as np
import torch
import torch.nn as nn
import MinkowskiEngine as ME

class RandSparseCNN(nn.Module):
    """
    Submanifold-sparse CNN (2D):
      [ SubMConv3x3 -> ReLU ] x depth
      GlobalAvgPool -> [B, width]
    Random init (He) and frozen.
    """
    def __init__(self, in_ch: int = 3, width: int = 256, depth: int = 2, D: int = 2):
        super().__init__()
        assert D == 2
        blocks = []
        c = in_ch
        depth = max(1, int(depth))
        for _ in range(depth):
            blocks.append(ME.MinkowskiSubmanifoldConvolution(
                in_channels=c, out_channels=width, kernel_size=3,
                stride=1, dilation=1, bias=False, dimension=D))
            blocks.append(ME.MinkowskiReLU(inplace=True))
            c = width
        blocks.append(ME.MinkowskiGlobalAvgPooling())
        self.net = nn.Sequential(*blocks)
        for m in self.modules():
            if isinstance(m, ME.MinkowskiSubmanifoldConvolution):
                nn.init.kaiming_normal_(m.kernel, nonlinearity="relu")
                for p in m.parameters():
                    p.requires_grad_(False)
        for p in self.parameters():
            p.requires_grad_(False)

    @torch.no_grad()
    def forward(self, st: ME.SparseTensor) -> torch.Tensor:
        return self.net(st).F

def sparse_from_chw(chw: np.ndarray, thr: float = 0.0, device: str = "cpu") -> ME.SparseTensor:
    assert chw.ndim == 3
    C, H, W = chw.shape
    mask = np.any(np.abs(chw) > thr, axis=0)
    rows, cols = np.nonzero(mask)
    if rows.size == 0:
        rows = np.array([0], dtype=np.int32)
        cols = np.array([0], dtype=np.int32)
    coords = np.stack([np.zeros_like(rows, dtype=np.int32), rows.astype(np.int32), cols.astype(np.int32)], axis=1)
    feats  = np.stack([chw[c, rows, cols] for c in range(C)], axis=1).astype(np.float32, copy=False)
    return ME.SparseTensor(features=torch.from_numpy(feats),
                           coordinates=torch.from_numpy(coords),
                           device=device)

@torch.no_grad()
def extract_features_from_chw(chw: np.ndarray, width: int = 256, depth: int = 2,
                              thr: float = 0.0, seed: int = 12345, device: str = "cpu") -> np.ndarray:
    torch.manual_seed(seed)
    np.random.seed(seed)
    st = sparse_from_chw(chw, thr=thr, device=device)
    model = RandSparseCNN(in_ch=chw.shape[0], width=width, depth=depth).to(device).eval()
    feat = model(st)
    return feat[0].detach().cpu().numpy().astype(np.float32, copy=False)
