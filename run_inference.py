import sys
import os
try:
  import MinkowskiEngine as ME
except ModuleNotFoundError:
  print("Error: MinkowskiEngine is required but not installed. Please install MinkowskiEngine to run this script.", file=sys.stderr)
  sys.exit(1)
import torch
import numpy as np
import argparse
import time
import torch.nn as nn

class ResidualBlock(nn.Module):
  def __init__(self, in_channels, out_channels, dimension):
    super().__init__()
    self.conv1 = ME.MinkowskiConvolution(in_channels, out_channels, kernel_size=3, dimension=dimension)
    self.bn1 = ME.MinkowskiBatchNorm(out_channels)
    self.conv2 = ME.MinkowskiConvolution(out_channels, out_channels, kernel_size=3, dimension=dimension)
    self.bn2 = ME.MinkowskiBatchNorm(out_channels)
    self.relu = ME.MinkowskiReLU(inplace=True)
    if in_channels != out_channels:
      self.shortcut = ME.MinkowskiConvolution(in_channels, out_channels, kernel_size=1, dimension=dimension)
    else:
      self.shortcut = None

  def forward(self, x):
    identity = x
    out = self.relu(self.bn1(self.conv1(x)))
    out = self.bn2(self.conv2(out))
    if self.shortcut is not None:
      identity = self.shortcut(identity)
    return self.relu(out + identity)

class MinkUNetClassifier(nn.Module):
  def __init__(self, in_channels=4, out_channels=1, dimension=2, base_filters=16, num_strides=3):
    super().__init__()
    self.conv0 = ME.MinkowskiConvolution(in_channels, base_filters, kernel_size=3, dimension=dimension)
    self.encoder = nn.ModuleList()
    ch = base_filters
    for _ in range(num_strides):
      self.encoder.append(ResidualBlock(ch, ch*2, dimension))
      self.encoder.append(ME.MinkowskiConvolution(ch*2, ch*2, kernel_size=2, stride=2, dimension=dimension))
      ch *= 2
    self.bottleneck = ResidualBlock(ch, ch, dimension)
    self.decoder = nn.ModuleList()
    for i in range(num_strides):
      up = ch//2
      self.decoder.append(ME.MinkowskiConvolutionTranspose(ch, up, kernel_size=2, stride=2, dimension=dimension))
      skip_ch = base_filters * (2 ** (num_strides - i))
      self.decoder.append(ResidualBlock(up + skip_ch, up, dimension))
      ch = up
    self.bn_relu = nn.Sequential(
      ME.MinkowskiBatchNorm(base_filters),
      ME.MinkowskiReLU(inplace=True))
    self.global_pool = ME.MinkowskiGlobalPooling()
    self.linear = ME.MinkowskiLinear(base_filters, out_channels)

  def forward(self, x):
    x = self.conv0(x)
    skips = []
    for i in range(0, len(self.encoder), 2):
      x = self.encoder[i](x)
      skips.append(x)
      x = self.encoder[i+1](x)
    x = self.bottleneck(x)
    for i in range(0, len(self.decoder), 2):
      x = self.decoder[i](x)
      skip = skips.pop()
      x = ME.cat(x, skip)
      x = self.decoder[i+1](x)
    x = self.bn_relu(x)
    x = self.global_pool(x)
    x = self.linear(x)
    return x

def extract_features(flat_img):
  img2d = flat_img.reshape(512, 512)
  rows, cols = np.nonzero(img2d > 0)
  if len(rows) == 0: return None, None
  coords = np.vstack([rows, cols]).T.astype(np.int32)
  adc = img2d[rows, cols]
  norm_rows = (rows - 256.0) / 256.0
  norm_cols = (cols - 256.0) / 256.0
  log_adc = np.log10(adc + 1e-6)
  from scipy.spatial import cKDTree
  tree = cKDTree(coords)
  density = tree.query_ball_point(coords, r=5, return_length=True).astype(np.float32)
  norm_density = np.log1p(density)
  feats = np.vstack([norm_rows, norm_cols, log_adc, norm_density]).T.astype(np.float32)
  return coords, feats

def load_imgs_from_npy(path):
  arr = np.load(path, allow_pickle=False, mmap_mode=None)
  if arr.ndim == 2 and arr.shape == (512, 512):
    arr = arr.reshape(1, -1)
  elif arr.ndim == 3 and arr.shape[1:] == (512, 512):
    arr = arr.reshape(arr.shape[0], -1)
  elif arr.ndim == 1 and arr.size == 512*512:
    arr = arr.reshape(1, -1)
  elif arr.ndim == 2 and arr.shape[1] == 512*512:
    pass
  else:
    raise RuntimeError(f"Unexpected .npy shape {arr.shape}; expected (512,512), (N,512,512), (512*512,), or (N,512*512)")
  if arr.dtype != np.float32:
    arr = arr.astype(np.float32, copy=False)
  return arr

def main():
  start_time = time.time()
  parser = argparse.ArgumentParser(description="Run inference (NPY input).")
  parser.add_argument("--npy", type=str, required=True, help="NPY file with image(s)")
  parser.add_argument("--output", type=str, required=True, help="Output text file")
  weight_group = parser.add_mutually_exclusive_group(required=True)
  weight_group.add_argument("--weights", type=str, help="Path to a weights file")
  weight_group.add_argument("--model", type=str, help="Model name under the 'weights' directory")
  args = parser.parse_args()

  if args.model:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    weights_path = os.path.join(script_dir, "weights", args.model)
    if not os.path.splitext(weights_path)[1]:
      weights_path += ".pth"
  else:
    weights_path = args.weights

  device = torch.device("cpu")
  model = MinkUNetClassifier().to(device)
  model.load_state_dict(torch.load(weights_path, map_location="cpu"))
  model.eval()

  imgs = load_imgs_from_npy(args.npy)
  n = len(imgs)
  logits = np.zeros(n, dtype=np.float32)

  torch.set_num_threads(1)
  with torch.no_grad():
    for i, flat in enumerate(imgs):
      if flat.size == 0:
        continue
      coords, feats = extract_features(flat)
      if coords is None:
        continue
      batched = ME.utils.batched_coordinates([coords])
      feats_t = torch.from_numpy(feats).to(device)
      x = ME.SparseTensor(feats_t, coordinates=batched, device=device)
      logits[i] = float(model(x).F.squeeze().cpu().item())

  score = float(np.mean(logits))
  with open(args.output, "w") as f:
    f.write(f"{score}\n")

  

if __name__ == "__main__":
  main()
