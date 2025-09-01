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
from models import load_model

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
  parser.add_argument("--arch", type=str, required=True, help="Model architecture module under 'models'")
  weight_group = parser.add_mutually_exclusive_group(required=True)
  weight_group.add_argument("--weights", type=str, help="Path to a weights file")
  weight_group.add_argument("--model", type=str, help="Model weights name under the 'weights' directory")
  args = parser.parse_args()

  if args.model:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    weights_path = os.path.join(script_dir, "weights", args.model)
    if not os.path.splitext(weights_path)[1]:
      weights_path += ".pth"
  else:
    weights_path = args.weights

  device = torch.device("cpu")
  model = load_model(args.arch).to(device)
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
