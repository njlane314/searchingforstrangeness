import argparse
import numpy as np
import uproot
import torch
import torch.nn as nn
import MinkowskiEngine as ME


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
    def __init__(self, in_channels=1, out_channels=1, dimension=2, base_filters=16, num_strides=3):
        super().__init__()
        self.conv0 = ME.MinkowskiConvolution(in_channels, base_filters, kernel_size=3, dimension=dimension)
        self.encoder = nn.ModuleList()
        ch = base_filters
        for _ in range(num_strides):
            self.encoder.append(ResidualBlock(ch, ch * 2, dimension))
            self.encoder.append(ME.MinkowskiConvolution(ch * 2, ch * 2, kernel_size=2, stride=2, dimension=dimension))
            ch *= 2
        self.bottleneck = ResidualBlock(ch, ch, dimension)
        self.decoder = nn.ModuleList()
        for i in range(num_strides):
            up = ch // 2
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
            x = self.encoder[i + 1](x)
        x = self.bottleneck(x)
        for i in range(0, len(self.decoder), 2):
            x = self.decoder[i](x)
            skip = skips.pop()
            x = ME.cat(x, skip)
            x = self.decoder[i + 1](x)
        x = self.bn_relu(x)
        x = self.global_pool(x)
        x = self.linear(x)
        return x


def build_model():
    return MinkUNetClassifier(in_channels=1, out_channels=1, dimension=2)


def image_to_sparse_tensor(img_flat, width, height, device):
    img_flat = np.asarray(img_flat, dtype=np.float32)
    if img_flat.size != width * height:
        raise RuntimeError(f"image size {img_flat.size} != width*height {width*height}")
    idx = np.nonzero(img_flat)[0]
    if idx.size == 0:
        return None
    rows = idx // width
    cols = idx % width
    feats = img_flat[idx].reshape(-1, 1)
    batch = np.zeros_like(rows, dtype=np.int32)
    coords = np.stack([batch, rows.astype(np.int32), cols.astype(np.int32)], axis=1)
    coords_t = torch.from_numpy(coords).to(device=device, dtype=torch.int32)
    feats_t = torch.from_numpy(feats).to(device=device, dtype=torch.float32)
    st = ME.SparseTensor(coordinates=coords_t, features=feats_t)
    return st


def predict_score(model, st, random_weights):
    if st is None:
        return 0.0
    if random_weights:
        return float(torch.rand(1).item())
    with torch.no_grad():
        out = model(st)
    if isinstance(out, ME.SparseTensor):
        feats = out.F
    elif isinstance(out, torch.Tensor):
        feats = out
    else:
        raise RuntimeError(f"Unsupported model output type: {type(out)}")
    logit = feats.view(-1)[0]
    prob = torch.sigmoid(logit).item()
    return float(prob)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--checkpoint", default="")
    parser.add_argument("--tree", default="imageprod/ImageDump")
    parser.add_argument("--random-weights", action="store_true")
    parser.add_argument("--max-entries", type=int, default=-1)
    args = parser.parse_args()

    device = torch.device("cpu")

    model = build_model()
    model.to(device)
    if not args.random_weights:
        if not args.checkpoint:
            raise RuntimeError("checkpoint path is required unless --random-weights is set")
        ckpt = torch.load(args.checkpoint, map_location=device)
        state = ckpt.get("state_dict", ckpt)
        model.load_state_dict(state)
    model.eval()

    with uproot.open(args.input) as fin:
        tree = fin[args.tree]
        runs = tree["run"].array(library="np")
        subruns = tree["subrun"].array(library="np")
        events = tree["event"].array(library="np")
        views = tree["view"].array(library="np")
        widths = tree["width"].array(library="np")
        heights = tree["height"].array(library="np")
        images = tree["image"].array(library="np")

    n_total = len(events)
    if args.max_entries >= 0:
        n_entries = min(n_total, args.max_entries)
    else:
        n_entries = n_total

    out_run = []
    out_subrun = []
    out_event = []
    out_view = []
    out_score = []

    for i in range(n_entries):
        img_flat = np.array(images[i], dtype=np.float32)
        w = int(widths[i])
        h = int(heights[i])
        st = image_to_sparse_tensor(img_flat, w, h, device)
        score = predict_score(model, st, args.random_weights)
        out_run.append(int(runs[i]))
        out_subrun.append(int(subruns[i]))
        out_event.append(int(events[i]))
        out_view.append(int(views[i]))
        out_score.append(float(score))

    with uproot.recreate(args.output) as fout:
        fout["ImageScores"] = {
            "run": np.asarray(out_run, dtype="i4"),
            "subrun": np.asarray(out_subrun, dtype="i4"),
            "event": np.asarray(out_event, dtype="i4"),
            "view": np.asarray(out_view, dtype="i4"),
            "score": np.asarray(out_score, dtype="f4"),
        }


if __name__ == "__main__":
    main()
