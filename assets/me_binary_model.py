import MinkowskiEngine as ME
import torch.nn as nn


class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, dimension):
        super().__init__()
        self.conv1 = ME.MinkowskiConvolution(
            in_channels, out_channels, kernel_size=3, dimension=dimension
        )
        self.bn1 = ME.MinkowskiBatchNorm(out_channels)
        self.conv2 = ME.MinkowskiConvolution(
            out_channels, out_channels, kernel_size=3, dimension=dimension
        )
        self.bn2 = ME.MinkowskiBatchNorm(out_channels)
        self.relu = ME.MinkowskiReLU(inplace=True)
        if in_channels != out_channels:
            self.shortcut = ME.MinkowskiConvolution(
                in_channels, out_channels, kernel_size=1, dimension=dimension
            )
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
        self.conv0 = ME.MinkowskiConvolution(
            in_channels, base_filters, kernel_size=3, dimension=dimension
        )
        self.encoder = nn.ModuleList()
        ch = base_filters
        for _ in range(num_strides):
            self.encoder.append(ResidualBlock(ch, ch * 2, dimension))
            self.encoder.append(
                ME.MinkowskiConvolution(
                    ch * 2, ch * 2, kernel_size=2, stride=2, dimension=dimension
                )
            )
            ch *= 2
        self.bottleneck = ResidualBlock(ch, ch, dimension)
        self.decoder = nn.ModuleList()
        for i in range(num_strides):
            up = ch // 2
            self.decoder.append(
                ME.MinkowskiConvolutionTranspose(
                    ch, up, kernel_size=2, stride=2, dimension=dimension
                )
            )
            skip_ch = base_filters * (2 ** (num_strides - i))
            self.decoder.append(ResidualBlock(up + skip_ch, up, dimension))
            ch = up
        self.bn_relu = nn.Sequential(
            ME.MinkowskiBatchNorm(base_filters),
            ME.MinkowskiReLU(inplace=True),
        )
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
    return MinkUNetClassifier()
