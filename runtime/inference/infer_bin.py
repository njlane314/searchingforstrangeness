#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence

import numpy as np

from inference_config import InferenceConfig, parse_arch, resolve_device
from preprocessing import prepare_sparse_plane, validate_feature_schema
from protocol import read_sparse_request, write_results


PLANE_NAMES = ("u", "v", "w")


@dataclass(frozen=True)
class CliArguments:
    input_path: str
    output_path: str
    metrics_path: str
    width: int
    height: int
    arch: str
    weights_path: str


def parse_arguments(argv: Optional[Sequence[str]] = None) -> CliArguments:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="input_path", required=True)
    parser.add_argument("--out", dest="output_path", required=True)
    parser.add_argument("--metrics", dest="metrics_path", required=True)
    parser.add_argument("--W", dest="width", type=int, required=True)
    parser.add_argument("--H", dest="height", type=int, required=True)
    parser.add_argument("--arch", default="")
    parser.add_argument("--weights", dest="weights_path", required=True)
    return CliArguments(**vars(parser.parse_args(argv)))


def _sync_if_cuda(torch: Any, device: str) -> None:
    if device.startswith("cuda"):
        torch.cuda.synchronize(device)


def _max_rss_mb() -> float:
    try:
        import resource

        # The production runtime is Linux, where ru_maxrss is reported in KiB.
        return float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / 1024.0
    except Exception:
        return 0.0


def _feature_schema_notice(payload_feature_dim: int, model_in_ch: int) -> Optional[str]:
    if payload_feature_dim == 2 and model_in_ch in (3, 4):
        return (
            f"[infer_bin] sparse payload feature_dim={payload_feature_dim}; padding "
            f"missing dead feature(s) to match model_in_ch={model_in_ch}"
        )
    return None


def _setup_log(
    *,
    device: str,
    config: InferenceConfig,
    payload_feature_dim: int,
    nnz: Sequence[int],
    availability: Sequence[float],
) -> str:
    return (
        f"[infer_bin] device={device} thr={config.threshold} "
        f"sparse_feature_dim={payload_feature_dim} model_in_ch={config.in_ch} "
        f"nnz(u,v,w)={tuple(nnz)} "
        f"avail={tuple(int(value) for value in availability)}"
    )


def _write_metrics(
    path: str,
    *,
    started: float,
    setup_finished: float,
    inference_finished: float,
    write_finished: float,
    config: InferenceConfig,
) -> None:
    with open(path, "w") as metrics:
        metrics.write(f"t_total_ms={(write_finished - started) * 1000.0:.3f}\n")
        metrics.write(f"t_setup_ms={(setup_finished - started) * 1000.0:.3f}\n")
        metrics.write(
            f"t_infer_ms={(inference_finished - setup_finished) * 1000.0:.3f}\n"
        )
        metrics.write(
            f"t_post_ms={(write_finished - inference_finished) * 1000.0:.3f}\n"
        )
        metrics.write(f"max_rss_mb={_max_rss_mb():.3f}\n")
        metrics.write(f"seed={config.seed}\n")
        metrics.write(f"thr={config.threshold}\n")


def _runtime_dependencies():
    """Import the production stack before starting the legacy metrics clock."""

    import torch

    from me_compat import make_sparse_tensor
    from model_loading import build_llr_model

    return torch, make_sparse_tensor, build_llr_model


def run(arguments: CliArguments) -> int:
    # HEAD imported torch/ME/model code before main() started its metrics clock.
    # Keep that timing boundary while retaining a lightweight module import.
    torch, make_sparse_tensor, build_llr_model = _runtime_dependencies()

    started = time.time()
    config = parse_arch(arguments.arch)
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    device = resolve_device(
        config.device,
        cuda_available=bool(torch.cuda.is_available()),
    )

    request = read_sparse_request(
        arguments.input_path,
        expected_planes=len(PLANE_NAMES),
        expected_height=arguments.height,
        expected_width=arguments.width,
    )
    validate_feature_schema(request.feature_dim, config.in_ch)

    schema_notice = _feature_schema_notice(request.feature_dim, config.in_ch)
    if schema_notice is not None:
        print(schema_notice, flush=True)

    model = build_llr_model(
        weights_path=arguments.weights_path,
        device=device,
        config=config,
        plane_names=PLANE_NAMES,
    )

    prepared_planes = [
        prepare_sparse_plane(
            plane.coordinates,
            plane.features,
            payload_feature_dim=request.feature_dim,
            model_in_ch=config.in_ch,
            height=request.height,
            width=request.width,
            threshold=config.threshold,
        )
        for plane in request.planes
    ]

    inputs: Dict[str, Any] = {
        name: make_sparse_tensor(plane, device=device)
        for name, plane in zip(PLANE_NAMES, prepared_planes)
    }
    availability = [plane.available for plane in prepared_planes]
    available_mask = torch.tensor(
        [availability],
        dtype=torch.float32,
        device=device,
    )
    setup_finished = time.time()

    print(
        _setup_log(
            device=device,
            config=config,
            payload_feature_dim=request.feature_dim,
            nnz=[plane.nnz for plane in prepared_planes],
            availability=availability,
        ),
        flush=True,
    )

    _sync_if_cuda(torch, device)
    with torch.no_grad():
        logits = model(inputs, available_mask=available_mask)
    _sync_if_cuda(torch, device)
    inference_finished = time.time()

    logit = float(logits.reshape(-1)[0].detach().cpu().item())
    if not math.isfinite(logit):
        raise RuntimeError(f"Model produced a non-finite logit: {logit}")
    write_results(arguments.output_path, [logit])
    write_finished = time.time()

    _write_metrics(
        arguments.metrics_path,
        started=started,
        setup_finished=setup_finished,
        inference_finished=inference_finished,
        write_finished=write_finished,
        config=config,
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    return run(parse_arguments(argv))


if __name__ == "__main__":
    raise SystemExit(main())
