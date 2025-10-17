#!/usr/bin/env bash
set -euo pipefail

# Minimal wrapper invoked by InferenceEngine::runInferenceDetailed
# Forwards all arguments to the Python ME feature extractor.
#
# Environment overrides:
#   PY      : python executable (default: python3)
#   SCRIPT  : path to infer_bin.py (default: /app/infer_bin.py)

PY="${PY:-python3}"
SCRIPT="${SCRIPT:-/app/infer_bin.py}"

exec "$PY" "$SCRIPT" "$@"
