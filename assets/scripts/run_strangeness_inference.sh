#!/usr/bin/env bash
set -euo pipefail
unset PYTHONHOME PYTHONPATH
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

echo "CWD: $(pwd)"
ls -al

fw_script=""
if [[ -n "${FW_SEARCH_PATH:-}" ]]; then
  IFS=':' read -ra dirs <<< "$FW_SEARCH_PATH"
  for d in "${dirs[@]}"; do
    if [[ -f "$d/run_inference.py" ]]; then
      fw_script="$d/run_inference.py"
      break
    fi
  done
fi
if [[ -z "$fw_script" ]]; then
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  candidate="${script_dir}/run_inference.py"
  if [[ -f "$candidate" ]]; then
    fw_script="$candidate"
  fi
fi
if [[ -z "$fw_script" || ! -f "$fw_script" ]]; then
  echo "run_inference.py not found" >&2
  echo "CWD: $(pwd)" >&2
  ls -al >&2
  exit 1
else
  echo "Using run_inference.py at $fw_script"
fi

python3 "$fw_script" "$@"

