#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  create_training_template_shards.sh [--plan <triples_file>] [--target-events <count>] [--full-shards] [--dry-run]

Creates orthogonal training/template SAM definition shards from a triples file.
The default plan is scripts/run1_detvar_cv_shards.txt.
By default the wrapper caps only the training shard at roughly 100000 events.

Each non-comment line in the plan file must contain:
  <source_def> <training_shard_def> <template_shard_def>

Use --dry-run to print the source counts and the samweb commands without
creating any definitions.
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if ! command -v samweb >/dev/null 2>&1; then
  echo "Error: samweb is not available in PATH." >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
plan_file="${script_dir}/run1_detvar_cv_shards.txt"
target_events=100000
full_shards=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plan)
      if [[ $# -lt 2 ]]; then
        echo "Error: --plan expects a filename." >&2
        exit 1
      fi
      plan_file="$2"
      shift 2
      ;;
    --target-events)
      if [[ $# -lt 2 ]]; then
        echo "Error: --target-events expects an integer." >&2
        exit 1
      fi
      target_events="$2"
      shift 2
      ;;
    --full-shards)
      full_shards=1
      shift
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    *)
      echo "Error: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ! -f "${plan_file}" ]]; then
  echo "Error: plan file not found: ${plan_file}" >&2
  exit 1
fi

if (( full_shards == 0 )) && { ! [[ ${target_events} =~ ^[0-9]+$ ]] || (( target_events <= 0 )); }; then
  echo "Error: --target-events must be a positive integer." >&2
  exit 1
fi

run_split_script="${script_dir}/split_detvar_stride.sh"
if [[ ! -f "${run_split_script}" ]]; then
  echo "Error: helper not found: ${run_split_script}" >&2
  exit 1
fi

split_args=()
if (( full_shards == 0 )); then
  split_args+=(--target-events "${target_events}")
fi
if (( dry_run )); then
  split_args+=(--dry-run)
fi

while read -r source_def training_shard_def template_shard_def extra; do
  if [[ -z ${source_def:-} || ${source_def} == "#"* ]]; then
    continue
  fi
  if [[ -n ${extra:-} || -z ${training_shard_def:-} || -z ${template_shard_def:-} ]]; then
    echo "Error: invalid line (expected 3 fields): ${source_def} ${training_shard_def:-} ${template_shard_def:-} ${extra:-}" >&2
    exit 1
  fi
  bash "${run_split_script}" "${split_args[@]}" "${source_def}" "${training_shard_def}" "${template_shard_def}"
done < "${plan_file}"
