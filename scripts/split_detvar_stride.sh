#!/usr/bin/env bash

usage() {
  cat <<'USAGE'
Usage:
  split_detvar_stride.sh <source_def> <training_def> <nominal_def>
  split_detvar_stride.sh --batch <triples_file>

Creates two SAM definitions from a nominal detector-variation sample using
stride-2 splitting:
  - <training_def>: defname:<source_def> with stride 2
  - <nominal_def>: defname:<source_def> with stride 2 with offset 1

Arguments:
  <source_def>   Source SAM definition (strangeness, beam, dirt, EXT, etc.).
  <training_def> Training subset definition name.
  <nominal_def>  Nominal (complement) subset definition name.

Batch mode:
  Provide a file with whitespace-separated triples per line:
    <source_def> <training_def> <nominal_def>
  Empty lines and lines starting with # are ignored.

Examples:
  ./split_detvar_stride.sh \
    prod_strange_resample_fhc_run2_fhc_reco2_reco2_detvar_nominal \
    nl_strange_detvar_nominal_training_stride2 \
    nl_strange_detvar_nominal

  ./split_detvar_stride.sh \
    prod_numi_fhc_beam_detvar_nominal \
    nl_beam_detvar_nominal_training_stride2 \
    nl_beam_detvar_nominal

  ./split_detvar_stride.sh --batch detvar_triples.txt
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

run_split() {
  local source_def="$1"
  local training_def="$2"
  local nominal_def="$3"

  set -x
  samweb create-definition "${training_def}" "defname: ${source_def} with stride 2"
  samweb create-definition "${nominal_def}" "defname: ${source_def} with stride 2 with offset 1"
  set +x

  echo "Created definitions:"
  echo "  source  : ${source_def}"
  echo "  training: ${training_def}"
  echo "  nominal : ${nominal_def}"

  echo
  echo "Summary (training):"
  samweb list-files --summary "defname: ${training_def}"

  echo
  echo "Summary (nominal):"
  samweb list-files --summary "defname: ${nominal_def}"
}

if [[ ${1:-} == "--batch" ]]; then
  if [[ $# -ne 2 ]]; then
    echo "Error: --batch expects a single filename." >&2
    usage >&2
    exit 1
  fi
  triples_file="$2"
  if [[ ! -f "${triples_file}" ]]; then
    echo "Error: batch file not found: ${triples_file}" >&2
    exit 1
  fi

  while read -r source_def training_def nominal_def extra; do
    if [[ -z ${source_def:-} || ${source_def} == "#"* ]]; then
      continue
    fi
    if [[ -n ${extra:-} || -z ${training_def:-} || -z ${nominal_def:-} ]]; then
      echo "Error: invalid line (expected 3 fields): ${source_def} ${training_def:-} ${nominal_def:-} ${extra:-}" >&2
      exit 1
    fi
    run_split "${source_def}" "${training_def}" "${nominal_def}"
  done < "${triples_file}"
  exit 0
fi

if [[ $# -ne 3 ]]; then
  echo "Error: expected 3 arguments, got $#" >&2
  usage >&2
  exit 1
fi

run_split "$1" "$2" "$3"
