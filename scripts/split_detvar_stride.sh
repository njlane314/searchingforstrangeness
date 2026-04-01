#!/usr/bin/env bash

usage() {
  cat <<'USAGE'
Usage:
  split_detvar_stride.sh <source_def> <training_shard_def> <template_shard_def>
  split_detvar_stride.sh --batch <triples_file>

Creates two orthogonal SAM definition shards from a source sample using
stride-2 splitting:
  - <training_shard_def>: defname:<source_def> with stride 2
  - <template_shard_def>: defname:<source_def> with stride 2 with offset 1

Arguments:
  <source_def>          Source SAM definition (strangeness, beam, dirt, EXT, etc.).
  <training_shard_def>  Training subset definition name.
  <template_shard_def>  Template subset definition name.

Batch mode:
  Provide a file with whitespace-separated triples per line:
    <source_def> <training_shard_def> <template_shard_def>
  Empty lines and lines starting with # are ignored.

Examples:
  ./split_detvar_stride.sh \
    detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2 \
    nl_run1_fhc_strangeness_detvar_cv_train_shard \
    nl_run1_fhc_strangeness_detvar_cv_template_shard

  ./split_detvar_stride.sh \
    prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2 \
    nl_run1_fhc_beam_detvar_cv_train_shard \
    nl_run1_fhc_beam_detvar_cv_template_shard

  ./split_detvar_stride.sh --batch scripts/run1_detvar_cv_shards.txt
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
  local training_shard_def="$2"
  local template_shard_def="$3"

  set -x
  samweb create-definition "${training_shard_def}" "defname: ${source_def} with stride 2"
  samweb create-definition "${template_shard_def}" "defname: ${source_def} with stride 2 with offset 1"
  set +x

  echo "Created definitions:"
  echo "  source         : ${source_def}"
  echo "  training shard : ${training_shard_def}"
  echo "  template shard : ${template_shard_def}"

  echo
  echo "Summary (training shard):"
  samweb list-files --summary "defname: ${training_shard_def}"

  echo
  echo "Summary (template shard):"
  samweb list-files --summary "defname: ${template_shard_def}"
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

  while read -r source_def training_shard_def template_shard_def extra; do
    if [[ -z ${source_def:-} || ${source_def} == "#"* ]]; then
      continue
    fi
    if [[ -n ${extra:-} || -z ${training_shard_def:-} || -z ${template_shard_def:-} ]]; then
      echo "Error: invalid line (expected 3 fields): ${source_def} ${training_shard_def:-} ${template_shard_def:-} ${extra:-}" >&2
      exit 1
    fi
    run_split "${source_def}" "${training_shard_def}" "${template_shard_def}"
  done < "${triples_file}"
  exit 0
fi

if [[ $# -ne 3 ]]; then
  echo "Error: expected 3 arguments, got $#" >&2
  usage >&2
  exit 1
fi

run_split "$1" "$2" "$3"
