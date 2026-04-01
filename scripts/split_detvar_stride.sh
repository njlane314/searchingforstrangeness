#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  split_detvar_stride.sh [--target-events <count>] [--dry-run] <source_def> <training_shard_def> <template_shard_def>
  split_detvar_stride.sh [--target-events <count>] [--dry-run] --batch <triples_file>

Creates two orthogonal SAM definition shards from a source sample using
stride-2 splitting:
  - <training_shard_def>: defname:<source_def> with stride 2
  - <template_shard_def>: defname:<source_def> with stride 2 with offset 1

When --target-events is set, each shard is additionally capped with a
file-level `with limit <N>` chosen from the shard's average events per file so
the output lands near the requested event count.

Arguments:
  <source_def>          Source SAM definition (strangeness, beam, dirt, EXT, etc.).
  <training_shard_def>  Training subset definition name.
  <template_shard_def>  Template subset definition name.

Options:
  --target-events <count>  Approximate target events per shard.
  --dry-run                Print the derived commands without creating definitions.

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
    --target-events 100000 \
    prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2 \
    nl_run1_fhc_beam_detvar_cv_train_shard \
    nl_run1_fhc_beam_detvar_cv_template_shard

  ./split_detvar_stride.sh --batch scripts/run1_detvar_cv_shards.txt
USAGE
}

extract_summary_value() {
  local summary="$1"
  local key="$2"

  awk -F'\t' -v label="${key}:" '
    $1 == label {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
      print $2
      exit
    }
  ' <<< "${summary}"
}

query_summary() {
  local query="$1"
  samweb list-files --summary "${query}"
}

build_shard_query() {
  local source_def="$1"
  local offset="$2"
  local base_query="defname: ${source_def} with stride 2"
  local base_summary
  local full_file_count
  local full_event_count
  local file_limit
  local final_query
  local final_summary
  local result_file_count
  local result_event_count

  if (( offset == 1 )); then
    base_query="${base_query} with offset 1"
  fi

  base_summary="$(query_summary "${base_query}")"
  full_file_count="$(extract_summary_value "${base_summary}" "File count")"
  full_event_count="$(extract_summary_value "${base_summary}" "Event count")"

  if ! [[ ${full_file_count} =~ ^[0-9]+$ && ${full_event_count} =~ ^[0-9]+$ ]]; then
    echo "Error: could not parse samweb summary for query: ${base_query}" >&2
    exit 1
  fi

  file_limit="${full_file_count}"
  final_query="${base_query}"
  result_file_count="${full_file_count}"
  result_event_count="${full_event_count}"

  if [[ -n ${target_events:-} && ${target_events} -gt 0 && ${full_event_count} -gt ${target_events} ]]; then
    file_limit=$(( (target_events * full_file_count + full_event_count - 1) / full_event_count ))
    if (( file_limit < 1 )); then
      file_limit=1
    elif (( file_limit > full_file_count )); then
      file_limit="${full_file_count}"
    fi

    if (( file_limit < full_file_count )); then
      final_query="${base_query} with limit ${file_limit}"
      final_summary="$(query_summary "${final_query}")"
      result_file_count="$(extract_summary_value "${final_summary}" "File count")"
      result_event_count="$(extract_summary_value "${final_summary}" "Event count")"

      if ! [[ ${result_file_count} =~ ^[0-9]+$ && ${result_event_count} =~ ^[0-9]+$ ]]; then
        echo "Error: could not parse samweb summary for query: ${final_query}" >&2
        exit 1
      fi
    fi
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${base_query}" \
    "${final_query}" \
    "${full_file_count}" \
    "${full_event_count}" \
    "${file_limit}" \
    "${result_file_count}" \
    "${result_event_count}"
}

describe_shard() {
  local label="$1"
  local shard_name="$2"
  local result_file_count="$3"
  local full_file_count="$4"
  local result_event_count="$5"
  local full_event_count="$6"

  if (( result_file_count == full_file_count )); then
    printf '%-17s %s (%s files, %s events)\n' "${label}" "${shard_name}" "${result_file_count}" "${result_event_count}"
  else
    printf '%-17s %s (%s/%s files, %s/%s events)\n' \
      "${label}" \
      "${shard_name}" \
      "${result_file_count}" \
      "${full_file_count}" \
      "${result_event_count}" \
      "${full_event_count}"
  fi
}

run_split() {
  local source_def="$1"
  local training_shard_def="$2"
  local template_shard_def="$3"
  local train_base_query
  local train_query
  local train_full_files
  local train_full_events
  local train_limit
  local train_result_files
  local train_result_events
  local template_base_query
  local template_query
  local template_full_files
  local template_full_events
  local template_limit
  local template_result_files
  local template_result_events
  local source_file_count
  local source_event_count

  IFS=$'\t' read -r train_base_query train_query train_full_files train_full_events train_limit train_result_files train_result_events <<< "$(build_shard_query "${source_def}" 0)"
  IFS=$'\t' read -r template_base_query template_query template_full_files template_full_events template_limit template_result_files template_result_events <<< "$(build_shard_query "${source_def}" 1)"

  source_file_count=$(( train_full_files + template_full_files ))
  source_event_count=$(( train_full_events + template_full_events ))

  echo "Source definition  : ${source_def}"
  echo "Source file count  : ${source_file_count}"
  echo "Source event count : ${source_event_count}"
  if [[ -n ${target_events:-} && ${target_events} -gt 0 ]]; then
    echo "Target events/shard: ${target_events}"
  fi
  describe_shard "Training shard :" "${training_shard_def}" "${train_result_files}" "${train_full_files}" "${train_result_events}" "${train_full_events}"
  describe_shard "Template shard :" "${template_shard_def}" "${template_result_files}" "${template_full_files}" "${template_result_events}" "${template_full_events}"

  if (( dry_run )); then
    printf 'samweb create-definition "%s" "%s"\n' "${training_shard_def}" "${train_query}"
    printf 'samweb create-definition "%s" "%s"\n' "${template_shard_def}" "${template_query}"
    echo
    return 0
  fi

  set -x
  samweb create-definition "${training_shard_def}" "${train_query}"
  samweb create-definition "${template_shard_def}" "${template_query}"
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
  echo
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if ! command -v samweb >/dev/null 2>&1; then
  echo "Error: samweb is not available in PATH." >&2
  exit 1
fi

target_events=""
dry_run=0
triples_file=""
positional=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --batch)
      if [[ $# -lt 2 ]]; then
        echo "Error: --batch expects a single filename." >&2
        usage >&2
        exit 1
      fi
      triples_file="$2"
      shift 2
      ;;
    --target-events)
      if [[ $# -lt 2 ]]; then
        echo "Error: --target-events expects an integer." >&2
        usage >&2
        exit 1
      fi
      target_events="$2"
      shift 2
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    *)
      positional+=("$1")
      shift
      ;;
  esac
done

if [[ -n ${target_events} ]] && { ! [[ ${target_events} =~ ^[0-9]+$ ]] || (( target_events <= 0 )); }; then
  echo "Error: --target-events must be a positive integer." >&2
  exit 1
fi

if [[ -n ${triples_file} ]]; then
  if [[ ${#positional[@]} -ne 0 ]]; then
    echo "Error: positional arguments are not allowed with --batch." >&2
    usage >&2
    exit 1
  fi
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

if [[ ${#positional[@]} -ne 3 ]]; then
  echo "Error: expected 3 arguments, got ${#positional[@]}" >&2
  usage >&2
  exit 1
fi

run_split "${positional[0]}" "${positional[1]}" "${positional[2]}"
