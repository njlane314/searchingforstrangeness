#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  query_campaign_numjobs.sh [--xml <xml_file>]... [--threshold <files>] [--include-derived-shards] [--tsv]

Counts files for the input SAM definitions declared in the staged campaign XMLs.
Because the checked-in XMLs use maxfilesperjob=1, the file count is the numjobs
value for that definition.

Defaults:
  XML files : xml/numi_reco2_*_campaign.xml
  Threshold : 5000 files

Recommendations:
  keep    file count is at or below the threshold
  shard   file count is above the threshold
  missing samweb could not resolve the definition

By default, derived train/template shard entities are skipped because they are
usually created later from the source definitions. Use --include-derived-shards
to count them as well.
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
repo_root="$(cd "${script_dir}/.." && pwd)"

threshold=5000
include_derived=0
tsv=0
xml_files=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --xml)
      if [[ $# -lt 2 ]]; then
        echo "Error: --xml expects a filename." >&2
        exit 1
      fi
      xml_files+=("$2")
      shift 2
      ;;
    --threshold)
      if [[ $# -lt 2 ]]; then
        echo "Error: --threshold expects an integer." >&2
        exit 1
      fi
      threshold="$2"
      shift 2
      ;;
    --include-derived-shards)
      include_derived=1
      shift
      ;;
    --tsv)
      tsv=1
      shift
      ;;
    *)
      echo "Error: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if ! [[ ${threshold} =~ ^[0-9]+$ ]] || (( threshold <= 0 )); then
  echo "Error: --threshold must be a positive integer." >&2
  exit 1
fi

if [[ ${#xml_files[@]} -eq 0 ]]; then
  xml_files=("${repo_root}"/xml/numi_reco2_*_campaign.xml)
fi

for xml_file in "${xml_files[@]}"; do
  if [[ ! -f "${xml_file}" ]]; then
    echo "Error: XML file not found: ${xml_file}" >&2
    exit 1
  fi
done

parsed_file="$(mktemp)"
grouped_file="$(mktemp)"
trap 'rm -f "${parsed_file}" "${grouped_file}"' EXIT

awk -v include_derived="${include_derived}" '
  match($0, /<!ENTITY[[:space:]]+(input_[^[:space:]]+)[[:space:]]+"([^"]+)">/, m) {
    entity = m[1]
    samdef = m[2]
    if (!include_derived && entity ~ /_(train|template)_shard$/) {
      next
    }
    print samdef "\t" FILENAME ":" entity
  }
' "${xml_files[@]}" | sort -u > "${parsed_file}"

if [[ ! -s "${parsed_file}" ]]; then
  echo "Error: no input SAM definitions were found in the requested XML files." >&2
  exit 1
fi

awk -F'\t' '
  NR == 1 {
    current_def = $1
    refs = $2
    next
  }
  $1 == current_def {
    refs = refs "," $2
    next
  }
  {
    print current_def "\t" refs
    current_def = $1
    refs = $2
  }
  END {
    if (NR > 0) {
      print current_def "\t" refs
    }
  }
' "${parsed_file}" > "${grouped_file}"

keep_count=0
shard_count=0
missing_count=0

if (( tsv )); then
  printf "numjobs\tdecision\tsuggested_shards\tsam_definition\treferences\n"
else
  printf "Threshold for sharding: %s files\n" "${threshold}"
  printf "numjobs  decision  shards  sam definition\n"
fi

while IFS=$'\t' read -r samdef refs; do
  if count="$(samweb count-files "defname:${samdef}" 2>/dev/null)"; then
    if [[ ${count} =~ ^[0-9]+$ ]]; then
      if (( count > threshold )); then
        decision="shard"
        suggested_shards=$(( (count + threshold - 1) / threshold ))
        shard_count=$(( shard_count + 1 ))
      else
        decision="keep"
        suggested_shards=1
        keep_count=$(( keep_count + 1 ))
      fi
    else
      count="missing"
      decision="missing"
      suggested_shards="-"
      missing_count=$(( missing_count + 1 ))
    fi
  else
    count="missing"
    decision="missing"
    suggested_shards="-"
    missing_count=$(( missing_count + 1 ))
  fi

  if (( tsv )); then
    printf "%s\t%s\t%s\t%s\t%s\n" "${count}" "${decision}" "${suggested_shards}" "${samdef}" "${refs}"
  else
    printf "%-7s  %-8s  %-6s  %s\n" "${count}" "${decision}" "${suggested_shards}" "${samdef}"
    printf "  refs: %s\n" "${refs}"
  fi
done < "${grouped_file}"

if (( ! tsv )); then
  echo
  printf "Summary: keep=%s shard=%s missing=%s\n" "${keep_count}" "${shard_count}" "${missing_count}"
fi
