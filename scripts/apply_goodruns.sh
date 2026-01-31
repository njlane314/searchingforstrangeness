#!/usr/bin/env bash

usage() {
  cat <<'USAGE'
Usage:
  apply_goodruns.sh <source_def> <goodruns_def> [condition]

Creates a new SAM definition by applying a good-runs condition to a source
definition.

Arguments:
  <source_def>   Source SAM definition.
  <goodruns_def> Output definition name.
  [condition]    Optional condition (default: "goodruns: 1").

Example:
  ./apply_goodruns.sh prod_numi_fhc_beam_run1_reco2 \
    nl_numi_fhc_beam_run1_reco2_goodruns
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

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Error: expected 2 or 3 arguments, got $#" >&2
  usage >&2
  exit 1
fi

SOURCE_DEF="$1"
GOODRUNS_DEF="$2"
CONDITION="${3:-goodruns: 1}"

samweb create-definition "${GOODRUNS_DEF}" "defname: ${SOURCE_DEF} and (${CONDITION})"

echo "Created definition: ${GOODRUNS_DEF}"
echo "Condition: ${CONDITION}"

echo
echo "Summary (goodruns):"
samweb list-files --summary "defname: ${GOODRUNS_DEF}"
