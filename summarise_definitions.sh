#!/usr/bin/env bash

set -u

TIMEOUT=60

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 definitions.txt" >&2
  exit 1
fi

defs_file="$1"
if [[ ! -f "$defs_file" ]]; then
  echo "Definitions file '$defs_file' not found" >&2
  exit 1
fi

while IFS= read -r def || [[ -n "$def" ]]; do
  [[ -z "$def" || "$def" =~ ^# ]] && continue

  echo "================================================================"
  echo "Definition: $def"

  desc="$(samweb describe-definition "$def" 2>/dev/null || true)"

  if [[ -z "$desc" ]]; then
    echo "  Definition summary:"
    echo "    (samweb describe-definition failed)"
    echo
    echo "  Representative file meta-data:"
    echo "    (skipped because describe-definition failed)"
    echo
    continue
  fi

  owner="$(awk -F: '/^Owner/ {gsub(/^[ \t]+/, "", $2); print $2}' <<<"$desc" | head -n1)"
  group="$(awk -F: '/^Group/ {gsub(/^[ \t]+/, "", $2); print $2}' <<<"$desc" | head -n1)"
  dims="$(awk -F: '/^Definition:/ {$1=""; sub(/^ /,""); print}' <<<"$desc" | head -n1)"

  echo "  Definition summary:"
  [[ -n "$owner" ]] && echo "    Owner   : $owner"
  [[ -n "$group" ]] && echo "    Group   : $group"

  ver="$(grep -o 'ub_project.version [^ ]*' <<<"$dims" 2>/dev/null | awk '{print $2}' | head -n1 || true)"
  [[ -n "$ver" ]] && echo "    ub_project.version : $ver"

  echo "    Dimensions:"
  echo "      $dims"
  echo

  echo "  Representative file meta-data:"

  if grep -q 'minus isparentof' <<<"$dims"; then
    echo "    Skipped (definition uses 'minus isparentof'; listing files can hang)."
    echo
    continue
  fi

  rep_file="$(timeout "$TIMEOUT" samweb list-files "defname: $def with limit 1" 2>/dev/null | head -n1 || true)"
  rc=$?

  if [[ $rc -eq 124 ]]; then
    echo "    Timed out after ${TIMEOUT}s trying to get a representative file; skipping."
    echo
    continue
  fi

  if [[ -z "$rep_file" ]]; then
    echo "    No files found (or samweb list-files failed)."
    echo
    continue
  fi

  summary="$(timeout "$TIMEOUT" samweb list-files --summary "defname: $def" 2>/dev/null || true)"
  if [[ $? -eq 124 ]]; then
    echo "    File summary        : (timed out after ${TIMEOUT}s)"
  elif [[ -n "$summary" ]]; then
    echo "    File summary        : $summary"
  fi

  echo "    File: $rep_file"

  meta="$(samweb get-metadata "$rep_file" 2>/dev/null || true)"
  if [[ -z "$meta" ]]; then
    echo "    (samweb get-metadata failed)"
    echo
    continue
  fi

  name="$(grep '^ub_project.name'    <<<"$meta" | awk -F= '{print $2}' | xargs || true)"
  stage="$(grep '^ub_project.stage'  <<<"$meta" | awk -F= '{print $2}' | xargs || true)"
  vers2="$(grep '^ub_project.version' <<<"$meta" | awk -F= '{print $2}' | xargs || true)"

  [[ -n "$name"  ]] && printf '    %-20s : %s\n' 'ub_project.name'    "$name"
  [[ -n "$stage" ]] && printf '    %-20s : %s\n' 'ub_project.stage'   "$stage"
  [[ -n "$vers2" ]] && printf '    %-20s : %s\n' 'ub_project.version' "$vers2"

  echo

done < "$defs_file"
