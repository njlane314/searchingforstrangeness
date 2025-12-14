#!/usr/bin/env bash
set -euo pipefail

DEFS=(
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_v2_sample0_goodruns"
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample1_goodruns"
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample2_goodruns"
  "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_sample3_goodruns"
  "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample0_goodruns"
  "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_sample1_goodruns"
  "nl_extnumi_run1_all_mcc9_goodruns_reco2_p0"
  "nl_extnumi_run1_all_mcc9_goodruns_reco2_p1"
  "nl_extnumi_run1_all_mcc9_goodruns_reco2_p2"
  "prod_strange_resample_fhc_run1_fhc_reco2_reco2_goodruns"
)

MAXFILESPERJOB=1

declare -A CACHE

ceil_div() { local a="$1" b="$2"; echo $(((a + b - 1) / b)); }

get_file_count() {
  local def="$1"
  if [[ -n "${CACHE[${def}]:-}" ]]; then
    printf '%s' "${CACHE[${def}]}"
    return 0
  fi

  local out fc
  if ! out="$(samweb list-files --summary "defname: ${def}" 2>&1)"; then
    echo "WARN: samweb failed for '${def}': ${out}" >&2
    CACHE["${def}"]=""
    printf ''
    return 0
  fi

  fc="$(awk '/^[[:space:]]*File count:/{print $NF; exit}' <<<"${out}" | tr -d $'\r')"
  [[ "${fc}" =~ ^[0-9]+$ ]] || fc=""
  CACHE["${def}"]="${fc}"
  printf '%s' "${fc}"
}

printf "%-70s  %12s  %10s\n" "samdef" "files" "jobs"
printf "%-70s  %12s  %10s\n" "----------------------------------------------------------------------" "------------" "----------"

for def in "${DEFS[@]}"; do
  fc="$(get_file_count "${def}")"
  if [[ -n "${fc}" ]]; then
    jobs="$(ceil_div "${fc}" "${MAXFILESPERJOB}")"
    printf "%-70s  %12s  %10s\n" "${def}" "${fc}" "${jobs}"
  else
    printf "%-70s  %12s  %10s\n" "${def}" "NA" "NA"
  fi

done
