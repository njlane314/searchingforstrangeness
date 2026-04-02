#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  goodruns-from-samples.sh

Creates standardized `nl_*_goodruns` SAM definitions for the campaign samples
listed in `SAMPLES`, using the hardcoded MCC9 good-runs definitions.

This script is intentionally explicit: every `goodruns.sh` invocation is listed
as a raw command below rather than generated algorithmically.
USAGE
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
apply_script="${script_dir}/goodruns.sh"

if ! command -v samweb >/dev/null 2>&1; then
  echo "Error: samweb is not available in PATH." >&2
  exit 1
fi

if [[ ! -x "${apply_script}" ]]; then
  echo "Error: ${apply_script} is not executable." >&2
  exit 1
fi

# Run 1A
"${apply_script}" "run1_beamon_opentrigger_numi_pandora_3v_a_reprocessing_v08_00_00_84_reco2_reco2_beam_good_good_runs" "nl_run1a_open_trigger_beam_pandora_3v_a_reprocessing_goodruns" "defname: goodruns_mcc9_run1_open_trigger_hardcoded"

# Run 1
"${apply_script}" "prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2" "nl_run1_ext_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prod_mcc9_v08_00_00_45_numi_reco2_run1_beam_good_reco2" "nl_run1_data_goodruns" "defname: goodruns_mcc9_run1_hardcoded"

# Run 1 FHC
"${apply_script}" "prodgenie_numi_uboone_overlay_fhc_mcc9_run1_v28_all_snapshot" "nl_run1_fhc_beam_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie__numi_nue_overlay_mcc9_v08_00_00_48_CV_reco2_run1_reco2" "nl_run1_fhc_nue_cv_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot" "nl_run1_fhc_dirt_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "run1_numi_fhc_CCNCPi0_overlay_v08_00_00_66_reco2_fv_reco2_reco2" "nl_run1_fhc_ccncpi0_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2" "nl_run1_fhc_detvar_cv_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2" "nl_run1_fhc_detvar_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2" "nl_run1_fhc_detvar_ly_down_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2" "nl_run1_fhc_detvar_sce_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2" "nl_run1_fhc_detvar_recomb2_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2" "nl_run1_fhc_detvar_wiremodx_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2" "nl_run1_fhc_detvar_wiremodyz_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2" "nl_run1_fhc_detvar_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2" "nl_run1_fhc_detvar_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_intrinsic_nue_Reco2_LYRayleigh_v08_00_00_48_rerun_reco2_run1_reco2" "nl_run1_fhc_detvar_nue_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "NuMI_run1_intrinsic_nue_reco2_LYDown_reco2_LYDown_reco2" "nl_run1_fhc_detvar_nue_ly_down_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_mcc9_v08_00_00_48_TPCDetVar_run1_reco2_nue_SCE_reco2" "nl_run1_fhc_detvar_nue_sce_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_intrinsic_nue_Reco2_run1_LYAttenuation_v08_00_00_48_reco2_run1_reco2" "nl_run1_fhc_detvar_nue_ly_attenuation_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_mcc9_v08_00_00_48_TPCDetVar_run1_reco2_nue_Recomb2_reco2" "nl_run1_fhc_detvar_nue_recomb2_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nue_Reco2_run1_WireModX_v08_00_00_48_fixedinput_reco2_run1_reco2" "nl_run1_fhc_detvar_nue_wiremodx_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_nue_intrinsic_Reco2_run1_WireModYZ_v08_00_00_48_fixed_input_reco2_run1_reco2" "nl_run1_fhc_detvar_nue_wiremodyz_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie__numi_nue_Reco2_run1_WireModThetaXZ_v08_00_00_48_reco2_reco2_run1_reco2" "nl_run1_fhc_detvar_nue_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie__numi_nue_Reco2_run1_WireModThetaYZ_with_sigma_splines_v08_00_00_48_reco2_reco2_run1_reco2" "nl_run1_fhc_detvar_nue_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run1_hardcoded"

# Run 1 RHC
"${apply_script}" "prodgenie_numi_rhc_nu_overlay_v08_00_00_54_run1_reco2_reco2" "nl_run1_rhc_beam_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prodgenie_numi_rhc_intrinsic_nue_overlay_v08_00_00_54_run1_reco2_reco2_reco2" "nl_run1_rhc_nue_goodruns" "defname: goodruns_mcc9_run1_hardcoded"
"${apply_script}" "prod_extunbiased_numi_rhc_dirt_overlay_run1_reco2_v08_00_00_67_reco2" "nl_run1_rhc_dirt_goodruns" "defname: goodruns_mcc9_run1_hardcoded"

# Run 2a FHC
"${apply_script}" "prodgenie_numi_overlay_detvar_CV_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_beam_cv_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "run2_numi_nu_overlay_pandora_unified_reco2_run2a_fhc_reco2" "nl_run2a_fhc_beam_pandora_unified_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_CV_run2_FHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2a_fhc_nue_cv_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_fhc_dirt_overlay_pandora_reco2_run2_reco2" "nl_run2a_fhc_dirt_pandora_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prod_extnumi_swizzle_inclusive_v4_run2a_reco2_run2a_all_reco2" "nl_run2a_fhc_ext_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prod_numi_swizzle_inclusive_v4_run2_reco2_run2a_beam_good_reco2" "nl_run2a_fhc_data_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_LYAttenuation_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_ly_attenuation_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_LYRayleigh_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_LYDown_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_ly_down_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_SCE_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_sce_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_Recomb2_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_recomb2_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_detvar_WireModX_standard_FHC_nu_overlay_400k_run2_reco2_reco2_reco2" "nl_run2a_fhc_detvar_wiremodx_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_WireModYZ_run2_FHC_standard_nu_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_wiremodyz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModThetaXZ_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_pandora_reco2_run2_FHC_WireModThetaYZ_withSplines_reco2" "nl_run2a_fhc_detvar_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_CV_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_cv_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_LYAttenuation_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_ly_attenuation_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_LYRayleigh_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_LYDown_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_ly_down_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_SCE_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_sce_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_pandora_reco2_run2_FHC_Recomb2_reco2" "nl_run2a_fhc_detvar_nue_recomb2_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_WireModX_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_wiremodx_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_CV_run2_FHC_intrinsic_nue_WireModYZ_reco2_run2_reco2" "nl_run2a_fhc_detvar_nue_wiremodyz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_WireModThetaXZ_run2_FHC_reco2_v08_00_00_55_run2_reco2" "nl_run2a_fhc_detvar_nue_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_pandora_reco2_run2_FHC_WireModThetaYZ_withSplines_reco2" "nl_run2a_fhc_detvar_nue_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"

# Run 2b RHC
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_CV_run2_RHC_reco2_v08_00_00_55_run2_reco2" "nl_run2b_rhc_beam_cv_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "run2_numi_nu_overlay_pandora_unified_reco2_run2b_rhc_reco2" "nl_run2b_rhc_beam_pandora_unified_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_CV_run2_RHC_nue_reco2_v08_00_00_55_run2_run2_reco2_goodruns" "nl_run2b_rhc_nue_cv_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_rhc_dirt_overlay_pandora_reco2_run2_reco2" "nl_run2b_rhc_dirt_pandora_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prod_extnumi_swizzle_crt_inclusive_v4b_offbeam_run2_reco2_new2_run2_reco2_all" "nl_run2b_rhc_ext_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prod_numi_swizzle_inclusive_v4b_run2_beamon_run2_reco2_run2_beam_good_reco2" "nl_run2b_rhc_data_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_LYAttenuation_run2_RHC_reco2_v08_00_00_55_run2_reco2" "nl_run2b_rhc_detvar_ly_attenuation_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_400k_run2_CV_reco2_detvar_LYRayleigh_run2_reco2" "nl_run2b_rhc_detvar_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_400k_run2_CV_reco2_detvar_LYDown_run2_reco2" "nl_run2b_rhc_detvar_ly_down_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_run2_RHC_standard_nu_reco2_detvar_SCE_run2_reco2" "nl_run2b_rhc_detvar_sce_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_Recomb2_run2_RHC_reco2_v08_00_00_55_run2_reco2" "nl_run2b_rhc_detvar_recomb2_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModX_run2_RHC_reco2_v08_00_00_55_run2_reco2" "nl_run2b_rhc_detvar_wiremodx_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModYZ_run2_RHC_reco2_v08_00_00_55_run2_reco2" "nl_run2b_rhc_detvar_wiremodyz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModThetaXZ_run2_RHC_reco2_v08_00_00_55_run2_reco2" "nl_run2b_rhc_detvar_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_wsplines_run2_RHC_reco2_v08_00_00_5_run2_reco2" "nl_run2b_rhc_detvar_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_overlay_detvar_LYAttenuation_run2_RHC_nue_reco2_v08_00_00_55_run2_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_ly_attenuation_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_LYRayleigh_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_LYDown_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_ly_down_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_SCE_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_sce_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_Recomb2_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_recomb2_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_WireModX_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_wiremodx_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_WireModYZ_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_wiremodyz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_detvar_WireModThetaXZ_run2_RHC_reco2_v08_00_00_55_run2_reco2_goodruns" "nl_run2b_rhc_detvar_nue_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_pandora_reco2_run2_RHC_WireModThetaYZ_withSplines_reco2_goodruns" "nl_run2b_rhc_detvar_nue_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run2_hardcoded"

# Run 3b RHC
"${apply_script}" "prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_all_snapshot" "nl_run3b_rhc_beam_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_intrinsic_nue_overlay_Reco2_CV_run3b_v08_00_00_48_reco2_run3b_reco2" "nl_run3b_rhc_nue_cv_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2" "nl_run3b_rhc_ext_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prod_numi_mcc9_v08_00_00_45_run3b_run3b_reco2_beam_good_reco2" "nl_run3b_rhc_data_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "run3_numi_rhc_CCNCPi0_overlay_v08_00_00_66_reco2_fv_fix_reco2_reco2" "nl_run3b_rhc_ccncpi0_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_cv_reco2_v08_00_00_54_run3_reco2" "nl_run3b_rhc_detvar_cv_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_LYattenuation_reco2_run3_reco2" "nl_run3b_rhc_detvar_ly_attenuation_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_LYrayleigh_reco2_run3_reco2" "nl_run3b_rhc_detvar_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_LYDown_reco2_v08_00_00_54_run3_reco2" "nl_run3b_rhc_detvar_ly_down_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_SCE_reco2_v08_00_00_54_run3_reco2" "nl_run3b_rhc_detvar_sce_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_Recomb2_reco2_v08_00_00_54_run3_reco2" "nl_run3b_rhc_detvar_recomb2_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_cv_v08_00_00_54_run3_reco2_WiremodX_run3b_reco2" "nl_run3b_rhc_detvar_wiremodx_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_WireModYZ_reco2_run3b_reco2" "nl_run3b_rhc_detvar_wiremodyz_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_WireModThetaXZ_reco2_v08_00_00_54_run3_reco2" "nl_run3b_rhc_detvar_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "numi_run3_standard_nu_overlay_WireModThetaYZwithsplines_reco2_run3_reco2" "nl_run3b_rhc_detvar_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_intrinsic_nue_mcc9_v08_00_00_48_LYRayleigh_run3b_run3b_reco2" "nl_run3b_rhc_detvar_nue_ly_rayleigh_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_LYDown_reco2_run3b_reco2" "nl_run3b_rhc_detvar_nue_ly_down_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_run3b_reco2_LY_Attenuation_intrinsic_nue_run_3b_reco2" "nl_run3b_rhc_detvar_nue_ly_attenuation_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_Reco2_SCE_run3b_v08_00_00_48_reco2_run3b_reco2" "nl_run3b_rhc_detvar_nue_sce_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_Recomb2_reco1_run3b_reco2_run3b_reco2" "nl_run3b_rhc_detvar_nue_recomb2_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_standard_intrinsic_nue_mcc9_v08_00_00_48_WireModX_reco2_run3b_run3b_reco2" "nl_run3b_rhc_detvar_nue_wiremodx_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_WireModYZ_reco2_run3b_reco2" "nl_run3b_rhc_detvar_nue_wiremodyz_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_nue_overlay_mcc9_v08_00_00_48_WireModThetaXZ_reco2_run3b_reco2" "nl_run3b_rhc_detvar_nue_wiremodthetaxz_goodruns" "defname: goodruns_mcc9_run3_hardcoded"
"${apply_script}" "prodgenie_numi_run3b_reco2_WireModThetaYZ_sigmaSplines_fixed2_intrinsic_nue_run_3b_reco2" "nl_run3b_rhc_detvar_nue_wiremodthetayz_goodruns" "defname: goodruns_mcc9_run3_hardcoded"

# Run 4a
"${apply_script}" "run4a_NuMI_RHC_nu_overlay_pandora_unified_reco2_reco2_reco2" "nl_run4a_rhc_beam_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "numi_run4a_rhc_nue_overlay_pandora_unified_reco2_run4a_rhc_reco2" "nl_run4a_rhc_nue_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "run_4a_numi_rhc_dirt_overlay_pandora_unified_reco2_run4a_rhc_reco2" "nl_run4a_rhc_dirt_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "NuMI_beam_off_reco_2_for_run_4a_run4a_reco2" "nl_run4a_rhc_ext_goodruns" "defname: goodruns_mcc9_run4_hardcoded"

# Run 4b
"${apply_script}" "numi_run4b_rhc_nu_overlay_pandora_unified_reco2_run4b_reco2" "nl_run4b_rhc_beam_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "numi_run4b_rhc_nue_overlay_pandora_unified_reco2_run4b_rhc_reco2" "nl_run4b_rhc_nue_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "numi_run4b_rhc_dirt_overlay_pandora_unified_reco2_run4b_reco2" "nl_run4b_rhc_dirt_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "Run_4b_NuMI_beam_off_Pandora_Reco2_run4b_reco2_all" "nl_run4b_rhc_ext_pandora_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "Run_4b_NuMI_beam_on_Pandora_Reco2_run4b_reco2_beam_good" "nl_run4b_rhc_data_pandora_goodruns" "defname: goodruns_mcc9_run4_hardcoded"

# Run 4c
"${apply_script}" "numi_run4c_fhc_nu_overlay_pandora_unified_reco2_run4c_reco2" "nl_run4c_fhc_beam_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "numi_run4c_fhc_nue_overlay_pandora_unified_reco2_run4c_fhc_reco2" "nl_run4c_fhc_nue_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "numi_run4c_fhc_dirt_overlay_pandora_unified_reco2_run4c_reco2" "nl_run4c_fhc_dirt_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "run4c_numi_beam_off_pandora_reco2_run4c_reco2_all" "nl_run4c_fhc_ext_pandora_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "Run_4c_NuMI_beam_on_Pandora_Reco2_run4c_reco2_beam_good" "nl_run4c_fhc_data_pandora_goodruns" "defname: goodruns_mcc9_run4_hardcoded"

# Run 4d
"${apply_script}" "run4d_numi_fhc_nu_overlay_pandora_unified_reco2_run4d_fhc_reco2" "nl_run4d_fhc_beam_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "run4d_numi_fhc_nue_overlay_pandora_unified_reco2_run4d_fhc_reco2" "nl_run4d_fhc_nue_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "numi_run4d_fhc_dirt_overlay_pandora_unified_reco2_run4d_reco2" "nl_run4d_fhc_dirt_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "run4d_numi_beam_off_pandora_reco2_run4d_reco2_all" "nl_run4d_fhc_ext_pandora_goodruns" "defname: goodruns_mcc9_run4_hardcoded"
"${apply_script}" "run4d_NuMI_beam_on_pandora_unified_reco2_run4d_reco2_beam_good" "nl_run4d_fhc_data_pandora_unified_goodruns" "defname: goodruns_mcc9_run4_hardcoded"

# Run 5
"${apply_script}" "NuMI_Run5_FHC_nu_overlay_pandora_Reco2_run5_reco2_nonzerolifetime_goodruns" "nl_run5_fhc_beam_nonzerolifetime_goodruns" "defname: goodruns_mcc9_run5_hardcoded"
"${apply_script}" "run5_numi_fhc_nue_overlay_pandora_unified_reco2_run5_fhc_reco2_goodruns" "nl_run5_fhc_nue_pandora_unified_goodruns" "defname: goodruns_mcc9_run5_hardcoded"
"${apply_script}" "numi_fhc_dirt_overlay_pandora_reco2_run5_run5_fhc_reco2" "nl_run5_fhc_dirt_pandora_goodruns" "defname: goodruns_mcc9_run5_hardcoded"
"${apply_script}" "run5_numi_fhc_CCNCPi0_overlay_v08_00_00_66_reco2_fv_fix_reco2_reco2" "nl_run5_fhc_ccncpi0_goodruns" "defname: goodruns_mcc9_run5_hardcoded"
"${apply_script}" "run5_numi_beam_off_pandora_reco2_all_run5_reco2_all" "nl_run5_fhc_ext_pandora_goodruns" "defname: goodruns_mcc9_run5_hardcoded"
"${apply_script}" "run5_numi_beam_on_pandora_reco2_all_run5_reco2_beam_good" "nl_run5_fhc_data_pandora_goodruns" "defname: goodruns_mcc9_run5_hardcoded"
