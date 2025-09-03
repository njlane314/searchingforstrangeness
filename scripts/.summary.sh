#!/bin/bash

input_definitions=(
    "New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2"
    "prod_strange_resample_fhc_run2_fhc_reco2_reco2"
    "detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2"
    "Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2"
    "Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2"
    "detvar_prod_strange_resample_fhc_run_respin_wiremodX_sce_reco2_reco2"
    "detvar_prod_strange_resample_fhc_run_respin_wiremodYZ_sce_reco2_reco2"
    "Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2"
    "Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2"
    "detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2"
    "detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2"
    "prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_LY_suppression75attenuation8m_run1_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_LY_Rayleigh_run1_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_LYDown_run1_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_v08_00_00_53_SCE_300k_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_Recomb2_run1_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_WireModX_run1_reco2_fixed_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_WireModYZ_run1_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_v08_00_00_53_WireModThetaXZ_300k_reco2_run1_reco2"
    "prodgenie_numi_nu_overlay_detvar_WireModThetaYZ_withSplines_run1_reco2_run1_reco2"
)

echo "--- Samweb List Summary Definitions ---"
echo ""

for def in "${input_definitions[@]}"; do
    echo "Running for definition: $def"
    samweb list-files --summary "defname: $def"
    echo ""
done

echo "--- Script Finished ---"