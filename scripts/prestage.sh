#!/bin/bash

LOG_DIR="./tmp"
mkdir -p "${LOG_DIR}"

prestage_dataset_def() {
    local defname="$1"
    local log_file="${LOG_DIR}/prestage_${defname}.log"

    echo "Prestaging dataset: ${defname}. Logging to ${log_file}"
    nohup samweb prestage-dataset --defname="${defname}" >& "${log_file}" &
    sleep 5
}

STRANGE_SAMPLES=(
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
)

BEAM_SAMPLES=(
    "New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2"
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

echo "Starting dataset prestaging..."

case "$1" in
    "strange")
        echo "Prestaging only strange samples..."
        for defname in "${STRANGE_SAMPLES[@]}"; do
            prestage_dataset_def "${defname}"
        done
        ;;
    "beam")
        echo "Prestaging only beam samples..."
        for defname in "${BEAM_SAMPLES[@]}"; do
            prestage_dataset_def "${defname}"
        done
        ;;
    *)
        echo "No specific sample type provided. Prestaging both strange and beam samples."
        echo "Prestaging strange samples first..."
        for defname in "${STRANGE_SAMPLES[@]}"; do
            prestage_dataset_def "${defname}"
        done
        echo ""
        echo "Prestaging beam samples now..."
        for defname in "${BEAM_SAMPLES[@]}"; do
            prestage_dataset_def "${defname}"
        done
        ;;
esac

echo "All selected prestage commands submitted. Check logs in ${LOG_DIR} for status."