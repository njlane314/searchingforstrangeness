#!/bin/bash

LOG_DIR=".tmp"
FILE_LIMIT=1000

PARENT_DEFINITIONS=(
    #"detvar_prod_strange_resample_fhc_run1_respin_cv_reco2_reco2"
    #"Run_1_MuMI_FHC_detvars_LY_Rayleigh_reco2_reco2_reco2"
    #"Run1_NuMI_FHC_detvars_LY_Down_Reco2_lydown_reco2"
    #"detvar_prod_strange_resample_fhc_run_respin_wiremodX_sce_reco2_reco2"
    #"detvar_prod_strange_resample_fhc_run_respin_wiremodYZ_sce_reco2_reco2"
    #"Run1_NuMI_nu_overlay_FHC_Strangeness_DetVar_WireMod_YZ_reco2_reco2_reco2"
    #"Run1_NuMI_FHC_detvars_wiremod_thetaYZ_Reco2_reco2_reco2"
    #"detvar_prod_strange_resample_fhc_run1_respin_sce_reco2_reco2"
    #"detvar_prod_strange_resample_fhc_run1_respin_recomb2_reco2_reco2"
    "prodgenie_numi_nu_overlay_v08_00_00_53_CV_300k_reco2_run1_reco2"
)

mkdir -p "${LOG_DIR}"

echo "Starting process to create and prestage smaller definitions."
echo "Each new definition will contain ${FILE_LIMIT} files."
echo "---"

for PARENT_DEF in "${PARENT_DEFINITIONS[@]}"; do
    CHILD_DEF="${PARENT_DEF}_${FILE_LIMIT}"
    
    echo "Processing Parent: ${PARENT_DEF}"

    samweb describe-definition "${CHILD_DEF}" &> /dev/null
    if [ $? -ne 0 ]; then
        echo "-> Definition '${CHILD_DEF}' not found. Creating..."
        samweb create-definition "${CHILD_DEF}" "defname: ${PARENT_DEF} with limit ${FILE_LIMIT}"
        if [ $? -ne 0 ]; then
            echo "   Error: Failed to create definition. Skipping prestage."
            echo "---"
            continue
        fi
        echo "   Successfully created definition."
    else
        echo "-> Definition '${CHILD_DEF}' already exists. Skipping creation."
    fi

    LOG_FILE="${LOG_DIR}/prestage_${CHILD_DEF}.log"
    echo "-> Starting prestage for '${CHILD_DEF}'."
    echo "   Logging to: ${LOG_FILE}"
    nohup samweb prestage-dataset --defname="${CHILD_DEF}" >& "${LOG_FILE}" &
    
    echo "---"
done

echo "All jobs have been launched in the background."
echo "Monitor progress by checking the log files in the '${LOG_DIR}' directory."