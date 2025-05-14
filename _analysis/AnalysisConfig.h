#include "DataLoader.h" // Assuming this header is named DataLoader.h

namespace AnalysisFramework {

// This function would be called by your main analysis executable
// to populate the ConfigurationManager.
void PopulateRunConfigurations(ConfigurationManager& config_mgr) {

    // --- BNB Run 1 Example ---
    RunConfig bnb_r1;
    bnb_r1.run_id_key = "run1_bnb_uboone"; // Unique identifier
    bnb_r1.beam_type = "BNB";
    bnb_r1.relative_path = "Run1/BNB/prod_v08_00_00_XX"; // Adjust to your actual path structure

    bnb_r1.datasets["data_bnb"] = {
        "data_bnb_run1_reco2",      // file_prefix
        1.0e20,                     // pot
        5000000,                    // triggers
        "",                         // filter_query (N/A for data)
        NtupleFileType::kDataBNB,
        "data_bnb"                  // name_in_code (key for maps)
    };
    bnb_r1.datasets["ext_bnb"]  = {
        "ext_bnb_run1_reco2", 0.0, 20000000, "", NtupleFileType::kExtBNB, "ext_bnb"
    };
    bnb_r1.datasets["mc_bnb"]   = { // Inclusive BNB MC
        "mc_bnb_run1_overlay_reco2", 6.0e20, 0, "", NtupleFileType::kMCBNB, "mc_bnb"
    };
    bnb_r1.datasets["dirt_bnb"] = {
        "mc_bnb_run1_dirt_reco2", 12.0e20, 0, 
        "isTrueDirtInteraction == true", // Example filter query for dirt sample
        NtupleFileType::kMCDirtBNB, "dirt_bnb"
    };
    // Dedicated Intrinsic Strangeness sample for BNB Run 1
    bnb_r1.datasets["intrinsic_strangeness_bnb"] = {
        "mc_bnb_run1_intrinsic_strangeness_reco2", // File prefix for IS sample
        2.0e21, // POT for this specific sample
        0,
        "truth_HasKaonPlus == 1 || truth_HasSigmaZero == 1", // Example filter defining IS events
        NtupleFileType::kMCIntrinsicStrangenessBNB,
        "intrinsic_strangeness_bnb"
    };
    // Example of a truth-filtered nue sample from the main BNB MC (if not using a dedicated IS sample for this)
    // bnb_r1.datasets["nue_filter_bnb"] = {
    //     "mc_bnb_run1_overlay_reco2", // Uses same file as mc_bnb
    //     6.0e20, // Same POT as mc_bnb
    //     0,
    //     "(abs(nu_pdg) == 12 && ccnc == 0)", // Filter for true nue CC events
    //     NtupleFileType::kMCBNB, // Still of type MCBNB, but filtered
    //     "nue_filter_bnb"
    // };
    config_mgr.AddRunConfig(bnb_r1);


    // --- NuMI Run 5 Example ---
    RunConfig numi_r5;
    numi_r5.run_id_key = "run5a_numi_fhc";
    numi_r5.beam_type = "NuMI";
    numi_r5.relative_path = "Run5/NuMI/FHC/prod_v09_00_00_YY";

    numi_r5.datasets["data_numi"] = {
        "data_numi_run5a_fhc", 5.0e20, 12000000, "", NtupleFileType::kDataNuMI, "data_numi"
    };
    numi_r5.datasets["ext_numi"] = {
         "ext_numi_run5a_fhc", 0.0, 30000000, "", NtupleFileType::kExtNuMI, "ext_numi"
    };
    numi_r5.datasets["mc_numi"] = {
        "mc_numi_run5a_fhc_overlay", 25.0e20, 0, "", NtupleFileType::kMCNuMI, "mc_numi"
    };
    // Dedicated Intrinsic Strangeness sample for NuMI Run 5
    numi_r5.datasets["intrinsic_strangeness_numi"] = {
        "mc_numi_run5a_fhc_intrinsic_strangeness", 
        50.0e20, 
        0, 
        "truth_IsKaonEventFromNu == 1", // Example truth filter for this IS sample
        NtupleFileType::kMCIntrinsicStrangenessNuMI,
        "intrinsic_strangeness_numi"
    };
    config_mgr.AddRunConfig(numi_r5);

    // ... Add more RunConfig objects for other runs/beams ...

    // --- Detector Variation Example for BNB Run 1 ---
    // DetVarRunConfig dv_bnb_r1;
    // dv_bnb_r1.run_id_key = "run1_bnb_uboone"; // Must match a RunConfig key
    // dv_bnb_r1.beam_type = "BNB";
    // dv_bnb_r1.relative_path_detvar_base = "Run1/BNB/prod_v08_00_00_XX_detvar"; // Or however it's structured
    //
    // DetVarDatasetInfo lydown_info;
    // lydown_info.file_prefix_variant_suffix = "_lydown"; // Appended to CV file prefix
    // lydown_info.base_sample_type = NtupleFileType::kMCBNB; // It's a variation of the BNB MC
    // lydown_info.name_in_code = "mc_bnb_lydown"; // Unique key for this varied sample
    // dv_bnb_r1.variations["lydown"]["mc_bnb"] = lydown_info; // "mc_bnb" is the CV dataset it varies
    //
    // DetVarDatasetInfo sce_info;
    // sce_info.file_prefix_variant_suffix = "_sce";
    // sce_info.base_sample_type = NtupleFileType::kMCBNB;
    // sce_info.name_in_code = "mc_bnb_sce";
    // dv_bnb_r1.variations["sce"]["mc_bnb"] = sce_info;
    //
    // config_mgr.AddDetVarRunConfig(dv_bnb_r1);
}

} // namespace AnalysisFramework



// --- Example: AnalysisConfig.cpp (Separate File) ---
/*
#include "DataLoader.h" // Assuming this header is named DataLoader.h

namespace AnalysisFramework {
void ConfigurationManager::PopulateConfigurations(const std::string& config_set_name) {
    if (config_set_name == "default") {
        RunConfig bnb_r1;
        bnb_r1.run_id_key = "run1_bnb_official"; bnb_r1.beam_type = "BNB"; bnb_r1.relative_path = "Run1/BNB_processed_v1";
        bnb_r1.datasets["data_bnb"] = {"data_bnb_run1", 6.6e20, 15000000, "", NTupleFileType::kDataBNB, "data_bnb"};
        bnb_r1.datasets["ext_bnb"]  = {"ext_bnb_run1", 0.0, 60000000, "", NTupleFileType::kBackgroundEXT_BNB, "ext_bnb"};
        bnb_r1.datasets["mc_bnb"]   = {"mc_bnb_overlay_run1", 20.0e20, 0, "", NTupleFileType::kMCInclusiveBNB, "mc_bnb"};
        bnb_r1.datasets["dirt_bnb"] = {"mc_bnb_dirt_run1", 40.0e20, 0, "isTrueDirtEvent == 1", NTupleFileType::kMCDirtBNB, "dirt_bnb"};
        // Example: Dedicated IS sample with a specific filter defining its signal content
        bnb_r1.datasets["intrinsic_strangeness_bnb"] = { 
            "mc_bnb_is_run1", // File prefix
            100.0e20,         // POT
            0,                // Triggers
            "true_event_has_kaon_plus == 1 && nu_pdg == 14 && ccnc == 0", // Filter to define IS signal within this file
            NTupleFileType::kMCStrangenessSignalBNB, 
            "intrinsic_strangeness_bnb"
        };
        this->AddRunConfig(bnb_r1);
        // Add other RunConfig and DetVarRunConfig objects here
    }
}
}
*/