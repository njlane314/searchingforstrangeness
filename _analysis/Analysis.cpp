
// --- Example: main.cpp (Illustrative Usage) ---
/*
#include "DataLoader.h"
#include <iostream>
// #include "TApplication.h" // For interactive ROOT graphics
// #include "TCanvas.h"    // For drawing

namespace AnalysisFramework {
    void PopulateRunConfigurations(ConfigurationManager& config_mgr, const std::string& config_set_name = "default");
}

int main(int argc, char* argv[]) {
    // ROOT::EnableImplicitMT(); 
    AnalysisFramework::FilePathConfig file_paths;
    // file_paths.ntuple_path_base = "/path/to/your/ntuples/";

    AnalysisFramework::ConfigurationManager config_mgr(file_paths);
    AnalysisFramework::VariableManager var_mgr;

    try {
        config_mgr.PopulateConfigurations("default");
    } catch (const std::exception& e) { std::cerr << "FATAL: Config population error: " << e.what() << std::endl; return 1;}

    std::vector<std::string> runs_to_load = {"run1_bnb_official"};
    
    AnalysisFramework::DataLoader::LoadedRunData loaded_collections;
    try {
        loaded_collections = AnalysisFramework::DataLoader::LoadRuns(
            config_mgr, var_mgr, runs_to_load,
            "data_%B%", "ext_%B%", "mc_%B%", {"dirt_%B%"}, 
            false, // blinded
            true,  // load_dedicated_intrinsic_strangeness_sample
            "intrinsic_strangeness_%B%", 
            true,  // load_reco_event_level
            true,  // load_reco_track_level 
            true,  // load_mc_truth_level
            true,  // load_mc_weights_and_systematics
            true   // apply_specific_is_event_weights (for "is_channel_correction_weight" branch)
        );
    } catch (const std::exception& e) { std::cerr << "FATAL: Data loading error: " << e.what() << std::endl; return 1;}

    std::cout << "Loaded data collections. Total Data POT: " << loaded_collections.pot << std::endl;

    std::string mc_bnb_key = "mc_bnb"; 
    if (loaded_collections.dfs.count(mc_bnb_key)) {
        auto& mc_df = loaded_collections.dfs.at(mc_bnb_key);
        if(mc_df.GetNRuns() > 0) { 
            std::cout << "Inclusive BNB MC (" << mc_bnb_key << ") entries: " << *mc_df.Count() << std::endl;
            // auto h_cat = mc_df.Histo1D<int>({"h_event_category_mc_bnb", "Event Categories (Inclusive BNB MC);Category ID;Events", 250, 0, 250}, "event_category", "final_event_weight");
            // if(h_cat->GetEntries() > 0) std::cout << "  Mean event_category in " << mc_bnb_key << ": " << h_cat->GetMean() << std::endl;
        }
    }
    return 0;
}
*/
