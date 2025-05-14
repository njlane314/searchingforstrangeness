#pragma once

#include <string>
#include <vector>
#include <stdexcept> // For std::runtime_error
#include <set>       // For std::set to handle unique column names
#include <algorithm> // For std::sort, std::unique

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx" // For RDFUtils::GetColumnNames

namespace AnalysisFramework {

class MLDatasetSaver {
public:
    MLDatasetSaver() = default;

    /**
     * @brief Saves the contents of an RDataFrame to a new ROOT file and TTree.
     * @param df The RDataFrame to save. (Passed by value as Snapshot is terminal)
     * @param output_file_name The name of the ROOT file to create.
     * @param tree_name The name of the TTree within the ROOT file.
     * @param columns_to_save User-specified list of columns to save. If empty, behavior depends on force_keep_columns.
     * @param force_keep_columns Columns to always include if columns_to_save is specified. Ignored if columns_to_save is empty (then all are saved).
     * @param snapshot_options Optional RSnapshotOptions.
     * @return True if successful (Snapshot throws on error, so this is more for conceptual flow).
     */
    bool SaveDataFrame(
        ROOT::RDF::RInterface df,
        const std::string& output_file_name,
        const std::string& tree_name,
        const std::vector<std::string>& columns_to_save = {},
        const std::vector<std::string>& force_keep_columns = {},
        const ROOT::RDF::RSnapshotOptions& snapshot_options = {}
    ) const {
        if (output_file_name.empty() || tree_name.empty()) {
            throw std::invalid_argument("Output file name and tree name cannot be empty.");
        }

        ROOT::Detail::RDF::ColumnNameVector_t final_columns_to_save_rdf;
        bool save_all_columns = true;

        if (!columns_to_save.empty()) {
            save_all_columns = false;
            std::set<std::string> unique_cols(columns_to_save.begin(), columns_to_save.end());
            unique_cols.insert(force_keep_columns.begin(), force_keep_columns.end());
            for (const auto& col_name : unique_cols) {
                final_columns_to_save_rdf.emplace_back(col_name);
            }
        } else if (!force_keep_columns.empty()) { 
            // If columns_to_save is empty, but force_keep_columns is not, save only force_keep_columns
            save_all_columns = false;
             for (const auto& col_name : force_keep_columns) {
                final_columns_to_save_rdf.emplace_back(col_name);
            }
        }
        // If both columns_to_save and force_keep_columns are empty, save_all_columns remains true.


        try {
            if (save_all_columns) {
                df.Snapshot(tree_name, output_file_name, ROOT::Detail::RDF::ColumnNameVector_t{}, snapshot_options);
            } else {
                if (final_columns_to_save_rdf.empty() && !columns_to_save.empty()){
                    // This case implies columns_to_save was not empty but force_keep_columns might have cleared it (not possible with current set logic)
                    // Or, more likely, user provided empty columns_to_save and empty force_keep_columns, handled by save_all_columns=true
                    // However, if somehow final_columns_to_save_rdf is empty when save_all_columns is false, it's an issue.
                    // For safety, if user specified columns but the final list is empty, it might be better to save all or error.
                    // RDataFrame's Snapshot with an empty explicit list might behave like saving all, or error.
                    // Let's assume if final_columns_to_save_rdf is empty here, it means save all (as per initial logic).
                    // This path should ideally not be hit if save_all_columns logic is correct.
                     df.Snapshot(tree_name, output_file_name, ROOT::Detail::RDF::ColumnNameVector_t{}, snapshot_options);
                } else {
                    df.Snapshot(tree_name, output_file_name, final_columns_to_save_rdf, snapshot_options);
                }
            }
        } catch (const std::exception& e) {
            throw; // Re-throw to let the caller handle it
        }
        return true;
    }

    /**
     * @brief Saves multiple RDataFrames to the same ROOT file, each in its own TTree.
     * @param data_frames Map of TTree name to RDataFrame.
     * @param output_file_name Name of the ROOT file.
     * @param columns_to_save_map Map of TTree name to specific columns for that TTree.
     * @param default_force_keep_columns Default list of columns to force-keep for each TTree if its specific list in columns_to_save_map is not empty.
     * @param snapshot_options Options for snapshotting.
     */
    bool SaveMultipleDataFrames(
        const std::map<std::string, ROOT::RDF::RInterface>& data_frames,
        const std::string& output_file_name,
        const std::map<std::string, std::vector<std::string>>& columns_to_save_map = {},
        const std::vector<std::string>& default_force_keep_columns = {},
        ROOT::RDF::RSnapshotOptions snapshot_options = {}
    ) const {
        bool first_tree = true;
        for (const auto& pair : data_frames) {
            const std::string& tree_name = pair.first;
            ROOT::RDF::RInterface df = pair.second; // Make a copy for Snapshot

            std::vector<std::string> specific_cols_for_this_tree;
            auto it_cols_map = columns_to_save_map.find(tree_name);
            if (it_cols_map != columns_to_save_map.end()) {
                specific_cols_for_this_tree = it_cols_map->second;
            }

            std::vector<std::string> current_force_keep_cols = default_force_keep_columns;
            // If specific_cols_for_this_tree is empty, force_keep_columns in SaveDataFrame will be default_force_keep_columns.
            // If specific_cols_for_this_tree is NOT empty, default_force_keep_columns will be merged into it by SaveDataFrame.

            ROOT::RDF::RSnapshotOptions current_options = snapshot_options;
            if (!first_tree) {
                current_options.fMode = "UPDATE";
            }

            // Pass specific_cols_for_this_tree as columns_to_save,
            // and default_force_keep_columns as the force_keep_columns to SaveDataFrame.
            SaveDataFrame(df, output_file_name, tree_name, specific_cols_for_this_tree, default_force_keep_columns, current_options);
            first_tree = false;
        }
        return true;
    }
};

} // namespace AnalysisFramework


// --- Example Usage (to be placed in your main.cpp or analysis driver) ---
/*
#include "DataLoader.h"     // Your main data loading header
#include "MLDatasetSaver.h" // This new header
#include "Selections.h"     // For selection strings
#include <iostream>

namespace AnalysisFramework {
    void PopulateRunConfigurations(ConfigurationManager& config_mgr, const std::string& config_set_name = "default");
}

int main() {
    AnalysisFramework::FilePathConfig file_paths;
    AnalysisFramework::ConfigurationManager config_mgr(file_paths);
    AnalysisFramework::VariableManager var_mgr;

    try {
        config_mgr.PopulateRunConfigurations("default");
    } catch (const std::exception& e) {
        std::cerr << "Config population error: " << e.what() << std::endl; return 1;
    }

    std::vector<std::string> runs_to_load = {"run1_bnb_official"};
    
    AnalysisFramework::DataLoader::LoadedRunData loaded_collections;
    try {
        loaded_collections = AnalysisFramework::DataLoader::LoadRuns(
            config_mgr, var_mgr, runs_to_load,
            "data_%B%", "ext_%B%", "mc_%B%", {"dirt_%B%"}, 
            true,  // blinded = true 
            true,  // load_dedicated_intrinsic_strangeness_sample
            "intrinsic_strangeness_%B%", 
            true, true, true, true, true 
        );
    } catch (const std::exception& e) {
        std::cerr << "Data loading error: " << e.what() << std::endl; return 1;
    }

    AnalysisFramework::MLDatasetSaver saver;

    // Define essential MC weight and truth columns you always want to save for ML
    std::vector<std::string> essential_mc_columns_for_ml = {
        "final_event_weight", // Already POT scaled
        "base_event_weight",  // Before POT scaling, might include specific IS weight
        "event_category",     // Truth category for labeling
        // Add any systematic knob branches you might want to use as features or for reweighting studies
        "knobRPAup", "knobRPAdown", 
        // Add multi-universe weight branches if you plan to use them directly in ML or for uncertainty studies
        // "weightsGenie", // Note: these are ROOT::RVec<float> or similar, ensure ML framework can handle
        // Truth information for context or alternative labeling
        "nu_pdg", "nu_e", "total_true_strange_particles", "inclusive_strangeness_multiplicity_type"
    };


    // Example 1: Save selected events from inclusive MC for ML training
    std::string mc_bnb_key = "mc_bnb"; 
    if (loaded_collections.dfs.count(mc_bnb_key)) {
        auto mc_df_orig = loaded_collections.dfs.at(mc_bnb_key);
        
        std::string selection_for_ml = AnalysisFramework::Selections::NuMuPreselectionCRT + 
                                       " && selected_muon_idx != -1"; 
                                       // Keep selection broader, event_category can be used for ML target
        
        auto df_for_ml = mc_df_orig.Filter(selection_for_ml, "ML Inclusive MC Preselection");

        // User-defined features for the model
        std::vector<std::string> ml_features = {
            "selected_muon_length", "selected_muon_momentum_range", "selected_muon_cos_theta",
            "n_reco_tracks_passing_basic_cuts"
            // "event_category" is already in essential_mc_columns_for_ml
        };
        
        std::cout << "Attempting to save inclusive MC for ML..." << std::endl;
        if (*df_for_ml.Count() > 0) {
            // SaveDataFrame will merge ml_features and essential_mc_columns_for_ml
            saver.SaveDataFrame(df_for_ml, "ml_inclusive_mc_dataset.root", "ml_tree", ml_features, essential_mc_columns_for_ml);
            std::cout << "Inclusive MC dataset for ML saved." << std::endl;
        } else {
            std::cout << "No events passed selection for inclusive MC ML dataset." << std::endl;
        }
    }

    // Example 2: Save dedicated Intrinsic Strangeness MC events, saving ALL columns
    std::string is_bnb_key = "intrinsic_strangeness_bnb";
    if (loaded_collections.dfs.count(is_bnb_key)) {
        auto is_df_orig = loaded_collections.dfs.at(is_bnb_key);
        
        std::string is_selection_for_ml = AnalysisFramework::Selections::NuMuPreselectionCRT + 
                                          " && selected_muon_idx != -1";

        auto df_is_for_ml = is_df_orig.Filter(is_selection_for_ml, "ML IS MC Preselection");

        std::cout << "Attempting to save IS MC for ML (all columns)..." << std::endl;
        if (*df_is_for_ml.Count() > 0) {
            // Pass empty columns_to_save and empty force_keep_columns to save all columns
            saver.SaveDataFrame(df_is_for_ml, "ml_intrinsic_strangeness_dataset_all_cols.root", "is_ml_tree");
            std::cout << "Intrinsic Strangeness MC (all columns) for ML saved." << std::endl;
        } else {
            std::cout << "No events passed selection for IS MC ML dataset." << std::endl;
        }
    }
    
    // Example 3: Save multiple RDataFrames (e.g. different event categories from one sample)
    // if (loaded_collections.dfs.count(mc_bnb_key)) {
    //    auto base_df = loaded_collections.dfs.at(mc_bnb_key);
    //    std::map<std::string, ROOT::RDF::RInterface> dfs_to_save;
    //    std::map<std::string, std::vector<std::string>> cols_map;

    //    std::vector<std::string> features_for_ml = {"selected_muon_momentum_range", "selected_muon_cos_theta"};

    //    dfs_to_save["numu_cc_single_strange_tree"] = base_df.Filter("event_category == 10");
    //    cols_map["numu_cc_single_strange_tree"] = features_for_ml; // Will also get essential_mc_columns_for_ml

    //    dfs_to_save["numu_cc_0pi0p_bkg_tree"] = base_df.Filter("event_category == 100");
    //    cols_map["numu_cc_0pi0p_bkg_tree"] = features_for_ml;

    //    if (!dfs_to_save.empty()) {
    //        saver.SaveMultipleDataFrames(dfs_to_save, "ml_multi_category_dataset.root", cols_map, essential_mc_columns_for_ml);
    //        std::cout << "Multi-category ML dataset saved." << std::endl;
    //    }
    // }

    return 0;
}
*/
