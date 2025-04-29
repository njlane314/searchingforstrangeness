#include "FilePropertiesManager.h"
#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TTree.h>

int main() {
    try {
        FilePropertiesManager& fpm = FilePropertiesManager::Instance();
        std::vector<std::string> input_files;

        // Extract files directly from ntuple_file_map_
        for (const auto& run_pair : fpm.ntuple_file_map()) {
            for (const auto& type_pair : run_pair.second) {
                for (const auto& file : type_pair.second) {
                    input_files.push_back(file);
                }
            }
        }

        // Process each file
        for (const auto& f_name : input_files) {
            try {
                NtupleFileType file_type = fpm.get_ntuple_file_type(f_name);
                std::string type_str = fpm.ntuple_type_to_string(file_type);

                TFile in_file(f_name.c_str(), "READ");
                if (!in_file.IsOpen()) {
                    std::cerr << "Error: Could not open file " << f_name << "\n";
                    continue;
                }

                TTree* tree = dynamic_cast<TTree*>(in_file.Get("compressed_tree"));
                if (!tree) {
                    std::cerr << "Error: Missing 'compressed_tree' in file " << f_name << "\n";
                    in_file.Close();
                    continue;
                }

                int nhits_u;
                float charge_u;
                float bdt_score;
                double weight;

                tree->SetBranchAddress("nhits_u", &nhits_u);
                tree->SetBranchAddress("charge_u", &charge_u);
                tree->SetBranchAddress("bdt_score", &bdt_score);
                tree->SetBranchAddress("weight", &weight);

                Long64_t nentries = tree->GetEntries();
                std::cout << "File: " << f_name << ", Type: " << type_str << ", Entries: " << nentries << "\n";

                in_file.Close();
            } catch (const std::runtime_error& e) {
                std::cerr << "Error processing file " << f_name << ": " << e.what() << "\n";
            }
        }
    } catch (const std::runtime_error& e) {
        std::cerr << "Error initializing FilePropertiesManager: " << e.what() << "\n";
        return 1;
    }

    return 0;
}