#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TDirectory.h"

constexpr int BOGUS_INT = 9999;

struct EventID {
    int run = BOGUS_INT;
    int sub = BOGUS_INT;
    int evt = BOGUS_INT;

    bool operator<(const EventID& other) const {
        if (run != other.run) return run < other.run;
        if (sub != other.sub) return sub < other.sub;
        return evt < other.evt;
    }
};

class SelectionAnalyser {
public:
    SelectionAnalyser(const std::vector<std::string>& in_file_names, const std::string& output_filename)
        : in_file_names_(in_file_names), output_filename_(output_filename) {}

    void process_selections();

private:
    std::vector<std::string> in_file_names_;
    std::string output_filename_;

    void set_input_branch_addresses(TTree& tree, EventID& id, bool& selected, float& weight);
    void set_output_branch_addresses(TTree& out_tree, EventID& id, bool& pass_fv, bool& pass_image, float& weight, bool create = false);
    float calculate_total_pot(TChain& subruns_ch);
    void check_branch_status(TTree& tree);
};

void SelectionAnalyser::set_input_branch_addresses(TTree& tree, EventID& id, bool& selected, float& weight) {
    tree.SetBranchAddress("run", &id.run);
    tree.SetBranchAddress("sub", &id.sub);
    tree.SetBranchAddress("evt", &id.evt);
    tree.SetBranchAddress("selected", &selected);
    tree.SetBranchAddress("weight", &weight);
}

void SelectionAnalyser::set_output_branch_addresses(TTree& out_tree, EventID& id, bool& pass_fv, bool& pass_image, float& weight, bool create) {
    if (create) {
        out_tree.Branch("run", &id.run, "run/I");
        out_tree.Branch("sub", &id.sub, "sub/I");
        out_tree.Branch("evt", &id.evt, "evt/I");
        out_tree.Branch("pass_fv_selection", &pass_fv, "pass_fv_selection/O");
        out_tree.Branch("pass_image_selection", &pass_image, "pass_image_selection/O");
        out_tree.Branch("weight", &weight, "weight/F");
    } else {
        out_tree.SetBranchAddress("run", &id.run);
        out_tree.SetBranchAddress("sub", &id.sub);
        out_tree.SetBranchAddress("evt", &id.evt);
        out_tree.SetBranchAddress("pass_fv_selection", &pass_fv);
        out_tree.SetBranchAddress("pass_image_selection", &pass_image);
        out_tree.SetBranchAddress("weight", &weight);
    }
}

float SelectionAnalyser::calculate_total_pot(TChain& subruns_ch) {
    float pot;
    subruns_ch.SetBranchAddress("pot", &pot);
    float summed_pot = 0.;
    Long64_t nentries = subruns_ch.GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        subruns_ch.GetEntry(i);
        summed_pot += pot;
    }
    return summed_pot;
}

void SelectionAnalyser::check_branch_status(TTree& tree) {
    if (!tree.GetBranch("run") || !tree.GetBranch("sub") || !tree.GetBranch("evt") || !tree.GetBranch("selected")) {
        std::cerr << "Error: Missing required branches in input tree\n";
        exit(1);
    }
}

void SelectionAnalyser::process_selections() {
    TChain subruns_ch("fiducialFilter/SubRun");
    for (const auto& f_name : in_file_names_) {
        subruns_ch.Add((f_name + "/fiducialFilter/SubRun").c_str());
    }

    float total_pot = 0.;
    if (subruns_ch.GetEntries() > 0) {
        total_pot = calculate_total_pot(subruns_ch);
    }

    std::unique_ptr<TFile> out_file(new TFile(output_filename_.c_str(), "RECREATE"));
    if (!out_file || out_file->IsZombie()) {
        std::cerr << "Error: Could not create output file " << output_filename_ << '\n';
        return;
    }
    out_file->cd();
    TTree out_tree("compressed_tree", "");

    TParameter<float> summed_pot_param("summed_pot", total_pot);
    summed_pot_param.Write();

    TParameter<int> n_input_files("n_input_files", in_file_names_.size());
    n_input_files.Write();

    std::map<EventID, std::tuple<bool, bool, float>> event_flags;

    std::vector<std::pair<std::string, int>> selections = {
        {"fiducialFilter", 0},
        {"imageFilter", 1}
    };

    for (const auto& sel : selections) {
        const std::string& dir_name = sel.first;
        int flag_index = sel.second;

        TChain events_ch((dir_name + "/EventSelectionFilter").c_str());
        for (const auto& f_name : in_file_names_) {
            events_ch.Add(f_name.c_str());
        }

        EventID id;
        bool selected = false;
        float weight = 1.0;
        set_input_branch_addresses(events_ch, id, selected, weight);
        check_branch_status(events_ch);

        Long64_t nentries = events_ch.GetEntries();
        for (Long64_t i = 0; i < nentries; ++i) {
            events_ch.GetEntry(i);
            if (i % 10000 == 0) {
                std::cout << "Processed " << i << " of " << nentries << " entries in " << dir_name << '\n';
            }
            if (flag_index == 0) {
                event_flags[id] = std::make_tuple(selected, std::get<1>(event_flags[id]), weight);
            } else {
                event_flags[id] = std::make_tuple(std::get<0>(event_flags[id]), selected, weight);
            }
        }
    }

    EventID out_id;
    bool pass_fv_selection = false;
    bool pass_image_selection = false;
    float out_weight = 1.0;
    set_output_branch_addresses(out_tree, out_id, pass_fv_selection, pass_image_selection, out_weight, true);

    for (const auto& pair : event_flags) {
        out_id = pair.first;
        pass_fv_selection = std::get<0>(pair.second);
        pass_image_selection = std::get<1>(pair.second);
        out_weight = std::get<2>(pair.second);
        out_tree.Fill();
    }

    out_tree.Write();
    std::cout << "Total events written: " << out_tree.GetEntries() << '\n';
}

void analyse(const std::vector<std::string>& in_file_names, const std::string& output_filename) {
    SelectionAnalyser analyser(in_file_names, output_filename);
    analyser.process_selections();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " OUTPUT_ROOT_FILE INPUT_ROOT_FILE1 [INPUT_ROOT_FILE2 ...]\n";
        return 1;
    }
    std::string output_file(argv[1]);
    std::vector<std::string> input_files;
    for (int i = 2; i < argc; ++i) {
        input_files.push_back(argv[i]);
    }
    analyse(input_files, output_file);
    return 0;
}