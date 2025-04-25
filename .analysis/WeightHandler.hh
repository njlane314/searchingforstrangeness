#pragma once

#include <algorithm>
#include <array>
#include <map>
#include <string>

#include "TTree.h"

#include "TreeUtils.hh"

class WeightHandler {
  public:
    WeightHandler() {};

    void set_branch_addresses(TTree& in_tree, const std::vector<std::string>* branch_names = nullptr);

    void add_branch(TTree& in_tree, const std::string& branch_name, bool throw_when_missing = true);

    void set_branch_addresses(TTree& in_tree, const std::string& branch_name);

    inline const auto& weight_map() const { return weight_map_; }
    inline auto& weight_map() { return weight_map_; }

  protected:
    std::map<std::string, MyPointer<std::vector<double>>> weight_map_;
};

void WeightHandler::set_branch_addresses(TTree& in_tree, const std::vector<std::string>* branch_names) {
    weight_map_.clear();
    auto* lob = in_tree.GetListOfBranches();
    for (int b = 0; b < lob->GetEntries(); ++b) {
        auto* branch = dynamic_cast<TBranch*>(lob->At(b));
        std::string br_name = branch->GetName();
        bool include_branch = false;
        if (branch_names) {
            auto iter = std::find(branch_names->begin(), branch_names->end(), br_name);
            include_branch = (iter != branch_names->end());
        }
        else {
            const std::string wgt_br_prefix = "weight_";
            int compare_result = br_name.compare(0, wgt_br_prefix.size(), wgt_br_prefix);
            include_branch = (compare_result == 0);
        }
        if (!include_branch) continue;
        weight_map_[br_name] = MyPointer<std::vector<double>>();
        auto& wgt_vec = weight_map_.at(br_name);
        set_object_input_branch_address(in_tree, br_name, wgt_vec);
    }
}

void WeightHandler::set_branch_addresses(TTree& in_tree, const std::string& branch_name) {
    std::vector<std::string> br_names = {branch_name};
    set_branch_addresses(in_tree, &br_names);
}

void WeightHandler::add_branch(TTree& in_tree, const std::string& branch_name, bool throw_when_missing) {
    auto iter = weight_map_.find(branch_name);
    if (iter != weight_map_.end()) return;
    TBranch* br = in_tree.GetBranch(branch_name.c_str());
    if (!br) {
        if (throw_when_missing) throw std::runtime_error("Missing TTree branch " + branch_name);
        return;
    }
    weight_map_[branch_name] = MyPointer<std::vector<double>>();
    auto& wgt_vec = weight_map_.at(branch_name);
    set_object_input_branch_address(in_tree, branch_name, wgt_vec);
}