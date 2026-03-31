#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

bool require_branch(TTree* tree, const char* branch_name) {
    if (tree->GetBranch(branch_name)) return true;
    std::cerr << "Missing required branch: " << branch_name << std::endl;
    return false;
}

void warn_size(Long64_t entry,
               const char* lhs_name,
               std::size_t lhs_size,
               const char* rhs_name,
               std::size_t rhs_size) {
    std::cerr << "entry " << entry << ": size mismatch for "
              << lhs_name << " (" << lhs_size << ") and "
              << rhs_name << " (" << rhs_size << ")" << std::endl;
}

void warn_missing_ref(Long64_t entry,
                      const char* label,
                      int value) {
    std::cerr << "entry " << entry << ": " << label
              << " references missing id " << value << std::endl;
}

void check_plane_refs(Long64_t entry,
                      const char* label,
                      const std::vector<int>& plane_uid,
                      const std::vector<int>& refs) {
    std::unordered_set<int> plane_ids(plane_uid.begin(), plane_uid.end());
    for (const int ref : refs) {
        if (!plane_ids.count(ref)) warn_missing_ref(entry, label, ref);
    }
}

void check_comp_refs(Long64_t entry,
                     const char* label,
                     const std::vector<int>& comp_ids,
                     const std::vector<int>& refs) {
    std::unordered_set<int> ids(comp_ids.begin(), comp_ids.end());
    for (const int ref : refs) {
        if (ref >= 0 && !ids.count(ref)) warn_missing_ref(entry, label, ref);
    }
}

} // namespace

void print_graph_summary(const char* filename = "graphanalysis_dev.root",
                         Long64_t max_entries_to_print = 10) {
    TFile file(filename, "READ");
    if (file.IsZombie()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    auto* tree = dynamic_cast<TTree*>(file.Get("EventSelectionFilter"));
    if (!tree) {
        std::cerr << "Could not find TTree 'EventSelectionFilter' in " << filename << std::endl;
        return;
    }

    const char* required[] = {
        "run",
        "sub",
        "evt",
        "selected",
        "da_plane_uid",
        "da_plane_has_primary_vtx",
        "da_dead_plane_uid",
        "da_hit_plane_uid",
        "da_hit_component_id",
        "da_hit_activity_id",
        "da_comp_plane_uid",
        "da_comp_id",
        "da_comp_label",
        "da_comp_activity_id",
        "da_comp_parent_id",
        "da_edge_plane_uid",
        "da_edge_src_component",
        "da_edge_dst_component"
    };
    for (const char* branch_name : required) {
        if (!require_branch(tree, branch_name)) return;
    }

    int run = -1;
    int sub = -1;
    int evt = -1;
    int selected = 0;

    std::vector<int>* da_plane_uid = nullptr;
    std::vector<int>* da_plane_has_primary_vtx = nullptr;
    std::vector<int>* da_dead_plane_uid = nullptr;
    std::vector<int>* da_hit_plane_uid = nullptr;
    std::vector<int>* da_hit_component_id = nullptr;
    std::vector<int>* da_hit_activity_id = nullptr;
    std::vector<int>* da_comp_plane_uid = nullptr;
    std::vector<int>* da_comp_id = nullptr;
    std::vector<int>* da_comp_label = nullptr;
    std::vector<int>* da_comp_activity_id = nullptr;
    std::vector<int>* da_comp_parent_id = nullptr;
    std::vector<int>* da_edge_plane_uid = nullptr;
    std::vector<int>* da_edge_src_component = nullptr;
    std::vector<int>* da_edge_dst_component = nullptr;

    tree->SetBranchAddress("run", &run);
    tree->SetBranchAddress("sub", &sub);
    tree->SetBranchAddress("evt", &evt);
    tree->SetBranchAddress("selected", &selected);
    tree->SetBranchAddress("da_plane_uid", &da_plane_uid);
    tree->SetBranchAddress("da_plane_has_primary_vtx", &da_plane_has_primary_vtx);
    tree->SetBranchAddress("da_dead_plane_uid", &da_dead_plane_uid);
    tree->SetBranchAddress("da_hit_plane_uid", &da_hit_plane_uid);
    tree->SetBranchAddress("da_hit_component_id", &da_hit_component_id);
    tree->SetBranchAddress("da_hit_activity_id", &da_hit_activity_id);
    tree->SetBranchAddress("da_comp_plane_uid", &da_comp_plane_uid);
    tree->SetBranchAddress("da_comp_id", &da_comp_id);
    tree->SetBranchAddress("da_comp_label", &da_comp_label);
    tree->SetBranchAddress("da_comp_activity_id", &da_comp_activity_id);
    tree->SetBranchAddress("da_comp_parent_id", &da_comp_parent_id);
    tree->SetBranchAddress("da_edge_plane_uid", &da_edge_plane_uid);
    tree->SetBranchAddress("da_edge_src_component", &da_edge_src_component);
    tree->SetBranchAddress("da_edge_dst_component", &da_edge_dst_component);

    Long64_t total_planes = 0;
    Long64_t total_dead = 0;
    Long64_t total_hits = 0;
    Long64_t total_components = 0;
    Long64_t total_edges = 0;
    Long64_t total_detached = 0;
    Long64_t total_unresolved = 0;
    Long64_t total_prompt = 0;
    Long64_t total_plane_pv = 0;

    const auto nentries = tree->GetEntries();
    const auto nprint = std::min(nentries, max_entries_to_print);

    std::cout << "Reading " << nentries << " entries from " << filename << std::endl;

    for (Long64_t entry = 0; entry < nentries; ++entry) {
        tree->GetEntry(entry);

        if (da_plane_uid->size() != da_plane_has_primary_vtx->size()) {
            warn_size(entry,
                      "da_plane_uid",
                      da_plane_uid->size(),
                      "da_plane_has_primary_vtx",
                      da_plane_has_primary_vtx->size());
        }

        check_plane_refs(entry, "da_dead_plane_uid", *da_plane_uid, *da_dead_plane_uid);
        check_plane_refs(entry, "da_hit_plane_uid", *da_plane_uid, *da_hit_plane_uid);
        check_plane_refs(entry, "da_comp_plane_uid", *da_plane_uid, *da_comp_plane_uid);
        check_plane_refs(entry, "da_edge_plane_uid", *da_plane_uid, *da_edge_plane_uid);

        check_comp_refs(entry, "da_hit_component_id", *da_comp_id, *da_hit_component_id);
        check_comp_refs(entry, "da_hit_activity_id", *da_comp_id, *da_hit_activity_id);
        check_comp_refs(entry, "da_comp_activity_id", *da_comp_id, *da_comp_activity_id);
        check_comp_refs(entry, "da_comp_parent_id", *da_comp_id, *da_comp_parent_id);
        check_comp_refs(entry, "da_edge_src_component", *da_comp_id, *da_edge_src_component);
        check_comp_refs(entry, "da_edge_dst_component", *da_comp_id, *da_edge_dst_component);

        const auto entry_prompt =
            std::count(da_comp_label->begin(), da_comp_label->end(), 0);
        const auto entry_detached =
            std::count(da_comp_label->begin(), da_comp_label->end(), 1);
        const auto entry_unresolved =
            std::count(da_comp_label->begin(), da_comp_label->end(), 2);
        const auto entry_plane_pv =
            std::count(da_plane_has_primary_vtx->begin(), da_plane_has_primary_vtx->end(), 1);

        total_planes += da_plane_uid->size();
        total_dead += da_dead_plane_uid->size();
        total_hits += da_hit_plane_uid->size();
        total_components += da_comp_id->size();
        total_edges += da_edge_plane_uid->size();
        total_prompt += entry_prompt;
        total_detached += entry_detached;
        total_unresolved += entry_unresolved;
        total_plane_pv += entry_plane_pv;

        if (entry < nprint) {
            std::cout
                << "entry " << entry
                << " run/sub/evt=" << run << "/" << sub << "/" << evt
                << " selected=" << selected
                << " planes=" << da_plane_uid->size()
                << " dead=" << da_dead_plane_uid->size()
                << " hits=" << da_hit_plane_uid->size()
                << " comps=" << da_comp_id->size()
                << " edges=" << da_edge_plane_uid->size()
                << " pv_planes=" << entry_plane_pv
                << " prompt=" << entry_prompt
                << " detached=" << entry_detached
                << " unresolved=" << entry_unresolved
                << std::endl;
        }
    }

    std::cout << "Totals:"
              << " planes=" << total_planes
              << " dead=" << total_dead
              << " hits=" << total_hits
              << " comps=" << total_components
              << " edges=" << total_edges
              << " pv_planes=" << total_plane_pv
              << " prompt=" << total_prompt
              << " detached=" << total_detached
              << " unresolved=" << total_unresolved
              << std::endl;
}
