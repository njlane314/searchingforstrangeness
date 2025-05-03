#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <limits>
#include <algorithm>
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "EventCategory.hh"
#include "FilePropertiesManager.hh"
#include "TMath.h"

using NFT = NtupleFileType;

void plot_topology_variables(const std::set<int>& runs, const std::string& output_file_name) {
    std::vector<std::string> var_list = {"nhit_u", "nhit_v", "nhit_w", "charge_u", "charge_v", "charge_w", "wirerange_u", "wirerange_v", "wirerange_w", "timerange_u", "timerange_v", "timerange_w"};
    std::set<std::string> integer_vars = {"nhit_u", "nhit_v", "nhit_w"};
    const EventCategoryInterpreter& eci = EventCategoryInterpreter::Instance();
    const FilePropertiesManager& fpm = FilePropertiesManager::Instance();
    constexpr std::array<NFT, 5> file_types = {NFT::kOnNuMI, NFT::kExtNuMI, NFT::kNumuMC, NFT::kStrangenessMC, NFT::kDirtMC};
    constexpr std::array<NFT, 3> mc_file_types = {NFT::kNumuMC, NFT::kStrangenessMC, NFT::kDirtMC};
    std::map<NFT, std::unique_ptr<TChain>> tchain_map;
    std::map<NFT, double> pot_map;
    std::map<NFT, long> trigger_map;
    for (const auto& type : file_types) {
        tchain_map.emplace(std::make_pair(type, std::make_unique<TChain>("analysis_tree")));
        pot_map[type] = 0.;
        trigger_map[type] = 0;
    }
    const auto& ntuple_map = fpm.ntuple_file_map();
    const auto& data_norm_map = fpm.data_norm_map();
    for (const int& run : runs) {
        const auto& run_map = ntuple_map.at(run);
        for (const auto& type : file_types) {
            const auto& ntuple_files = run_map.at(type);
            auto* tchain = tchain_map.at(type).get();
            double& total_pot = pot_map.at(type);
            long& total_triggers = trigger_map.at(type);
            for (const auto& file_name : ntuple_files) {
                tchain->Add(file_name.c_str());
                if (type == NFT::kOnNuMI || type == NFT::kExtNuMI) {
                    const auto& norm_info = data_norm_map.at(file_name);
                    total_triggers += norm_info.trigger_count_;
                    total_pot += norm_info.pot_;
                } else {
                    TFile temp_file(file_name.c_str(), "read");
                    TParameter<float>* temp_pot = nullptr;
                    temp_file.GetObject("summed_pot", temp_pot);
                    if (temp_pot) total_pot += temp_pot->GetVal();
                }
            }
        }
    }
    double trigs_on = trigger_map.at(NFT::kOnNuMI);
    double trigs_off = trigger_map.at(NFT::kExtNuMI);
    double ext_scale_factor = (trigs_off > 0) ? trigs_on / trigs_off : 1.0;
    std::string ext_scale_factor_str = std::to_string(ext_scale_factor);
    double on_pot = pot_map.at(NFT::kOnNuMI);
    std::map<std::string, double> overall_min;
    std::map<std::string, double> overall_max;
    for (const auto& var : var_list) {
        overall_min[var] = std::numeric_limits<double>::max();
        overall_max[var] = std::numeric_limits<double>::min();
        for (const auto& type : file_types) {
            TChain* chain = tchain_map.at(type).get();
            chain->Draw(var.c_str(), "sel_reco_vertex_in_fv", "goff");
            int n = chain->GetSelectedRows();
            if (n > 0) {
                double* values = chain->GetV1();
                double min_val = *std::min_element(values, values + n);
                double max_val = *std::max_element(values, values + n);
                if (min_val < overall_min[var]) overall_min[var] = min_val;
                if (max_val > overall_max[var]) overall_max[var] = max_val;
            }
        }
    }
    std::map<std::string, TH1D*> hist_on;
    std::map<std::string, TH1D*> hist_off;
    std::map<std::string, std::map<EventCategory, TH1D*>> hist_mc;
    for (const auto& var : var_list) {
        double min_val = overall_min[var];
        double max_val = overall_max[var];
        int nbins;
        double low_edge, high_edge;
        if (integer_vars.count(var)) {
            int imin = std::floor(min_val);
            int imax = std::ceil(max_val);
            nbins = imax - imin + 1;
            low_edge = imin - 0.5;
            high_edge = imax + 0.5;
        } else {
            nbins = 100;
            low_edge = min_val;
            high_edge = max_val;
        }
        hist_on[var] = new TH1D(("hist_on_" + var).c_str(), (var + ";Events").c_str(), nbins, low_edge, high_edge);
        hist_off[var] = new TH1D(("hist_off_" + var).c_str(), (var + ";Events").c_str(), nbins, low_edge, high_edge);
        hist_mc[var] = std::map<EventCategory, TH1D*>();
        for (const auto& ec_pair : eci.label_map()) {
            EventCategory ec = ec_pair.first;
            std::string ec_str = std::to_string(static_cast<int>(ec));
            hist_mc[var][ec] = new TH1D(("hist_mc_" + var + "_" + ec_str).c_str(), (var + ";Events").c_str(), nbins, low_edge, high_edge);
        }
    }
    std::string selection = "sel_reco_vertex_in_fv";
    for (const auto& var : var_list) {
        TChain* on_chain = tchain_map.at(NFT::kOnNuMI).get();
        on_chain->Draw((var + " >> hist_on_" + var).c_str(), selection.c_str(), "goff");
        TChain* off_chain = tchain_map.at(NFT::kExtNuMI).get();
        off_chain->Draw((var + " >> hist_off_" + var).c_str(), (selection + " * " + ext_scale_factor_str).c_str(), "goff");
        for (const auto& type : mc_file_types) {
            TChain* mc_ch = tchain_map.at(type).get();
            double mc_pot = pot_map.at(type);
            double mc_scale_factor = (mc_pot > 0) ? on_pot / mc_pot : 1.0;
            std::string mc_scale_factor_str = std::to_string(mc_scale_factor);
            for (const auto& ec_pair : eci.label_map()) {
                EventCategory ec = ec_pair.first;
                if (type != NFT::kStrangenessMC && ec == EventCategory::k_nu_mu_cc_with_strange) continue;
                std::string ec_str = std::to_string(static_cast<int>(ec));
                std::string weight = mc_scale_factor_str + " * (" + selection + " && event_category == " + ec_str + ")";
                mc_ch->Draw((var + " >> hist_mc_" + var + "_" + ec_str).c_str(), weight.c_str(), "goff");
            }
        }
    }
    TFile* out_file = new TFile(output_file_name.c_str(), "recreate");
    for (const auto& var : var_list) {
        THStack* stack = new THStack(("stack_" + var).c_str(), (var + ";Events").c_str());
        TH1D* h_off = hist_off[var];
        h_off->SetFillColor(kGray);
        stack->Add(h_off);
        for (const auto& ec_pair : eci.label_map()) {
            EventCategory ec = ec_pair.first;
            TH1D* h_mc = hist_mc[var][ec];
            int color = eci.color_code(ec);
            h_mc->SetFillColor(color);
            h_mc->SetLineColor(color);
            stack->Add(h_mc);
        }
        TCanvas* c = new TCanvas(("c_" + var).c_str(), var.c_str(), 800, 600);
        stack->Draw("hist");
        TH1D* h_on = hist_on[var];
        h_on->SetMarkerStyle(20);
        h_on->Draw("same E");
        TLegend* lg = new TLegend(0.7, 0.7, 0.9, 0.9);
        lg->AddEntry(h_on, "Beam-on Data", "lep");
        lg->AddEntry(h_off, "Beam-off Data", "f");
        for (const auto& ec_pair : eci.label_map()) {
            EventCategory ec = ec_pair.first;
            lg->AddEntry(hist_mc[var][ec], eci.label(ec).c_str(), "f");
        }
        lg->Draw();
        c->Write();
    }
    out_file->Close();
    delete out_file;
    for (auto& var_pair : hist_on) delete var_pair.second;
    for (auto& var_pair : hist_off) delete var_pair.second;
    for (auto& var_pair : hist_mc) {
        for (auto& ec_pair : var_pair.second) delete ec_pair.second;
    }
}

int main() {
    std::set<int> runs = {1};
    std::string output_file = "topology_plots.root";
    plot_topology_variables(runs, output_file);
    std::cout << "Plots saved to " << output_file << std::endl;
    return 0;
}