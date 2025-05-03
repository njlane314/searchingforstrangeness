#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TParameter.h"
#include "TStyle.h"
#include "TPad.h"
#include "EventCategory.hh"
#include "FilePropertiesManager.hh"

using NFT = NtupleFileType;

TH2D* make_tally_histogram(const std::string& hist_name_prefix,
                           const std::string& sel,
                           const std::set<int>& runs,
                           const std::string& mc_event_weight = 1)
{
    std::string selection = '(' + sel + ')';
    int num_category_bins = static_cast<int>(EventCategory::kOther) + 3;
    std::string tally_hist_name = hist_name_prefix + "_tally_hist";
    TH2D* event_tally_hist = new TH2D(tally_hist_name.c_str(),
                                      "; event category; selected; events",
                                      num_category_bins, -2., num_category_bins - 2,
                                      2, 0., 2);

    const EventCategoryInterpreter& eci = EventCategoryInterpreter::Instance();
    const FilePropertiesManager& fpm = FilePropertiesManager::Instance();

    constexpr std::array<NFT, 5> file_types = {NFT::kOnNuMI, NFT::kExtNuMI,
                                               NFT::kNumuMC, NFT::kStrangenessMC, NFT::kDirtMC};

    constexpr std::array<NFT, 3> mc_file_types = {NFT::kNumuMC, NFT::kStrangenessMC, NFT::kDirtMC};

    std::map<NFT, std::unique_ptr<TChain>> tchain_map;
    std::map<NFT, double> pot_map;
    std::map<NFT, long> trigger_map;
    for (const auto& type : file_types) {
        tchain_map.emplace(std::make_pair(type, new TChain("stv_tree")));
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
                    double pot = temp_pot->GetVal();
                    total_pot += pot;
                }
            }
        }
    }

    double trigs_on = trigger_map.at(NFT::kOnNuMI);
    double trigs_off = trigger_map.at(NFT::kExtNuMI);
    double ext_scale_factor = trigs_on / trigs_off;
    std::string ext_scale_factor_str = std::to_string(ext_scale_factor);

    TChain* off_chain = tchain_map.at(NFT::kExtNuMI).get();
    off_chain->Draw((selection + " : -1 >>+" + tally_hist_name).c_str(),
                    ext_scale_factor_str.c_str(), "goff");

    TChain* on_chain = tchain_map.at(NFT::kOnNuMI).get();
    on_chain->Draw((selection + " : -2 >>+" + tally_hist_name).c_str(),
                   "1", "goff");

    for (const auto& type : mc_file_types) {
        TChain* mc_ch = tchain_map.at(type).get();
        double on_pot = pot_map.at(NFT::kOnNuMI);
        double mc_pot = pot_map.at(type);
        double mc_scale_factor = on_pot / mc_pot;
        std::string mc_scale_factor_str = std::to_string(mc_scale_factor);

        for (const auto& pair : eci.label_map()) {
            EventCategory ec = pair.first;
            std::string ec_str = std::to_string(static_cast<int>(ec));
            std::string weight_str = mc_event_weight + " * " + mc_scale_factor_str;

            if (type != NFT::kStrangenessMC && ec == EventCategory::k_nu_mu_cc_with_strange) {
                continue;  
            }

            std::string draw_cmd = selection + " && category == " + ec_str + " : " + ec_str + " >>+" + tally_hist_name;
            mc_ch->Draw(draw_cmd.c_str(), weight_str.c_str(), "goff");
        }
    }

    return event_tally_hist;
}

void calc_effpur(const TH2D& tally_hist, double& efficiency, double& purity, double& mc_purity)
{
    const int category_bin_offset = 3;
    const int signal_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_strange);
    const int first_signal_bin = signal_category + category_bin_offset;
    const int last_signal_bin = first_signal_bin;
    const int first_nonNuMI_category_bin = 2;
    const int first_mc_category_bin = category_bin_offset;
    int last_category_bin = tally_hist.GetNbinsX();
    const int missed_bin = 1;
    const int selected_bin = 2;

    double all_signal = tally_hist.Integral(first_signal_bin, last_signal_bin, missed_bin, selected_bin);
    double selected_signal = tally_hist.Integral(first_signal_bin, last_signal_bin, selected_bin, selected_bin);
    double all_selected = tally_hist.Integral(first_nonNuMI_category_bin, last_category_bin, selected_bin, selected_bin);
    double all_mc_selected = tally_hist.Integral(first_mc_category_bin, last_category_bin, selected_bin, selected_bin);

    efficiency = (all_signal > 0) ? selected_signal / all_signal : 0;
    purity = (all_selected > 0) ? selected_signal / all_selected : 0;
    mc_purity = (all_mc_selected > 0) ? selected_signal / all_mc_selected : 0;
}

void effpur_cutflow()
{
    const std::vector<std::string> selection_defs = {
        "sel_reco_vertex_in_fv"
    };

    const std::vector<std::string> selection_flags = {
        "sel\\char 95reco\\char 95vertex\\char 95in\\char 95FV"
    };

    size_t num_points = selection_defs.size();
    TGraph* eff_graph = new TGraph(num_points);
    TGraph* pur_graph = new TGraph(num_points);
    TGraph* mc_pur_graph = new TGraph(num_points);

    double eff, pur, mc_pur;
    for (size_t k = 0u; k < num_points; ++k) {
        const auto& selection = selection_defs.at(k);
        std::string k_str = std::to_string(k);
        auto* tally_hist = make_tally_histogram(k_str, selection, {1});
        calc_effpur(*tally_hist, eff, pur, mc_pur);

        eff_graph->SetPoint(k, k + 1, eff);
        pur_graph->SetPoint(k, k + 1, pur);
        mc_pur_graph->SetPoint(k, k + 1, mc_pur);

        std::cout << "selection = " << selection << '\n';
        std::cout << "eff = " << eff << '\n';
        std::cout << "pur = " << pur << '\n';
        std::cout << "mc_pur = " << mc_pur << '\n';
        std::cout << "\n\n";

        delete tally_hist;
    }

    const std::vector<std::string> bin_labels = {
        "in FV"
    };

    const std::vector<std::string> latex_bin_labels = {
        "in FV"
    };

    TCanvas* c1 = new TCanvas;
    c1->SetBottomMargin(0.21);
    eff_graph->SetTitle(";; efficiency or purity");
    eff_graph->SetLineColor(kBlue);
    eff_graph->SetMarkerColor(kBlue);
    eff_graph->SetLineWidth(3);
    eff_graph->SetMarkerStyle(20);
    eff_graph->GetYaxis()->SetRangeUser(0., 1.);
    eff_graph->Draw("alp");

    for (int b = 1; b <= bin_labels.size(); ++b) {
        eff_graph->GetHistogram()->GetXaxis()->SetBinLabel(
            eff_graph->GetHistogram()->FindBin(b), bin_labels.at(b - 1).c_str());
    }
    eff_graph->Draw("same");

    pur_graph->SetLineColor(kRed);
    pur_graph->SetMarkerColor(kRed);
    pur_graph->SetLineWidth(3);
    pur_graph->SetMarkerStyle(20);
    pur_graph->Draw("same lp");

    mc_pur_graph->SetLineColor(kBlack);
    mc_pur_graph->SetMarkerColor(kBlack);
    mc_pur_graph->SetLineWidth(3);
    mc_pur_graph->SetMarkerStyle(20);
    mc_pur_graph->Draw("same lp");

    TLegend* lg = new TLegend(0.65, 0.4, 0.85, 0.6);
    lg->AddEntry(eff_graph, "efficiency", "lp");
    lg->AddEntry(pur_graph, "purity", "lp");
    lg->AddEntry(mc_pur_graph, "MC purity", "lp");
    lg->Draw("same");

    std::map<std::string, std::vector<std::string>> pgfplots_table;
    dump_tgraph("eff", *eff_graph, pgfplots_table);
    dump_tgraph("pur", *pur_graph, pgfplots_table);
    dump_tgraph("mcpur", *mc_pur_graph, pgfplots_table);

    pgfplots_table["x_label"] = std::vector<std::string>();
    for (const auto& label : latex_bin_labels) {
        std::string entry = '{' + label + '}';
        pgfplots_table.at("x_label").push_back(entry);
    }

    pgfplots_table["sel_flag"] = std::vector<std::string>();
    for (const auto& flag : selection_flags) {
        std::string entry = '{' + flag + '}';
        pgfplots_table.at("sel_flag").push_back(entry);
    }

    write_pgfplots_file("cutflow.txt", pgfplots_table);
}