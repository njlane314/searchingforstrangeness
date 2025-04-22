#include "Config.h"
#include "Sample.h"
#include "Event.h"
#include "Plotter.h"

#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <utility>
#include <functional>
#include <array>

int main() {
    std::string configName = "default";
    auto config = analysis::Config::GetConfig(configName);
    std::vector<std::unique_ptr<analysis::Sample>> samples;

    if (config.empty())
        return 1;

    for (const auto& [type, files] : config) {
        auto sample = std::make_unique<analysis::Sample>(type);
        for (const auto& file : files) {
            try {
                sample->AddFile(file);
            } catch (const std::runtime_error& e) {
                std::cerr << e.what() << std::endl;
                continue;
            }
        }
        samples.push_back(std::move(sample));
    }

    std::vector<analysis::InteractionGroup> allGroups = {
        analysis::InteractionGroup::QE,
        analysis::InteractionGroup::Res,
        analysis::InteractionGroup::DIS,
        analysis::InteractionGroup::Coh,
        analysis::InteractionGroup::MEC,
        analysis::InteractionGroup::Other
    };

    std::array<std::string, 3> view_names = {"U", "V", "W"};
    std::map<analysis::InteractionGroup, TH1F*> completeness_hit_hists_U;
    // std::array<std::map<analysis::InteractionGroup, TH1F*>, 3> completeness_hit_hists;
    // std::array<std::map<analysis::InteractionGroup, TH1F*>, 3> completeness_total_hits_hists;
    // std::array<std::map<analysis::InteractionGroup, TH1F*>, 3> exclusivity_ratio_hists;
    // std::array<std::map<analysis::InteractionGroup, TH1F*>, 3> hit_exclusivity_fraction_hists;
    // std::array<std::map<analysis::InteractionGroup, TH1F*>, 3> particle_start_active_hists;
    // std::array<std::map<analysis::InteractionGroup, TH1F*>, 3> particle_end_active_hists;
    // std::map<analysis::InteractionGroup, TH1F*> vertex_distance_hists;

    for (const auto& group : allGroups) {
        std::string ch_name = "completeness_hit_U_" + analysis::GetInteractionGroupName(group);
        completeness_hit_hists_U[group] = new TH1F(ch_name.c_str(), "", 20, 0, 1);
    }
    /*
    for (int view = 0; view < 3; ++view) {
        for (const auto& group : allGroups) {
            std::string ch_name = "completeness_hit_" + view_names[view] + "_" + analysis::GetInteractionGroupName(group);
            completeness_hit_hists[view][group] = new TH1F(ch_name.c_str(), "", 20, 0, 1);
            std::string cth_name = "completeness_total_hits_" + view_names[view] + "_" + analysis::GetInteractionGroupName(group);
            completeness_total_hits_hists[view][group] = new TH1F(cth_name.c_str(), "", 50, 0, 5000);
            std::string er_name = "exclusivity_ratio_" + view_names[view] + "_" + analysis::GetInteractionGroupName(group);
            exclusivity_ratio_hists[view][group] = new TH1F(er_name.c_str(), "", 20, 0, 1);
            std::string hef_name = "hit_exclusivity_fraction_" + view_names[view] + "_" + analysis::GetInteractionGroupName(group);
            hit_exclusivity_fraction_hists[view][group] = new TH1F(hef_name.c_str(), "", 20, 0, 1);
            std::string psa_name = "particle_start_active_" + view_names[view] + "_" + analysis::GetInteractionGroupName(group);
            particle_start_active_hists[view][group] = new TH1F(psa_name.c_str(), "", 2, 0, 2);
            std::string pea_name = "particle_end_active_" + view_names[view] + "_" + analysis::GetInteractionGroupName(group);
            particle_end_active_hists[view][group] = new TH1F(pea_name.c_str(), "", 2, 0, 2);
        }
    }

    for (const auto& group : allGroups) {
        std::string vd_name = "vertex_distance_" + analysis::GetInteractionGroupName(group);
        vertex_distance_hists[group] = new TH1F(vd_name.c_str(), "", 50, 0, 100);
    }
    */

    double target_pot = 1e20;
    plot::Plotter<TH1F> plotter(samples, target_pot);

    auto filler = [&completeness_hit_hists_U](const analysis::Event& event, double scale_factor) {
        if (event.is_signal) {
            analysis::InteractionGroup group = analysis::GetInteractionGroup(event.interaction);
            size_t n_signatures = event.signature_types->size();
            for (size_t i = 0; i < n_signatures; ++i) {
                if (i < (*event.completeness_hit[0]).size()) {
                    completeness_hit_hists_U[group]->Fill((*event.completeness_hit[0])[i], scale_factor);
                }
            }
        }
    };
    /*
    auto filler = [&completeness_hit_hists, &completeness_total_hits_hists,
                   &exclusivity_ratio_hists, &hit_exclusivity_fraction_hists, &particle_start_active_hists,
                   &particle_end_active_hists, &vertex_distance_hists](const analysis::Event& event, double scale_factor) {
        if (event.is_signal) {
            analysis::InteractionGroup group = analysis::GetInteractionGroup(event.interaction);
            size_t n_signatures = event.signature_types->size();
            for (size_t i = 0; i < n_signatures; ++i) {
                for (int view = 0; view < 3; ++view) {
                    if (i < (*event.completeness_hit[view]).size())
                        completeness_hit_hists[view][group]->Fill((*event.completeness_hit[view])[i], scale_factor);
                    if (i < (*event.completeness_total_hits[view]).size())
                        completeness_total_hits_hists[view][group]->Fill((*event.completeness_total_hits[view])[i], scale_factor);
                    if (i < (*event.exclusivity_ratio[view]).size())
                        exclusivity_ratio_hists[view][group]->Fill((*event.exclusivity_ratio[view])[i], scale_factor);
                    if (i < (*event.hit_exclusivity_fraction[view]).size())
                        hit_exclusivity_fraction_hists[view][group]->Fill((*event.hit_exclusivity_fraction[view])[i], scale_factor);
                }
                if (i < (*event.vertex_distance).size())
                    vertex_distance_hists[group]->Fill((*event.vertex_distance)[i], scale_factor);
            }
            for (size_t i = 0; i < (*event.particle_start_active[0]).size(); ++i) {
                for (int view = 0; view < 3; ++view) {
                    if (i < (*event.particle_start_active[view]).size())
                        particle_start_active_hists[view][group]->Fill((*event.particle_start_active[view])[i], scale_factor);
                    if (i < (*event.particle_end_active[view]).size())
                        particle_end_active_hists[view][group]->Fill((*event.particle_end_active[view])[i], scale_factor);
                }
            }
        }
    };
    */

    plotter.FillHistogram(filler);

    std::map<analysis::InteractionGroup, int> groupColors = {
        {analysis::InteractionGroup::QE, kBlue},
        {analysis::InteractionGroup::Res, kRed},
        {analysis::InteractionGroup::DIS, kGreen},
        {analysis::InteractionGroup::Coh, kMagenta},
        {analysis::InteractionGroup::MEC, kCyan},
        {analysis::InteractionGroup::Other, kGray}
    };

    TCanvas* ch_canvas = new TCanvas("ch_canvas_U", "", 800, 600);
    THStack* ch_stack = new THStack("completeness_hit_U", "Completeness Hit - View U");
    TLegend* ch_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
    ch_legend->SetBorderSize(0);
    ch_legend->SetNColumns(3);

    for (const auto& group : allGroups) {
        TH1F* hist = completeness_hit_hists_U[group];
        if (hist->GetEntries() > 0) {
            hist->SetFillColor(groupColors[group]);
            hist->SetLineColor(kBlack);
            ch_stack->Add(hist);
            ch_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
        }
    }
    ch_stack->Draw("HIST");
    ch_stack->GetXaxis()->SetTitle("Completeness Hit");
    ch_stack->GetYaxis()->SetTitle("Entries");
    ch_legend->Draw();
    ch_canvas->SaveAs("plots/completeness_hit_U.png");

    /*
    for (int view = 0; view < 3; ++view) {
        TCanvas* ch_canvas = new TCanvas(("ch_canvas_" + view_names[view]).c_str(), "", 800, 600);
        THStack* ch_stack = new THStack(("completeness_hit_" + view_names[view]).c_str(), ("Completeness Hit - View " + view_names[view]).c_str());
        TLegend* ch_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
        ch_legend->SetBorderSize(0);
        ch_legend->SetNColumns(3);
        for (const auto& group : allGroups) {
            TH1F* hist = completeness_hit_hists[view][group];
            if (hist->GetEntries() > 0) {
                hist->SetFillColor(groupColors[group]);
                hist->SetLineColor(kBlack);
                ch_stack->Add(hist);
                ch_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
            }
        }
        ch_stack->Draw("HIST");
        ch_stack->GetXaxis()->SetTitle("Completeness Hit");
        ch_stack->GetYaxis()->SetTitle("Entries");
        ch_legend->Draw();
        ch_canvas->SaveAs(("plots/completeness_hit_" + view_names[view] + ".png").c_str());

        TCanvas* cth_canvas = new TCanvas(("cth_canvas_" + view_names[view]).c_str(), "", 800, 600);
        THStack* cth_stack = new THStack(("completeness_total_hits_" + view_names[view]).c_str(), ("Completeness Total Hits - View " + view_names[view]).c_str());
        TLegend* cth_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
        cth_legend->SetBorderSize(0);
        cth_legend->SetNColumns(3);
        for (const auto& group : allGroups) {
            TH1F* hist = completeness_total_hits_hists[view][group];
            if (hist->GetEntries() > 0) {
                hist->SetFillColor(groupColors[group]);
                hist->SetLineColor(kBlack);
                cth_stack->Add(hist);
                cth_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
            }
        }
        cth_stack->Draw("HIST");
        cth_stack->GetXaxis()->SetTitle("Completeness Total Hits");
        cth_stack->GetYaxis()->SetTitle("Entries");
        cth_legend->Draw();
        cth_canvas->SaveAs(("plots/completeness_total_hits_" + view_names[view] + ".png").c_str());

        TCanvas* er_canvas = new TCanvas(("er_canvas_" + view_names[view]).c_str(), "", 800, 600);
        THStack* er_stack = new THStack(("exclusivity_ratio_" + view_names[view]).c_str(), ("Exclusivity Ratio - View " + view_names[view]).c_str());
        TLegend* er_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
        er_legend->SetBorderSize(0);
        er_legend->SetNColumns(3);
        for (const auto& group : allGroups) {
            TH1F* hist = exclusivity_ratio_hists[view][group];
            if (hist->GetEntries() > 0) {
                hist->SetFillColor(groupColors[group]);
                hist->SetLineColor(kBlack);
                er_stack->Add(hist);
                er_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
            }
        }
        er_stack->Draw("HIST");
        er_stack->GetXaxis()->SetTitle("Exclusivity Ratio");
        er_stack->GetYaxis()->SetTitle("Entries");
        er_legend->Draw();
        er_canvas->SaveAs(("plots/exclusivity_ratio_" + view_names[view] + ".png").c_str());

        TCanvas* hef_canvas = new TCanvas(("hef_canvas_" + view_names[view]).c_str(), "", 800, 600);
        THStack* hef_stack = new THStack(("hit_exclusivity_fraction_" + view_names[view]).c_str(), ("Hit Exclusivity Fraction - View " + view_names[view]).c_str());
        TLegend* hef_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
        hef_legend->SetBorderSize(0);
        hef_legend->SetNColumns(3);
        for (const auto& group : allGroups) {
            TH1F* hist = hit_exclusivity_fraction_hists[view][group];
            if (hist->GetEntries() > 0) {
                hist->SetFillColor(groupColors[group]);
                hist->SetLineColor(kBlack);
                hef_stack->Add(hist);
                hef_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
            }
        }
        hef_stack->Draw("HIST");
        hef_stack->GetXaxis()->SetTitle("Hit Exclusivity Fraction");
        hef_stack->GetYaxis()->SetTitle("Entries");
        hef_legend->Draw();
        hef_canvas->SaveAs(("plots/hit_exclusivity_fraction_" + view_names[view] + ".png").c_str());

        TCanvas* psa_canvas = new TCanvas(("psa_canvas_" + view_names[view]).c_str(), "", 800, 600);
        THStack* psa_stack = new THStack(("particle_start_active_" + view_names[view]).c_str(), ("Particle Start Active - View " + view_names[view]).c_str());
        TLegend* psa_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
        psa_legend->SetBorderSize(0);
        psa_legend->SetNColumns(3);
        for (const auto& group : allGroups) {
            TH1F* hist = particle_start_active_hists[view][group];
            if (hist->GetEntries() > 0) {
                hist->SetFillColor(groupColors[group]);
                hist->SetLineColor(kBlack);
                psa_stack->Add(hist);
                psa_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
            }
        }
        psa_stack->Draw("HIST");
        psa_stack->GetXaxis()->SetTitle("Particle Start Active (0=Inactive, 1=Active)");
        psa_stack->GetYaxis()->SetTitle("Entries");
        psa_legend->Draw();
        psa_canvas->SaveAs(("plots/particle_start_active_" + view_names[view] + ".png").c_str());

        TCanvas* pea_canvas = new TCanvas(("pea_canvas_" + view_names[view]).c_str(), "", 800, 600);
        THStack* pea_stack = new THStack(("particle_end_active_" + view_names[view]).c_str(), ("Particle End Active - View " + view_names[view]).c_str());
        TLegend* pea_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
        pea_legend->SetBorderSize(0);
        pea_legend->SetNColumns(3);
        for (const auto& group : allGroups) {
            TH1F* hist = particle_end_active_hists[view][group];
            if (hist->GetEntries() > 0) {
                hist->SetFillColor(groupColors[group]);
                hist->SetLineColor(kBlack);
                pea_stack->Add(hist);
                pea_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
            }
        }
        pea_stack->Draw("HIST");
        pea_stack->GetXaxis()->SetTitle("Particle End Active (0=Inactive, 1=Active)");
        pea_stack->GetYaxis()->SetTitle("Entries");
        pea_legend->Draw();
        pea_canvas->SaveAs(("plots/particle_end_active_" + view_names[view] + ".png").c_str());
    }

    TCanvas* vd_canvas = new TCanvas("vd_canvas", "", 800, 600);
    THStack* vd_stack = new THStack("vertex_distance", "Vertex Distance");
    TLegend* vd_legend = new TLegend(0.1, 0.91, 0.9, 0.99);
    vd_legend->SetBorderSize(0);
    vd_legend->SetNColumns(3);
    for (const auto& group : allGroups) {
        TH1F* hist = vertex_distance_hists[group];
        if (hist->GetEntries() > 0) {
            hist->SetFillColor(groupColors[group]);
            hist->SetLineColor(kBlack);
            vd_stack->Add(hist);
            vd_legend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
        }
    }
    vd_stack->Draw("HIST");
    vd_stack->GetXaxis()->SetTitle("Vertex Distance [cm]");
    vd_stack->GetYaxis()->SetTitle("Entries");
    vd_legend->Draw();
    vd_canvas->SaveAs("plots/vertex_distance.png");
    */

    for (auto& [group, hist] : completeness_hit_hists_U) {
        delete hist;
    }
    /*
    for (int view = 0; view < 3; ++view) {
        for (auto& [group, hist] : completeness_hit_hists[view]) delete hist;
        for (auto& [group, hist] : completeness_total_hits_hists[view]) delete hist;
        for (auto& [group, hist] : exclusivity_ratio_hists[view]) delete hist;
        for (auto& [group, hist] : hit_exclusivity_fraction_hists[view]) delete hist;
        for (auto& [group, hist] : particle_start_active_hists[view]) delete hist;
        for (auto& [group, hist] : particle_end_active_hists[view]) delete hist;
    }
    for (auto& [group, hist] : vertex_distance_hists) delete hist;
    */

    return 0;
}