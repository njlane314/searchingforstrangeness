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

    float min_charge = std::numeric_limits<float>::max();
    float max_charge = std::numeric_limits<float>::min();
    float slice_total_charge; 

    for (const auto& sample : samples) {
        TChain* event_chain = sample->GetEventChain();
        if (!event_chain || event_chain->GetEntries() == 0) 
            continue;

        if (event_chain->SetBranchAddress("slice_total_charge", &slice_total_charge) < 0) 
            continue;

        Long64_t n_entries = event_chain->GetEntries();
        for (Long64_t i = 0; i < n_entries; ++i) {
            event_chain->GetEntry(i);
            if (slice_total_charge < 0) continue;
            if (slice_total_charge < min_charge) min_charge = slice_total_charge;
            if (slice_total_charge > max_charge) max_charge = slice_total_charge;
        }
    }

    std::unordered_map<analysis::InteractionType, TH1F*> histograms;
    for (const auto& type : analysis::allInteractionTypes) {
        std::string name = analysis::GetInteractionType(type);
        std::string histName = "hist_" + name;
        std::string histTitle = "Slice Total Charge - " + name;
        histograms[type] = new TH1F(histName.c_str(), histTitle.c_str(), 20, min_charge, max_charge);  
    }

    std::vector<analysis::InteractionGroup> allGroups = {
        analysis::InteractionGroup::QE,
        analysis::InteractionGroup::Res,
        analysis::InteractionGroup::DIS,
        analysis::InteractionGroup::Coh,
        analysis::InteractionGroup::MEC,
        analysis::InteractionGroup::Other
    };

    std::map<analysis::InteractionGroup, TH1F*> groupHistograms;
    for (const auto& group : allGroups) {
        std::string name = analysis::GetInteractionGroupName(group);
        std::string histName = "hist_group_" + name;
        std::string histTitle = "Slice Total Charge - " + name;
        groupHistograms[group] = new TH1F(histName.c_str(), histTitle.c_str(), 20, min_charge, max_charge);  
    }

    double target_pot = 1e21;
    plot::Plotter<TH1F> plotter(samples, target_pot);

    auto filler = [&histograms, &groupHistograms](const analysis::Event& event, double scale_factor) {
        if (event.slice_total_charge < 0) return;
        auto typeIt = histograms.find(event.interaction);
        if (typeIt != histograms.end()) {
            typeIt->second->Fill(event.slice_total_charge, scale_factor);
        }
        analysis::InteractionGroup group = analysis::GetInteractionGroup(event.interaction);
        auto groupIt = groupHistograms.find(group);
        if (groupIt != groupHistograms.end()) {
            groupIt->second->Fill(event.slice_total_charge, scale_factor);
        }
    };

    plotter.FillHistogram(filler);

    std::map<analysis::InteractionGroup, int> groupColors = {
        {analysis::InteractionGroup::QE, kBlue},
        {analysis::InteractionGroup::Res, kRed},
        {analysis::InteractionGroup::DIS, kGreen},
        {analysis::InteractionGroup::Coh, kMagenta},
        {analysis::InteractionGroup::MEC, kCyan},
        {analysis::InteractionGroup::Other, kGray}
    };

    THStack* stack = new THStack("slice_charge_stack", "");
    stack->SetTitle("");

    TLegend* legend = new TLegend(0.1, 0.91, 0.9, 0.99);
    legend->SetBorderSize(0);
    legend->SetNColumns(4);
    int color = kBlue;
    for (const auto& [type, hist] : histograms) {
        if (hist->GetEntries() == 0) continue;
        hist->SetFillColor(color);
        hist->SetLineColor(kBlack);
        stack->Add(hist);
        legend->AddEntry(hist, analysis::GetInteractionType(type).c_str(), "f");
        color++;
    }

    TObjArray* stackList = stack->GetStack();
    TH1* sumHist = nullptr;
    if (stackList && stackList->GetSize() > 0) {
        sumHist = (TH1*)stackList->Last()->Clone("sumHist");
        sumHist->SetFillStyle(0);
        sumHist->SetLineColor(kBlack);
        sumHist->SetMarkerStyle(0);
    }

    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    stack->Draw("HIST");
    if (sumHist) {
        sumHist->Draw("E SAME");
    }
    stack->GetXaxis()->SetTitle("Slice Total Charge");
    stack->GetYaxis()->SetTitle("Events/bin");
    legend->Draw();

    canvas->SaveAs("plots/slice_total_charge_by_interaction.png");

    THStack* groupStack = new THStack("slice_charge_group_stack", "");
    groupStack->SetTitle("");

    TLegend* groupLegend = new TLegend(0.1, 0.91, 0.9, 0.99);
    groupLegend->SetBorderSize(0);
    groupLegend->SetNColumns(5);
    for (const auto& [group, hist] : groupHistograms) {
        if (hist->GetEntries() == 0) continue;
        int groupColor = groupColors[group];
        hist->SetFillColor(groupColor);
        hist->SetLineColor(kBlack);
        groupStack->Add(hist);
        groupLegend->AddEntry(hist, analysis::GetInteractionGroupName(group).c_str(), "f");
    }

    TObjArray* groupStackList = groupStack->GetStack();
    TH1* groupSumHist = nullptr;
    if (groupStackList && groupStackList->GetSize() > 0) {
        groupSumHist = (TH1*)groupStackList->Last()->Clone("groupSumHist");
        groupSumHist->SetFillStyle(0);
        groupSumHist->SetLineColor(kBlack);
        groupSumHist->SetMarkerStyle(0);
    }

    TCanvas* groupCanvas = new TCanvas("groupCanvas", "", 800, 600);
    groupCanvas->SetLeftMargin(0.15);
    groupStack->Draw("HIST");
    if (groupSumHist) {
        groupSumHist->Draw("E SAME");
    }
    groupStack->GetXaxis()->SetTitle("Slice Total Charge");
    groupStack->GetYaxis()->SetTitle("Events/bin");
    groupLegend->Draw();

    groupCanvas->SaveAs("plots/slice_total_charge_by_interaction_group.png");

    for (auto& [type, hist] : histograms) {
        delete hist;
    }
    for (auto& [group, hist] : groupHistograms) {
        delete hist;
    }
    delete stack;
    delete legend;
    delete canvas;
    delete groupStack;
    delete groupLegend;
    delete groupCanvas;

    return 0;
}