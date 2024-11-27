#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>

void plot_radial_density(const char* filename = "../prod_strange_resample_fhc_run2_fhc_reco2_reco2_trainingregionanalyser_10_new_analysis.root") {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("regana/TrainingRegionTree");
    if (!tree) {
        std::cerr << "Failed to retrieve tree from file." << std::endl;
        file->Close();
        return;
    }

    std::vector<float>* radii = nullptr;
    std::vector<float>* densities = nullptr;

    tree->SetBranchAddress("radii", &radii);
    tree->SetBranchAddress("radial_densities", &densities);

    TMultiGraph* mg = new TMultiGraph();

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (!radii || !densities || radii->empty() || densities->empty()) {
            continue;
        }

        float cumulative_density = 0.0;
        TGraph* g_cumulative = new TGraph();
        for (size_t j = 0; j < radii->size(); ++j) {
            cumulative_density += (*densities)[j];
            g_cumulative->SetPoint(g_cumulative->GetN(), (*radii)[j], cumulative_density);
        }

        g_cumulative->SetLineColor(i + 1);
        g_cumulative->SetLineWidth(2);
        mg->Add(g_cumulative);
    }

    TCanvas* c_multigraph = new TCanvas("c_multigraph", "Cumulative Charge Density", 800, 600);
    mg->SetTitle("Cumulative Charge Density;Radius (cm);Cumulative Charge Density (e^{-})");
    mg->Draw("AL");
    c_multigraph->SaveAs("cumulative_charge_density_multigraph.png");

    delete mg;
    file->Close();
}
