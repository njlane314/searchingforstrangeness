#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

void plot_region_ana(const char* filename = "prod_strange_resample_fhc_run2_fhc_reco2_reco2_trainingregionanalyser_100_new_analysis.root") {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "-- Error: Could not open file " << filename << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("regana/TrainingRegionTree");
    if (!tree) {
        std::cerr << "-- Error: Could not find TTree 'TrainingRegionTree' in file " << filename << std::endl;
        file->Close();
        return;
    }

    TH1F* h_hits_registered_u = new TH1F("h_hits_registered_u", ";Hits Registered;Entries/gin", 25, 0, -1);
    TH1F* h_hits_registered_v = new TH1F("h_hits_registered_v", ";Hits Registered;Entries/gin", 25, 0, -1);
    TH1F* h_hits_registered_w = new TH1F("h_hits_registered_w", ";Hits Registered;Entries/gin", 25, 0, -1);

    TH1F* h_charge_collected_u = new TH1F("h_charge_collected_u", ";Charged Collected;Entries/gin", 22, 0, -1);
    TH1F* h_charge_collected_v = new TH1F("h_charge_collected_v", ";Charged Collected;Entries/gin", 22, 0, -1);
    TH1F* h_charge_collected_w = new TH1F("h_charge_collected_w", ";Charged Collected;Entries/gin", 22, 0, -1);

    int hit_count_u, hit_count_v, hit_count_w;
    float total_charge_u, total_charge_v, total_charge_w;

    tree->SetBranchAddress("hit_count_u", &hit_count_u);
    tree->SetBranchAddress("hit_count_v", &hit_count_v);
    tree->SetBranchAddress("hit_count_w", &hit_count_w);
    tree->SetBranchAddress("total_charge_u", &total_charge_u);
    tree->SetBranchAddress("total_charge_v", &total_charge_v);
    tree->SetBranchAddress("total_charge_w", &total_charge_w);

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        h_hits_registered_u->Fill(hit_count_u);
        h_hits_registered_v->Fill(hit_count_v);
        h_hits_registered_w->Fill(hit_count_w);

        h_charge_collected_u->Fill(total_charge_u);
        h_charge_collected_v->Fill(total_charge_v);
        h_charge_collected_w->Fill(total_charge_w);
    }

    TCanvas* c_hits_registered = new TCanvas("hits", "", 800, 600);
    h_hits_registered_u->SetLineColor(kRed);
    h_hits_registered_v->SetLineColor(kBlue);
    h_hits_registered_w->SetLineColor(kGreen);
    h_hits_registered_u->SetLineWidth(2);
    h_hits_registered_v->SetLineWidth(2);
    h_hits_registered_w->SetLineWidth(2);
    h_hits_registered_u->SetStats(0);
    h_hits_registered_v->SetStats(0);
    h_hits_registered_w->SetStats(0);
    double max_val = std::max({h_hits_registered_u->GetMaximum(), h_hits_registered_v->GetMaximum(), h_hits_registered_w->GetMaximum()});
    h_hits_registered_u->SetMaximum(1.1 * max_val);
    h_hits_registered_u->Draw("HIST E");
    h_hits_registered_v->Draw("HIST E SAME");
    h_hits_registered_w->Draw("HIST E SAME");

    TLegend* legend = new TLegend(0.1, 0.91, 0.9, 0.99);
    legend->SetNColumns(3);
    legend->SetBorderSize(0);
    legend->AddEntry(h_hits_registered_u, "U", "l");
    legend->AddEntry(h_hits_registered_v, "V", "l");
    legend->AddEntry(h_hits_registered_w, "W", "l");
    legend->Draw();

    TCanvas* c_charged_collected = new TCanvas("charge", "", 800, 600);
    h_charge_collected_u->SetLineColor(kRed);
    h_charge_collected_v->SetLineColor(kBlue);
    h_charge_collected_w->SetLineColor(kGreen);
    h_charge_collected_u->SetLineWidth(2);
    h_charge_collected_v->SetLineWidth(2);
    h_charge_collected_w->SetLineWidth(2);
    h_charge_collected_u->SetStats(0);
    h_charge_collected_v->SetStats(0);
    h_charge_collected_w->SetStats(0);
    double max_val2 = std::max({h_charge_collected_u->GetMaximum(), h_charge_collected_v->GetMaximum(), h_charge_collected_w->GetMaximum()});
    h_charge_collected_u->SetMaximum(1.1 * max_val2);
    h_charge_collected_u->Draw("HIST E");
    h_charge_collected_v->Draw("HIST E SAME");
    h_charge_collected_w->Draw("HIST E SAME");

    TLegend* legend2 = new TLegend(0.1, 0.91, 0.9, 0.99);
    legend2->SetNColumns(3);
    legend2->SetBorderSize(0);
    legend2->AddEntry(h_charge_collected_u, "U", "l");
    legend2->AddEntry(h_charge_collected_v, "V", "l");
    legend2->AddEntry(h_charge_collected_w, "W", "l");
    legend2->Draw();

    c_hits_registered->SaveAs("hits_registered.png");
    c_charged_collected->SaveAs("charge_collected.png");

    delete h_hits_registered_u;
    delete h_hits_registered_v;
    delete h_hits_registered_w;
    delete h_charge_collected_u;
    delete h_charge_collected_v;
    delete h_charge_collected_w;

    file->Close();
}
