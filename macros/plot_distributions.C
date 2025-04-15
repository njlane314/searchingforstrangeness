#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <iostream>

void plot_distributions(const char* filename) {
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("selectionfilter/SelectionFilter"));
    if (!tree) {
        std::cerr << "Error: Tree 'selectionfilter/SelectionFilter' not found" << std::endl;
        file->Close();
        delete file;
        return;
    }

    TCanvas* c1 = new TCanvas("c1", "Event Type", 800, 600);
    TH1F* h1 = new TH1F("h_event_type", "Event Type Distribution", 100, 1, -1);
    tree->Draw("event_type>>h_event_type");
    h1->SetTitle("Event Type Distribution");
    h1->GetXaxis()->SetTitle("Event Type");
    h1->GetYaxis()->SetTitle("Events");
    c1->SaveAs("event_type.png");
    delete c1;

    TCanvas* c2 = new TCanvas("c2", "Signature Types", 800, 600);
    TH1F* h2 = new TH1F("h_signature_types", "Signature Types Distribution", 100, 1, -1);
    tree->Draw("signature_types>>h_signature_types");
    h2->SetTitle("Signature Types Distribution");
    h2->GetXaxis()->SetTitle("Signature Type");
    h2->GetYaxis()->SetTitle("Counts");
    c2->SaveAs("signature_types.png");
    delete c2;

    TCanvas* c3 = new TCanvas("c3", "Clarity U", 800, 600);
    TH1F* h3 = new TH1F("h_clarity_U", "Clarity U Distribution", 100, 1, -1);
    tree->Draw("clarity_U>>h_clarity_U");
    h3->SetTitle("Clarity U Distribution");
    h3->GetXaxis()->SetTitle("Clarity U");
    h3->GetYaxis()->SetTitle("Counts");
    c3->SaveAs("clarity_U.png");
    delete c3;

    TCanvas* c4 = new TCanvas("c4", "Completeness Hit U", 800, 600);
    TH1F* h4 = new TH1F("h_completeness_hit_U", "Completeness Hit U Distribution", 100, 1, -1);
    tree->Draw("completeness_hit_U>>h_completeness_hit_U");
    h4->SetTitle("Completeness Hit U Distribution");
    h4->GetXaxis()->SetTitle("Completeness Hit U");
    h4->GetYaxis()->SetTitle("Counts");
    c4->SaveAs("completeness_hit_U.png");
    delete c4;

    file->Close();
    delete file;
}