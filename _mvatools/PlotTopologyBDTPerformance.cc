#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

int main() {
    TFile* file = TFile::Open("tmva_output.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open output file." << std::endl;
        return 1;
    }
    TGraph* roc = (TGraph*)file->Get("dataset/MVA_XGBoost_rejBvsS");
    if (!roc) {
        std::cerr << "Error: Cannot find ROC curve." << std::endl;
        return 1;
    }
    TCanvas* c1 = new TCanvas("c1", "ROC Curve", 800, 600);
    roc->SetTitle("ROC Curve");
    roc->GetXaxis()->SetTitle("Signal Efficiency");
    roc->GetYaxis()->SetTitle("Background Rejection");
    roc->Draw("AL");
    c1->SaveAs("roc_curve.png");
    TH1* h_sig = (TH1*)file->Get("dataset/MVA_XGBoost_S");
    TH1* h_bkg = (TH1*)file->Get("dataset/MVA_XGBoost_B");
    if (!h_sig || !h_bkg) {
        std::cerr << "Error: Cannot find MVA output histograms." << std::endl;
        return 1;
    }
    TCanvas* c2 = new TCanvas("c2", "MVA Output", 800, 600);
    h_sig->Scale(1.0 / h_sig->Integral());
    h_bkg->Scale(1.0 / h_bkg->Integral());
    h_sig->SetLineColor(kBlue);
    h_bkg->SetLineColor(kRed);
    h_sig->SetTitle("BDT Output");
    h_sig->GetXaxis()->SetTitle("BDT Score");
    h_sig->GetYaxis()->SetTitle("Normalized Events");
    h_sig->Draw("HIST");
    h_bkg->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_sig, "Signal", "l");
    leg->AddEntry(h_bkg, "Background", "l");
    leg->Draw();
    c2->SaveAs("mva_output.png");
    TH1* var_imp = (TH1*)file->Get("dataset/Method_XGBoost/XGBoost/VariableImportance");
    if (!var_imp) {
        std::cerr << "Error: Cannot find variable importance histogram." << std::endl;
        return 1;
    }
    TCanvas* c3 = new TCanvas("c3", "Variable Importance", 1000, 600);
    var_imp->SetTitle("Variable Importance");
    var_imp->GetXaxis()->SetTitle("Variable");
    var_imp->GetYaxis()->SetTitle("Importance");
    var_imp->Draw("HIST");
    c3->SaveAs("variable_importance.png");
    delete c1;
    delete c2;
    delete c3;
    file->Close();
    delete file;
    return 0;
}