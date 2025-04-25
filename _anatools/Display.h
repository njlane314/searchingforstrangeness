#ifndef DISPLAY_H
#define DISPLAY_H

#include "Config.h"
#include "Sample.h"
#include "Event.h"
#include "Plotter.h"

#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <memory>
#include <utility>
#include <functional>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TLine.h>
#include <TBox.h>

namespace plot {
    class Display {
    public:
        Display(int img_size, const std::string& output_dir);
        void VisualiseInput(const analysis::Sample& sample, int event_index);
        void VisualiseLabels(const analysis::Sample& sample, int event_index);
        void VisualiseTruth(const analysis::Sample& sample, int event_index);
        void PlotLabelsLegend();
        void PlotTruthLegend();

    private:
        int img_size_;
        std::string output_dir_;
        std::vector<std::string> plane_names_ = {"U", "V", "W"};
        void SetHistogramStyle(TH2F* hist);
    };

    Display::Display(int img_size, const std::string& output_dir)
        : img_size_(img_size), output_dir_(output_dir) {}

    void Display::VisualiseInput(const analysis::Sample& sample, int event_index) {
        TChain* chain = sample.GetEventChain();
        analysis::Event event;
        event.SetBranches(chain);

        if (event_index < 0 || event_index >= chain->GetEntries()) 
            return;

        chain->GetEntry(event_index);

        std::vector<std::vector<float>*> plane_data = {&(*event.calo_pixels_u), &(*event.calo_pixels_v), &(*event.calo_pixels_w)};

        for (int plane = 0; plane < 3; ++plane) {
            std::vector<float>& data = *plane_data[plane];
            if (data.size() != static_cast<size_t>(img_size_ * img_size_)) 
                continue;

            std::string hist_name = "h_input_" + plane_names_[plane];
            std::string title = "Plane " + plane_names_[plane] + " Input (Run " + std::to_string(event.run) +
                                ", Subrun " + std::to_string(event.sub) + ", Event " + std::to_string(event.evt) + ")";
            TH2F* h_input = new TH2F(hist_name.c_str(), title.c_str(), img_size_, 0, img_size_, img_size_, 0, img_size_);

            float threshold = 8;
            float min_value = 8;
            float max_value = min_value;

            for (int i = 0; i < img_size_; ++i) {
                for (int j = 0; j < img_size_; ++j) {
                    float value = data[j * img_size_ + i];
                    if (value > threshold) {
                        h_input->SetBinContent(i + 1, j + 1, value);
                        if (value > max_value) max_value = value;
                    } else {
                        h_input->SetBinContent(i + 1, j + 1, min_value);
                    }
                }
            }

            this->SetHistogramStyle(h_input);
            h_input->SetMinimum(min_value);
            h_input->SetMaximum(max_value);

            TCanvas* c_input = new TCanvas(("c_input_" + plane_names_[plane]).c_str(), title.c_str(), 1200, 1200);
            c_input->SetFillColor(kWhite);
            c_input->SetLeftMargin(0.08);
            c_input->SetRightMargin(0.08);
            c_input->SetBottomMargin(0.08);
            c_input->SetTopMargin(0.08);
            c_input->SetLogz(1);

            h_input->Draw("COL");

            std::string output_path = output_dir_ + "/input_" + std::to_string(event.run) + "_" +
                                      std::to_string(event.sub) + "_" + std::to_string(event.evt) +
                                      "_plane_" + plane_names_[plane] + ".png";
            c_input->SaveAs(output_path.c_str());

            delete c_input;
            delete h_input;
        }
    }

    void Display::VisualiseLabels(const analysis::Sample& sample, int event_index) {
        TChain* chain = sample.GetEventChain();
        analysis::Event event;
        event.SetBranches(chain);

        if (event_index < 0 || event_index >= chain->GetEntries()) 
            return;

        chain->GetEntry(event_index);

        std::vector<std::vector<float>*> plane_data = {&(*event.reco_pixels_u), &(*event.reco_pixels_v), &(*event.reco_pixels_w)};

        for (int plane = 0; plane < 3; ++plane) {
            std::vector<float>& data = *plane_data[plane];
            if (data.size() != static_cast<size_t>(img_size_ * img_size_)) 
                continue;

            std::string hist_name = "h_label_" + plane_names_[plane];
            std::string title = "Plane " + plane_names_[plane] + " Label (Run " + std::to_string(event.run) +
                                ", Subrun " + std::to_string(event.sub) + ", Event " + std::to_string(event.evt) + ")";
            TH2F* h_label = new TH2F(hist_name.c_str(), title.c_str(), img_size_, 0, img_size_, img_size_, 0, img_size_);

            for (int i = 0; i < img_size_; ++i) {
                for (int j = 0; j < img_size_; ++j) {
                    float value = data[j * img_size_ + i];
                    h_label->SetBinContent(i + 1, j + 1, value);
                }
            }

            this->SetHistogramStyle(h_label);

            const int n_labels = 8;
            int colors[n_labels] = {kGray, kBlue, kGreen, kRed, kYellow, kMagenta, kCyan, kOrange};
            gStyle->SetNumberContours(n_labels);
            gStyle->SetPalette(n_labels, colors);
            h_label->SetMinimum(-0.5);
            h_label->SetMaximum(n_labels - 0.5);

            TCanvas* c_label = new TCanvas(("c_label_" + plane_names_[plane]).c_str(), title.c_str(), 1200, 1200);
            c_label->SetFillColor(kWhite);
            c_label->SetLeftMargin(0.08);
            c_label->SetRightMargin(0.08);
            c_label->SetBottomMargin(0.08);
            c_label->SetTopMargin(0.08);

            h_label->Draw("COL");

            std::string output_path = output_dir_ + "/label_" + std::to_string(event.run) + "_" +
                                      std::to_string(event.sub) + "_" + std::to_string(event.evt) +
                                      "_plane_" + plane_names_[plane] + ".png";
            c_label->SaveAs(output_path.c_str());

            delete c_label;
            delete h_label;
        }
    }

    void Display::VisualiseTruth(const analysis::Sample& sample, int event_index) {
        TChain* chain = sample.GetEventChain();
        analysis::Event event;
        event.SetBranches(chain);

        if (event_index < 0 || event_index >= chain->GetEntries()) 
            return;

        chain->GetEntry(event_index);

        std::vector<std::vector<float>*> plane_data = {&(*event.label_pixels_u), &(*event.label_pixels_v), &(*event.label_pixels_w)};

        for (int plane = 0; plane < 3; ++plane) {
            std::vector<float>& data = *plane_data[plane];
            if (data.size() != static_cast<size_t>(img_size_ * img_size_)) 
                continue;

            std::string hist_name = "h_truth_" + plane_names_[plane];
            std::string title = "Plane " + plane_names_[plane] + " Truth (Run " + std::to_string(event.run) +
                                ", Subrun " + std::to_string(event.sub) + ", Event " + std::to_string(event.evt) + ")";
            TH2F* h_truth = new TH2F(hist_name.c_str(), title.c_str(), img_size_, 0, img_size_, img_size_, 0, img_size_);

            for (int i = 0; i < img_size_; ++i) {
                for (int j = 0; j < img_size_; ++j) {
                    float value = data[j * img_size_ + i];
                    h_truth->SetBinContent(i + 1, j + 1, value);
                }
            }

            this->SetHistogramStyle(h_truth);

            const int n_labels = 18;
            int colors[18] = {kGray, kBlue, kRed, kGreen, kYellow, kCyan, kBlack, kMagenta, kOrange, kSpring, kTeal, kAzure, kViolet, kPink, kBlue+1, kRed+1, kGreen+1, kOrange+1};
            
            gStyle->SetNumberContours(n_labels);
            gStyle->SetPalette(n_labels, colors);
            h_truth->SetMinimum(-0.5);
            h_truth->SetMaximum(n_labels - 0.5);

            TCanvas* c_truth = new TCanvas(("c_truth_" + plane_names_[plane]).c_str(), title.c_str(), 1200, 1200);
            c_truth->SetFillColor(kWhite);
            c_truth->SetLeftMargin(0.08);
            c_truth->SetRightMargin(0.08);
            c_truth->SetBottomMargin(0.08);
            c_truth->SetTopMargin(0.08);

            c_truth->DrawFrame(0, 0, 512, 512, title.c_str());

            h_truth->Draw("COL");

            std::string output_path = output_dir_ + "/truth_" + std::to_string(event.run) + "_" +
                                      std::to_string(event.sub) + "_" + std::to_string(event.evt) +
                                      "_plane_" + plane_names_[plane] + ".png";
            c_truth->SaveAs(output_path.c_str());

            delete c_truth;
            delete h_truth;
        }
    }

    void Display::PlotLabelsLegend() {
        const std::array<std::string, 8> legend_label_names = {
            "empty", "cosmic", "MIP", "HIP", "shower", "michel", "diffuse", "invisible"
        };

        const std::array<int, 8> label_colors = {
            kGray, kBlue, kGreen, kRed, kYellow, kMagenta, kCyan, kOrange
        };

        TCanvas* c_legend = new TCanvas("legend", "Label Legend", 800, 600);
        c_legend->SetFillColor(kWhite);

        TLegend* leg = new TLegend(0.64, 0.32, 0.94, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->SetTextFont(42);
        leg->SetHeader("Labels", "C");

        for (size_t i = 0; i < legend_label_names.size(); ++i) {
            TBox* box = new TBox(0, 0, 1, 1);
            box->SetFillColor(label_colors[i]);
            box->SetLineColor(kBlack);
            leg->AddEntry(box, legend_label_names[i].c_str(), "f");
        }

        leg->Draw();

        c_legend->SaveAs((output_dir_ + "/labels_legend.png").c_str());
        delete c_legend;
    }

    void Display::PlotTruthLegend() {
        const std::array<std::string, 13> truth_label_names = {
            "empty", "cosmic", "muon", "electron", "photon",
            "proton", "pion", "neutral_pion", "kaon", "neutral_kaon", 
            "lambda", "sigma", "other"
        };

        const std::array<int, 13> label_colors = {
            kGray, kBlue, kRed, kGreen, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink
        };

        TCanvas* c_legend = new TCanvas("truth_legend", "Truth Label Legend", 800, 600);
        c_legend->SetFillColor(kWhite);

        TLegend* leg = new TLegend(0.64, 0.32, 0.94, 0.85);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        leg->SetTextFont(42);
        leg->SetHeader("Truth Labels", "C");

        for (size_t i = 0; i < truth_label_names.size(); ++i) {
            TBox* box = new TBox(0, 0, 1, 1);
            box->SetFillColor(label_colors[i]);
            box->SetLineColor(kBlack);
            leg->AddEntry(box, truth_label_names[i].c_str(), "f");
        }

        leg->Draw();

        c_legend->SaveAs((output_dir_ + "/truth_legend.png").c_str());
        delete c_legend;
    }

    void Display::SetHistogramStyle(TH2F* hist) {
        hist->GetXaxis()->SetTitle("Local Drift Time");
        hist->GetYaxis()->SetTitle("Local Wire Coordinate");
        hist->GetXaxis()->SetTitleOffset(1.0);
        hist->GetYaxis()->SetTitleOffset(1.0);
        hist->GetXaxis()->SetLabelSize(0.03);
        hist->GetYaxis()->SetLabelSize(0.03);
        hist->GetXaxis()->SetTitleSize(0.03);
        hist->GetYaxis()->SetTitleSize(0.03);
        hist->GetXaxis()->SetLabelColor(kBlack);
        hist->GetYaxis()->SetLabelColor(kBlack);
        hist->GetXaxis()->SetTitleColor(kBlack);
        hist->GetYaxis()->SetTitleColor(kBlack);
        hist->GetXaxis()->SetNdivisions(1);
        hist->GetYaxis()->SetNdivisions(1);
        hist->GetXaxis()->SetTickLength(0);
        hist->GetYaxis()->SetTickLength(0);
        hist->GetXaxis()->CenterTitle();
        hist->GetYaxis()->CenterTitle();
        hist->SetStats(0);
    }
}

#endif