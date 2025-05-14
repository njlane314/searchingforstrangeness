#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory> // For std::shared_ptr or std::unique_ptr if managing TH1s

// Forward declare ROOT classes to reduce header dependencies if possible,
// otherwise include them directly. For a library-like header, forward is better.
// However, for direct use and simplicity, includes are fine.
#include "TH1.h" // Base for TH1D
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"

#include "EventCategories.h" // From previous generation

namespace AnalysisFramework {

// Configuration structure for plotting options
struct PlotterConfig {
    // General plot appearance
    std::string output_plot_name_prefix = "hist_plot"; // Will append .png, .pdf, etc.
    std::string plot_title_main = "";       // Main title part, selection info added
    std::string x_axis_title = "X-axis";
    std::string y_axis_title = "Events";
    bool show_data_mc_ratio = true;
    double ratio_plot_y_min = 0.0;
    double ratio_plot_y_max = 2.0;
    bool log_y = false;
    double y_axis_max_scale = 1.4; // Scale factor for Y-axis maximum

    // Data
    bool show_data = true;
    int data_marker_style = 20; // kFullCircle
    double data_marker_size = 1.0;
    int data_color = 1; // kBlack
    std::string data_legend_label = "Data";

    // MC
    bool mc_stacked = true;
    bool show_total_mc_line = true;
    int total_mc_line_color = 1; // kBlack
    int total_mc_line_style = 1; // kSolid
    double total_mc_line_width = 2.0;
    bool show_mc_error_band = true;
    int mc_error_band_color = 15; // kGray
    int mc_error_band_style = 3004; // Hatch style for error band
    std::string total_mc_legend_label = "Total Pred.";

    // Signal (if any)
    bool show_signal_on_top = true; // If true, signal is added on top of stacked backgrounds
    std::string signal_default_legend_label = "Signal";
    int signal_line_color = 632; // kRed
    int signal_line_style = 2;   // kDashed
    double signal_line_width = 2.5;
    double signal_scale_factor = 1.0; // To scale signal independently if needed for visualization

    // Legend
    bool draw_legend = true;
    int legend_n_columns = 2;
    double legend_x1_ndc = 0.55; double legend_y1_ndc = 0.65;
    double legend_x2_ndc = 0.92; double legend_y2_ndc = 0.92;
    bool show_counts_in_legend = true;
    double legend_text_size = 0.035;


    // Labels & Text
    std::string experiment_preliminary_label = "My Experiment Preliminary";
    std::string pot_label_override = ""; // If empty, will be auto-generated from data_pot
    std::string extra_text_left = "";    // Additional text on the left
    std::string extra_text_right = "";   // Additional text on the right
    double top_label_text_size = 0.04;


    // Chi2
    bool show_chi_square = false;

    PlotterConfig() = default;
};


class HistogramPlotter {
public:
    HistogramPlotter(EventCategories event_cats);
    ~HistogramPlotter(); // To handle TObject cleanup if needed

    // Main plot function
    // Takes pre-filled histograms.
    // mc_background_hists: map from category_id (int) to TH1D* (owned by caller or managed by this class)
    // total_mc_pred_with_syst: TH1D* where bin content is CV and bin error is total systematic uncertainty
    void PlotHistograms(
        TH1D* data_hist,                                 // Can be nullptr if blinded
        TH1D* ext_hist,                                  // Can be nullptr
        const std::map<int, TH1D*>& mc_background_hists,
        TH1D* signal_hist,                               // Can be nullptr (if no dedicated signal to overlay)
        TH1D* total_mc_pred_with_syst,                   // Total MC prediction with systematic uncertainties
        double data_pot,                                 // For POT label
        const PlotterConfig& config = PlotterConfig(),   // Plotting options
        const std::string& specific_selection_title = "" // Title for this specific plot
    ) const;

private:
    EventCategories event_categories_;

    // Helper methods for drawing various components of the plot
    void SetupCanvasPads(TCanvas*& canvas, TPad*& main_pad, TPad*& ratio_pad, bool do_ratio) const;
    void DrawStackedMCHistograms(TPad* pad, THStack* mc_stack, const std::map<int, TH1D*>& mc_hists, TLegend* legend, const PlotterConfig& config, double& total_mc_sum_for_legend) const;
    void DrawTotalMCWithErrorBand(TPad* pad, TH1D* total_mc_hist, TLegend* legend, const PlotterConfig& config) const;
    void DrawDataHistogram(TPad* pad, TH1D* data_hist, TLegend* legend, const PlotterConfig& config) const;
    void DrawSignalHistogram(TPad* pad, TH1D* signal_hist, TH1D* total_mc_for_stacking_signal, TLegend* legend, const PlotterConfig& config) const;
    void DrawRatioPlot(TPad* pad, TH1D* data_hist, TH1D* total_mc_hist, TH1D* signal_hist_optional, const PlotterConfig& config) const;
    void AddPlotLabels(TCanvas* canvas, TPad* main_pad, double data_pot, const PlotterConfig& config, const std::string& specific_selection_title, double chi2_val = -1, int ndf = -1) const;
    std::string GetPOTLabel(double data_pot, const PlotterConfig& config) const;
    double CalculateChiSquare(TH1D* data_hist, TH1D* total_mc_hist, int& ndf) const; // Simple chi2
    void StyleAxis(TH1* hist, const PlotterConfig& config, bool is_ratio_plot = false) const;
};

} // namespace AnalysisFramework


#include "HistogramPlotter.h" // Assuming the header is in the same directory or include path

#include "TStyle.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TF1.h" // For fitting line in ratio plot if needed

#include <iostream>     // For std::cout, std::cerr (debugging or info)
#include <iomanip>      // For std::fixed, std::setprecision
#include <numeric>      // For std::accumulate
#include <algorithm>    // For std::sort, std::max_element

namespace AnalysisFramework {

HistogramPlotter::HistogramPlotter(EventCategories event_cats) : event_categories_(std::move(event_cats)) {
    // Set some global ROOT styles if desired, or do it in main.cpp
    // gROOT->SetStyle("Plain"); // Example
    // gStyle->SetOptStat(0);    // No stats box by default
    // gStyle->SetErrorX(0.5);   // For data error bars
}

HistogramPlotter::~HistogramPlotter() {
    // Cleanup of ROOT objects created by the plotter if they are not managed by TCanvas/TPad
}

std::string HistogramPlotter::GetPOTLabel(double data_pot, const PlotterConfig& config) const {
    if (!config.pot_label_override.empty()) {
        return config.pot_label_override;
    }
    if (data_pot <= 0) {
        return "MC (Unscaled)";
    }
    char pot_buffer[100];
    sprintf(pot_buffer, "%.2e POT", data_pot);
    std::string pot_str = pot_buffer;
    size_t e_pos = pot_str.find('e');
    if (e_pos != std::string::npos) {
        std::string base = pot_str.substr(0, e_pos);
        std::string exponent = pot_str.substr(e_pos + 1);
        // Remove leading '+' from exponent if present
        if (!exponent.empty() && exponent[0] == '+') {
            exponent = exponent.substr(1);
        }
        // Remove leading zeros from exponent unless it's just "0"
        size_t first_digit = exponent.find_first_not_of('0');
        if (std::string::npos == first_digit) { // all zeros
             if (exponent.length() > 1) exponent = "0";
        } else {
            exponent = exponent.substr(first_digit);
        }
        if (exponent.empty()) exponent = "0";


        return base + " #times 10^{" + exponent + "} POT";
    }
    return pot_str; // Fallback
}


void HistogramPlotter::SetupCanvasPads(TCanvas*& canvas, TPad*& main_pad, TPad*& ratio_pad, bool do_ratio) const {
    canvas = new TCanvas( ("c_" + std::to_string(rand()%100000)).c_str(), "Histogram Plot", 800, (do_ratio ? 800 : 600));
    canvas->cd();

    if (do_ratio) {
        main_pad = new TPad("main_pad", "Main Plot", 0.0, 0.3, 1.0, 1.0);
        main_pad->SetBottomMargin(0.02); // Reduce space for x-axis label, it will be on ratio
        main_pad->Draw();

        ratio_pad = new TPad("ratio_pad", "Ratio Plot", 0.0, 0.0, 1.0, 0.3);
        ratio_pad->SetTopMargin(0.05);
        ratio_pad->SetBottomMargin(0.35); // More space for x-axis label
        ratio_pad->SetGridy();
        ratio_pad->Draw();
    } else {
        main_pad = new TPad("main_pad", "Main Plot", 0.0, 0.0, 1.0, 1.0);
        main_pad->SetBottomMargin(0.12);
        main_pad->Draw();
        ratio_pad = nullptr;
    }
}

void HistogramPlotter::StyleAxis(TH1* hist, const PlotterConfig& config, bool is_ratio_plot) const {
    if (!hist) return;
    hist->GetXaxis()->SetTitle(config.x_axis_title.c_str());
    hist->GetYaxis()->SetTitle(config.y_axis_title.c_str());

    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.04);
    hist->GetYaxis()->SetTitleOffset(1.2);


    if (is_ratio_plot) {
        hist->GetYaxis()->SetTitle("Data/MC");
        hist->GetYaxis()->SetNdivisions(505); // Fewer divisions for ratio plot
        hist->GetYaxis()->SetRangeUser(config.ratio_plot_y_min, config.ratio_plot_y_max);
        hist->GetXaxis()->SetTitleSize(0.12);
        hist->GetXaxis()->SetLabelSize(0.11);
        hist->GetYaxis()->SetTitleSize(0.10);
        hist->GetYaxis()->SetLabelSize(0.10);
        hist->GetYaxis()->SetTitleOffset(0.5);
        hist->GetYaxis()->CenterTitle();
    }
}


void HistogramPlotter::DrawStackedMCHistograms(TPad* pad, THStack* mc_stack, const std::map<int, TH1D*>& mc_hists, TLegend* legend, const PlotterConfig& config, double& total_mc_sum_for_legend) const {
    pad->cd();
    total_mc_sum_for_legend = 0.0;

    // Order histograms for stacking (optional, but can make legend nicer)
    // For now, iterate as is. Could sort by integral or predefined order.
    std::vector<int> category_ids;
    for(const auto& pair : mc_hists) {
        if (pair.second && pair.second->Integral() > 1e-6) { // Only consider non-empty
            category_ids.push_back(pair.first);
        }
    }
    // Example sorting: by category ID (can be customized)
    std::sort(category_ids.begin(), category_ids.end());

    for (int cat_id : category_ids) {
        TH1D* hist = mc_hists.at(cat_id);
        if (hist) {
            CategoryDisplayInfo cat_info = event_categories_.GetCategoryInfo(cat_id);
            hist->SetFillColor(cat_info.color);
            hist->SetLineColor(kBlack);
            hist->SetLineWidth(1);
            mc_stack->Add(hist);
            double integral = hist->Integral();
            total_mc_sum_for_legend += integral;
            if (legend && config.draw_legend) {
                std::string label = config.show_counts_in_legend ?
                                    cat_info.short_label + " (" + std::to_string((int)round(integral)) + ")" :
                                    cat_info.short_label;
                legend->AddEntry(hist, label.c_str(), "f");
            }
        }
    }
    mc_stack->Draw("HIST SAME"); // "HIST" draws filled histograms
}

void HistogramPlotter::DrawTotalMCWithErrorBand(TPad* pad, TH1D* total_mc_hist, TLegend* legend, const PlotterConfig& config) const {
    if (!total_mc_hist) return;
    pad->cd();

    if (config.show_mc_error_band) {
        total_mc_hist->SetFillColorAlpha(config.mc_error_band_color, 0.5); // Semi-transparent
        total_mc_hist->SetFillStyle(config.mc_error_band_style);
        total_mc_hist->SetMarkerSize(0); // No markers for error band
        total_mc_hist->SetLineColor(config.mc_error_band_color); // Edge color for hatch
        total_mc_hist->Draw("E2 SAME"); // "E2" draws error band
        if (legend && config.draw_legend) {
            legend->AddEntry(total_mc_hist, "Stat #oplus Syst Unc.", "f");
        }
    }

    if (config.show_total_mc_line) {
        TH1D* total_line = (TH1D*)total_mc_hist->Clone("total_mc_line");
        total_line->SetLineColor(config.total_mc_line_color);
        total_line->SetLineStyle(config.total_mc_line_style);
        total_line->SetLineWidth(config.total_mc_line_width);
        total_line->SetFillColor(0); // No fill for the line
        total_line->Draw("HIST SAME L"); // "L" for line
        // Legend entry for total line is usually handled by the stack or if no stack
        if (legend && config.draw_legend && !config.mc_stacked && config.show_total_mc_line) {
             std::string label = config.show_counts_in_legend ?
                                config.total_mc_legend_label + " (" + std::to_string((int)round(total_line->Integral())) + ")" :
                                config.total_mc_legend_label;
            legend->AddEntry(total_line, label.c_str(), "l");
        }
    }
}

void HistogramPlotter::DrawDataHistogram(TPad* pad, TH1D* data_hist, TLegend* legend, const PlotterConfig& config) const {
    if (!data_hist || !config.show_data) return;
    pad->cd();
    data_hist->SetMarkerStyle(config.data_marker_style);
    data_hist->SetMarkerSize(config.data_marker_size);
    data_hist->SetMarkerColor(config.data_color);
    data_hist->SetLineColor(config.data_color);
    data_hist->SetLineWidth(2);
    data_hist->Draw("E1 X0 SAME"); // "E1" for error bars, "X0" to suppress horizontal error bars
    if (legend && config.draw_legend) {
        std::string label = config.show_counts_in_legend ?
                            config.data_legend_label + " (" + std::to_string((int)round(data_hist->Integral())) + ")" :
                            config.data_legend_label;
        legend->AddEntry(data_hist, label.c_str(), "lep");
    }
}

void HistogramPlotter::DrawSignalHistogram(TPad* pad, TH1D* signal_hist, TH1D* total_mc_for_stacking_signal, TLegend* legend, const PlotterConfig& config) const {
    if (!signal_hist || !config.show_signal || signal_hist->Integral() < 1e-6) return;
    pad->cd();

    TH1D* hist_to_draw = (TH1D*)signal_hist->Clone("signal_to_draw");
    hist_to_draw->Scale(config.signal_scale_factor);

    hist_to_draw->SetLineColor(config.signal_line_color);
    hist_to_draw->SetLineStyle(config.signal_line_style);
    hist_to_draw->SetLineWidth(config.signal_line_width);
    hist_to_draw->SetFillColor(0); // No fill for signal line usually

    if (config.mc_stacked && config.show_signal_on_top && total_mc_for_stacking_signal) {
        // Add signal on top of the total background stack for visualization
        TH1D* stacked_signal = (TH1D*)total_mc_for_stacking_signal->Clone("stacked_signal_temp");
        stacked_signal->Add(hist_to_draw);
        stacked_signal->SetLineColor(config.signal_line_color);
        stacked_signal->SetLineStyle(config.signal_line_style);
        stacked_signal->SetLineWidth(config.signal_line_width);
        stacked_signal->SetFillColor(0);
        stacked_signal->Draw("HIST SAME L");
        // Optionally draw the signal itself non-stacked as well for clarity if desired
        // hist_to_draw->Draw("HIST SAME L");
    } else {
        hist_to_draw->Draw("HIST SAME L");
    }

    if (legend && config.draw_legend) {
        std::string label = config.signal_default_legend_label;
        if (config.show_counts_in_legend) {
            label += " (" + std::to_string((int)round(hist_to_draw->Integral())) + ")";
            if (config.signal_scale_factor != 1.0) {
                 label += " #times" + std::to_string(config.signal_scale_factor);
            }
        }
        legend->AddEntry(hist_to_draw, label.c_str(), "l");
    }
    // delete hist_to_draw; // If cloned
}


void HistogramPlotter::DrawRatioPlot(TPad* pad, TH1D* data_hist, TH1D* total_mc_hist, TH1D* signal_hist_optional, const PlotterConfig& config) const {
    if (!data_hist || !total_mc_hist) return;
    pad->cd();

    TH1D* h_ratio = (TH1D*)data_hist->Clone("h_ratio_data_mc");
    h_ratio->Divide(total_mc_hist);
    StyleAxis(h_ratio, config, true); // Styles Y axis for ratio
    h_ratio->Draw("E1 X0 P"); // P for points, E1 for errors

    // MC uncertainty band on ratio (ratio of (MC+-err)/MC_CV = 1 +- err/MC_CV)
    TH1D* h_mc_ratio_band = (TH1D*)total_mc_hist->Clone("h_mc_ratio_band");
    for (int i = 1; i <= total_mc_hist->GetNbinsX(); ++i) {
        if (total_mc_hist->GetBinContent(i) > 1e-9) {
            h_mc_ratio_band->SetBinContent(i, 1.0);
            h_mc_ratio_band->SetBinError(i, total_mc_hist->GetBinError(i) / total_mc_hist->GetBinContent(i));
        } else {
            h_mc_ratio_band->SetBinContent(i, 1.0);
            h_mc_ratio_band->SetBinError(i, 0.0); // Or a large error
        }
    }
    h_mc_ratio_band->SetFillColorAlpha(config.mc_error_band_color, 0.5);
    h_mc_ratio_band->SetFillStyle(config.mc_error_band_style);
    h_mc_ratio_band->SetMarkerSize(0);
    h_mc_ratio_band->SetLineColor(config.mc_error_band_color);
    h_mc_ratio_band->Draw("E2 SAME");

    // Line at 1 for reference
    TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0, h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw("SAME");

    // Optional: Signal / MC on ratio plot
    if (config.show_signal_on_top && signal_hist_optional && signal_hist_optional->Integral() > 1e-6) {
        TH1D* h_signal_ratio = (TH1D*)signal_hist_optional->Clone("h_signal_ratio");
        h_signal_ratio->Scale(config.signal_scale_factor);
        h_signal_ratio->Add(total_mc_hist); // (Signal+MC)/MC
        h_signal_ratio->Divide(total_mc_hist);
        h_signal_ratio->SetLineColor(config.signal_line_color);
        h_signal_ratio->SetLineStyle(config.signal_line_style);
        h_signal_ratio->SetLineWidth(config.signal_line_width);
        h_signal_ratio->SetFillColor(0);
        h_signal_ratio->Draw("HIST SAME L");
    }
}

void HistogramPlotter::AddPlotLabels(TCanvas* canvas, TPad* main_pad, double data_pot, const PlotterConfig& config,
                                     const std::string& specific_selection_title, double chi2_val, int ndf) const {
    main_pad->cd();
    float left_x = 0.15; // NDC
    float right_x = 0.92; // NDC
    float top_y = 0.90;  // NDC
    float text_size = config.top_label_text_size;

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42); // Helvetica
    latex.SetTextSize(text_size);

    // Experiment Label + POT
    std::string exp_pot_label = config.experiment_preliminary_label;
    if (data_pot > 0) {
        exp_pot_label += " (" + GetPOTLabel(data_pot, config) + ")";
    }
    latex.SetTextAlign(11); // Left-top aligned
    latex.DrawLatex(left_x, top_y, exp_pot_label.c_str());

    // Selection Title (right aligned, below experiment label)
    std::string full_title = config.plot_title_main;
    if (!specific_selection_title.empty()){
        if(!full_title.empty()) full_title += ": ";
        full_title += specific_selection_title;
    }
    if (!full_title.empty()) {
        latex.SetTextAlign(31); // Right-top aligned
        latex.DrawLatex(right_x, top_y, full_title.c_str());
    }
    
    float current_y_pos_left = top_y - text_size * 1.5;
    float current_y_pos_right = top_y - text_size * 1.5;


    // Chi-square
    if (config.show_chi_square && chi2_val >= 0 && ndf > 0) {
        char chi2_buffer[100];
        sprintf(chi2_buffer, "#chi^{2}/ndf = %.1f/%d", chi2_val, ndf);
        latex.SetTextAlign(11); // Left-top
        latex.DrawLatex(left_x, current_y_pos_left, chi2_buffer);
        current_y_pos_left -= text_size * 1.2;
    }

    // Extra Text
    if (!config.extra_text_left.empty()) {
        latex.SetTextAlign(11);
        latex.DrawLatex(left_x, current_y_pos_left, config.extra_text_left.c_str());
    }
    if (!config.extra_text_right.empty()) {
        latex.SetTextAlign(31);
        latex.DrawLatex(right_x, current_y_pos_right, config.extra_text_right.c_str());
    }
}

double HistogramPlotter::CalculateChiSquare(TH1D* data_hist, TH1D* total_mc_hist, int& ndf) const {
    if (!data_hist || !total_mc_hist || data_hist->GetNbinsX() != total_mc_hist->GetNbinsX()) {
        ndf = 0;
        return -1.0;
    }
    double chi2 = 0.0;
    ndf = 0;
    for (int i = 1; i <= data_hist->GetNbinsX(); ++i) {
        double data_val = data_hist->GetBinContent(i);
        double mc_val = total_mc_hist->GetBinContent(i);
        double err_data_sq = data_hist->GetBinError(i) * data_hist->GetBinError(i);
        double err_mc_sq = total_mc_hist->GetBinError(i) * total_mc_hist->GetBinError(i); // This is stat+syst from input

        if (mc_val > 1e-6 && (err_data_sq + err_mc_sq) > 1e-9) { // Avoid division by zero and bins with no MC
            chi2 += std::pow(data_val - mc_val, 2) / (err_data_sq + err_mc_sq); // Simple chi2, no correlations
            ndf++;
        }
    }
    return chi2;
}


void HistogramPlotter::PlotHistograms(
    TH1D* data_hist,
    TH1D* ext_hist,
    const std::map<int, TH1D*>& mc_background_hists,
    TH1D* signal_hist,
    TH1D* total_mc_pred_with_syst,
    double data_pot,
    const PlotterConfig& config,
    const std::string& specific_selection_title
) const {
    TCanvas* canvas = nullptr;
    TPad* main_pad = nullptr;
    TPad* ratio_pad = nullptr;

    SetupCanvasPads(canvas, main_pad, ratio_pad, config.show_data_mc_ratio);

    main_pad->cd();
    if (config.log_y) {
        main_pad->SetLogy();
    }

    // Determine overall Y range for the main plot
    double max_y = 0;
    if (data_hist && config.show_data) {
        for (int i = 1; i <= data_hist->GetNbinsX(); ++i) {
            max_y = std::max(max_y, data_hist->GetBinContent(i) + data_hist->GetBinError(i));
        }
    }
    if (total_mc_pred_with_syst) {
         for (int i = 1; i <= total_mc_pred_with_syst->GetNbinsX(); ++i) {
            max_y = std::max(max_y, total_mc_pred_with_syst->GetBinContent(i) + total_mc_pred_with_syst->GetBinError(i));
        }
    }
    if (signal_hist && config.show_signal_on_top) {
        TH1D* temp_total_plus_signal = (TH1D*)total_mc_pred_with_syst->Clone();
        temp_total_plus_signal->Add(signal_hist, config.signal_scale_factor);
         for (int i = 1; i <= temp_total_plus_signal->GetNbinsX(); ++i) {
            max_y = std::max(max_y, temp_total_plus_signal->GetBinContent(i));
        }
        delete temp_total_plus_signal;
    } else if (signal_hist && !config.show_signal_on_top) {
         for (int i = 1; i <= signal_hist->GetNbinsX(); ++i) {
            max_y = std::max(max_y, signal_hist->GetBinContent(i) * config.signal_scale_factor);
        }
    }


    // Create a frame histogram to set axes ranges correctly before drawing stacks/graphs
    TH1* h_frame = nullptr;
    if (data_hist) h_frame = data_hist;
    else if (total_mc_pred_with_syst) h_frame = total_mc_pred_with_syst;
    else if (!mc_background_hists.empty()) h_frame = mc_background_hists.begin()->second;
    else if (ext_hist) h_frame = ext_hist;

    if (h_frame) {
        TH1* frame_clone = (TH1*)h_frame->Clone("frame_clone");
        StyleAxis(frame_clone, config, false);
        frame_clone->GetYaxis()->SetRangeUser(config.log_y ? 0.1 : 0, max_y * config.y_axis_max_scale);
        frame_clone->Draw("AXIS"); // Draw only axes
    }


    TLegend* legend = nullptr;
    if (config.draw_legend) {
        legend = new TLegend(config.legend_x1_ndc, config.legend_y1_ndc, config.legend_x2_ndc, config.legend_y2_ndc);
        legend->SetNColumns(config.legend_n_columns);
        legend->SetFillStyle(0); // Transparent
        legend->SetBorderSize(0);
        legend->SetTextFont(42);
        legend->SetTextSize(config.legend_text_size);
    }

    THStack* mc_stack = nullptr;
    double total_mc_sum_for_legend = 0; // Used if stack is drawn and legend needs total sum
    if (config.mc_stacked) {
        mc_stack = new THStack("mc_stack", "");
        // Add EXT to the stack if it exists
        std::map<int, TH1D*> hists_for_stacking = mc_background_hists;
        if (ext_hist) {
            hists_for_stacking[-1] = ext_hist; // Use a temporary ID for EXT for ordering/coloring
                                               // Or handle EXT separately if it shouldn't be stacked
        }
        DrawStackedMCHistograms(main_pad, mc_stack, hists_for_stacking, legend, config, total_mc_sum_for_legend);
        if (h_frame == nullptr && mc_stack->GetHistogram()) { // If no data/total, use stack to set axis titles
             StyleAxis(mc_stack->GetHistogram(), config, false);
             mc_stack->GetHistogram()->GetYaxis()->SetRangeUser(config.log_y ? 0.1 : 0, max_y * config.y_axis_max_scale);
        }
    } else {
        // Plot MC categories unstacked
        for (const auto& pair : mc_background_hists) {
            // DrawTotalMCWithErrorBand(main_pad, pair.second, legend, config); // This would draw error band for each
            pair.second->SetLineColor(event_categories_.GetColor(pair.first));
            pair.second->SetLineWidth(2);
            pair.second->Draw("HIST SAME L");
            if(legend) legend->AddEntry(pair.second, event_categories_.GetShortLabel(pair.first).c_str(), "l");
        }
        if(ext_hist){
            ext_hist->SetLineColor(event_categories_.GetColor(1)); // Assuming cat 1 is EXT
            ext_hist->SetLineWidth(2);
            ext_hist->Draw("HIST SAME L");
            if(legend) legend->AddEntry(ext_hist, event_categories_.GetShortLabel(1).c_str(), "l");
        }
    }

    // Draw total MC prediction with its systematic error band
    if (total_mc_pred_with_syst) {
        DrawTotalMCWithErrorBand(main_pad, total_mc_pred_with_syst, legend, config);
        if (config.mc_stacked && legend && config.draw_legend && config.show_total_mc_line) { // Add total from stack to legend
             std::string label = config.show_counts_in_legend ?
                                config.total_mc_legend_label + " (" + std::to_string((int)round(total_mc_sum_for_legend)) + ")" :
                                config.total_mc_legend_label;
            legend->AddEntry(total_mc_pred_with_syst, label.c_str(), "l"); // "l" for line, "f" if also showing band in legend
        }
    }
    
    // Draw signal histogram
    // If stacked, total_mc_pred_with_syst is used as the base for stacking signal on top.
    // If not stacked, signal is drawn independently.
    TH1D* base_for_signal = (config.mc_stacked && config.show_signal_on_top) ? total_mc_pred_with_syst : nullptr;
    DrawSignalHistogram(main_pad, signal_hist, base_for_signal, legend, config);


    // Draw data on top
    DrawDataHistogram(main_pad, data_hist, legend, config);

    if (legend && config.draw_legend) {
        legend->Draw();
    }
    main_pad->RedrawAxis(); // Redraw axis on top of filled histograms

    // Calculate Chi2
    double chi2 = -1; int ndf = 0;
    if (config.show_chi_square && data_hist && total_mc_pred_with_syst) {
        chi2 = CalculateChiSquare(data_hist, total_mc_pred_with_syst, ndf);
    }

    // Add top labels (Experiment, POT, Title, Chi2)
    AddPlotLabels(canvas, main_pad, data_pot, config, specific_selection_title, chi2, ndf);


    // Ratio Plot
    if (config.show_data_mc_ratio && ratio_pad && data_hist && total_mc_pred_with_syst) {
        DrawRatioPlot(ratio_pad, data_hist, total_mc_pred_with_syst, (config.show_signal_on_top ? signal_hist : nullptr), config);
        ratio_pad->RedrawAxis();
    }

    // Save canvas
    canvas->SaveAs((config.output_plot_name_prefix + ".png").c_str());
    canvas->SaveAs((config.output_plot_name_prefix + ".pdf").c_str());
    // canvas->SaveAs((config.output_plot_name_prefix + ".C").c_str()); // Optional

    // Cleanup
    delete legend;
    delete mc_stack; // THStack owns the histograms added to it if not told otherwise
                     // If hists are managed outside, ensure THStack doesn't delete them.
                     // For safety, if mc_hists are pointers to objects managed elsewhere,
                     // one might clone them before adding to THStack.
                     // Here, assuming THStack takes ownership or they are temporary.
    // delete canvas; // Caller should manage canvas if it's returned or used further
}


double HistogramPlotter::CalculateChiSquare(TH1D* data_hist, TH1D* total_mc_hist, int& ndf) const {
    double chi2_val = 0.0;
    int n_bins_for_chi2 = 0;
    for (int i = 1; i <= data_hist->GetNbinsX(); ++i) {
        double obs = data_hist->GetBinContent(i);
        double exp = total_mc_hist->GetBinContent(i);
        // Use data stat error for data, total (stat+syst) error for MC prediction
        double err_obs_sq = data_hist->GetBinError(i) * data_hist->GetBinError(i);
        double err_exp_sq = total_mc_hist->GetBinError(i) * total_mc_hist->GetBinError(i);

        if (exp > 1e-6) { // Only include bins with some MC expectation
            double denominator = err_obs_sq + err_exp_sq; // Assuming data and MC errors are uncorrelated for this simple chi2
            if (denominator > 1e-9) {
                chi2_val += std::pow(obs - exp, 2) / denominator;
                n_bins_for_chi2++;
            }
        }
    }
    ndf = n_bins_for_chi2; // Simplistic NDF, assumes no fit parameters constrained by this data
    return chi2_val;
}


} // namespace AnalysisFramework
