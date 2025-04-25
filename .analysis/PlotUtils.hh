#pragma once

#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TGraph.h"
#include "TH1D.h"

std::string get_legend_title(double bnb_pot) {
    std::stringstream temp_ss;
    temp_ss << std::scientific << std::setprecision(3) << bnb_pot;
    std::string pot_str = temp_ss.str();
    size_t e_pos = pot_str.find('e');
    std::string pot_digits_str = pot_str.substr(0, e_pos);
    pot_str.erase(0, e_pos + 1u);
    if (pot_str.front() == '+') pot_str.erase(0, 1u);
    std::string legend_title = "MicroBooNE " + pot_digits_str + " #times 10^{" + pot_str + "} POT, INTERNAL";
    return legend_title;
}

void dump_1d_histogram(const std::string& hist_col_prefix, const TH1D& hist,
    std::map<std::string, std::vector<std::string>>& pgf_plots_hist_table,
    bool include_yerror, bool include_x_coords)
{
    const std::string bin_col_name = "bin";
    if (include_x_coords) {
        pgf_plots_hist_table[bin_col_name] = std::vector<std::string>();
    }
    std::string y_col_name = hist_col_prefix;
    pgf_plots_hist_table[y_col_name] = std::vector<std::string>();
    std::string yerror_col_name;
    if (include_yerror) {
        yerror_col_name = hist_col_prefix + "_error";
        pgf_plots_hist_table[yerror_col_name] = std::vector<std::string>();
    }
    int overflow_bin_idx = hist.GetNbinsX() + 1;
    for (int b = 1; b <= overflow_bin_idx; ++b) {
        double y = hist.GetBinContent(b);
        pgf_plots_hist_table.at(y_col_name).push_back(std::to_string(y));
        if (include_yerror) {
            double yerror = hist.GetBinError(b);
            pgf_plots_hist_table.at(yerror_col_name).push_back(std::to_string(yerror));
        }
        if (include_x_coords) {
            pgf_plots_hist_table.at(bin_col_name).push_back(std::to_string(b));
            double low_edge = hist.GetBinLowEdge(b);
            double half_width = 0.5 * hist.GetBinWidth(b);
            std::string x_col_name = hist_col_prefix;
            if (!hist_col_prefix.empty()) x_col_name += '_';
            x_col_name += 'x';
            auto end = pgf_plots_hist_table.end();
            auto iter = pgf_plots_hist_table.find(x_col_name);
            if (iter == end) {
                pgf_plots_hist_table[x_col_name] = std::vector<std::string>{std::to_string(low_edge)};
            }
            else iter->second.push_back(std::to_string(low_edge));
            std::string x_hw_col_name = x_col_name + "_halfwidth";
            auto end2 = pgf_plots_hist_table.end();
            auto iter2 = pgf_plots_hist_table.find(x_hw_col_name);
            if (iter2 == end2) {
                pgf_plots_hist_table[x_hw_col_name] = std::vector<std::string>{std::to_string(half_width)};
            }
            else iter2->second.push_back(std::to_string(half_width));
        }
    }
}

void dump_tgraph(const std::string& col_prefix, const TGraph& graph,
    std::map<std::string, std::vector<std::string>>& pgf_plots_table)
{
    std::string x_col_name = col_prefix + "_x";
    pgf_plots_table[x_col_name] = std::vector<std::string>();
    std::string y_col_name = col_prefix + "_y";
    pgf_plots_hist_table[y_col_name] = std::vector<std::string>();
    int num_points = graph.GetN();
    for (int p = 0; p < num_points; ++p) {
        double x = graph.GetX()[p];
        double y = graph.GetY()[p];
        pgf_plots_table.at(x_col_name).push_back(std::to_string(x));
        pgf_plots_table.at(y_col_name).push_back(std::to_string(y));
    }
}

void write_pgfplots_file(const std::string& out_filename,
    std::map<std::string, std::vector<std::string>>& pgfplots_table)
{
    std::ofstream out_file(out_filename);
    for (const auto& pair : pgfplots_table) {
        const std::string& col_name = pair.first;
        out_file << "  " << col_name;
    }
    out_file << '\n';
    size_t num_rows = pgfplots_table.cbegin()->second.size();
    for (size_t r = 0u; r < num_rows; ++r) {
        for (const auto& pair : pgfplots_table) {
            out_file << "  " << pair.second.at(r);
        }
        out_file << '\n';
    }
}