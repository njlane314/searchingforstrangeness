#ifndef EVENT_CATEGORY_HH
#define EVENT_CATEGORY_HH

#include <map>
#include <string>
#include <memory>
#include "TH1.h"

enum class EventCategory {
    k_data = 0,
    k_external = 1,
    k_nc = 2,
    k_nu_e_cc = 3,
    k_nu_mu_cc_with_strange = 4,
    k_nu_mu_cc_with_protons = 5,
    k_nu_mu_cc_with_pions = 6,
    k_nu_mu_cc_with_protons_pions = 7,
    k_nu_mu_cc_other = 8,
    k_other = 9
};

class EventCategoryInterpreter {
public:
    static const EventCategoryInterpreter& Instance() {
        static std::unique_ptr<EventCategoryInterpreter> the_instance(new EventCategoryInterpreter());
        return *the_instance;
    }

    const std::map<EventCategory, std::string>& label_map() const {
        return event_category_to_label_map_;
    }

    std::string label(EventCategory ec) const {
        return event_category_to_label_map_.at(ec);
    }

    int color_code(EventCategory ec) const {
        return event_category_to_color_map_.at(ec);
    }

    void set_mc_histogram_style(EventCategory ec, TH1* mc_hist) const {
        int color = color_code(ec);
        mc_hist->SetFillColor(color);
        mc_hist->SetLineColor(color);
        mc_hist->SetStats(false);
    }

    void set_external_histogram_style(TH1* ext_hist) const {
        ext_hist->SetFillColor(28);
        ext_hist->SetLineColor(28);
        ext_hist->SetLineWidth(2);
        ext_hist->SetFillStyle(3005);
        ext_hist->SetStats(false);
    }

    void set_data_histogram_style(TH1* data_hist) const {
        data_hist->SetLineColor(kBlack);
        data_hist->SetLineWidth(3);
        data_hist->SetMarkerStyle(kFullCircle);
        data_hist->SetMarkerSize(0.8);
        data_hist->SetStats(false);
        data_hist->GetXaxis()->SetTitleOffset(0.0);
        data_hist->GetXaxis()->SetTitleSize(0.0);
        data_hist->GetYaxis()->SetTitleSize(0.05);
        data_hist->GetYaxis()->CenterTitle(true);
        data_hist->GetXaxis()->SetLabelSize(0.0);
        data_hist->SetMinimum(1e-3);
    }

private:
    EventCategoryInterpreter() {
        event_category_to_label_map_ = {
            { EventCategory::k_data, "Data" },
            { EventCategory::k_external, "External" },
            { EventCategory::k_nc, "NC" },
            { EventCategory::k_nu_e_cc, "#nu_{e} CC" },
            { EventCategory::k_nu_mu_cc_with_strange, "#nu_{#mu} CC with Strange" },
            { EventCategory::k_nu_mu_cc_with_protons, "#nu_{#mu} CC with Protons" },
            { EventCategory::k_nu_mu_cc_with_pions, "#nu_{#mu} CC with Pions" },
            { EventCategory::k_nu_mu_cc_with_protons_pions, "#nu_{#mu} CC with Protons and Pions" },
            { EventCategory::k_nu_mu_cc_other, "Other #nu_{#mu} CC" },
            { EventCategory::k_other, "Other" }
        };
        event_category_to_color_map_ = {
            { EventCategory::k_data, kBlack },
            { EventCategory::k_external, kGray },
            { EventCategory::k_nc, kOrange },
            { EventCategory::k_nu_e_cc, kViolet },
            { EventCategory::k_nu_mu_cc_with_strange, kBlue },
            { EventCategory::k_nu_mu_cc_with_protons, kGreen },
            { EventCategory::k_nu_mu_cc_with_pions, kCyan },
            { EventCategory::k_nu_mu_cc_with_protons_pions, kMagenta },
            { EventCategory::k_nu_mu_cc_other, kAzure },
            { EventCategory::k_other, kRed }
        };
    }

    std::map<EventCategory, std::string> event_category_to_label_map_;
    std::map<EventCategory, int> event_category_to_color_map_;
};

#endif