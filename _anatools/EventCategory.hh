#pragma once

#include "TH1.h"
#include <map>
#include <memory>

enum class EventCategory {
    k_unknown = 0,
    k_nc = 1,
    k_nu_e_cc = 2,
    k_nu_mu_cc_with_strange = 3,
    k_nu_mu_cc_with_protons = 4,
    k_nu_mu_cc_with_pions = 5,
    k_nu_mu_cc_with_protons_pions = 6,
    k_nu_mu_cc_other = 7,
    k_other = 8
};

class EventCategoryInterpreter {
public:
    EventCategoryInterpreter(const EventCategoryInterpreter&) = delete;
    EventCategoryInterpreter(EventCategoryInterpreter&&) = delete;
    EventCategoryInterpreter& operator=(const EventCategoryInterpreter&) = delete;
    EventCategoryInterpreter& operator=(EventCategoryInterpreter&&) = delete;

    inline static const EventCategoryInterpreter& Instance() {
        static std::unique_ptr<EventCategoryInterpreter> the_instance(new EventCategoryInterpreter());
        return *the_instance;
    }

    inline const std::map<EventCategory, std::string>& label_map() const {
        return event_category_to_label_map_;
    }

    inline std::string label(EventCategory ec) const {
        return event_category_to_label_map_.at(ec);
    }

    inline int color_code(EventCategory ec) const {
        return event_category_to_color_map_.at(ec);
    }

    inline void set_histogram_style(TH1* hist, EventCategory ec) const {
        if (ec == EventCategory::k_unknown) {
            hist->SetLineColor(kBlack);
            hist->SetLineWidth(3);
            hist->SetMarkerStyle(kFullCircle);
            hist->SetMarkerSize(0.8);
            hist->SetStats(false);
            hist->GetXaxis()->SetTitleOffset(0.0);
            hist->GetXaxis()->SetTitleSize(0.0);
            hist->GetYaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->CenterTitle(true);
            hist->GetXaxis()->SetLabelSize(0.0);
            hist->SetMinimum(1e-3);
        } else {
            int color = color_code(ec);
            hist->SetFillColor(color);
            hist->SetLineColor(color);
            hist->SetStats(false);
        }
    }

private:
    EventCategoryInterpreter() {}

    std::map<EventCategory, std::string> event_category_to_label_map_ = {
        { EventCategory::k_unknown, "Unknown" },
        { EventCategory::k_nc, "NC" },
        { EventCategory::k_nu_e_cc, "#nu_{e} CC" },
        { EventCategory::k_nu_mu_cc_with_strange, "#nu_{#mu} CC with strange hadrons" },
        { EventCategory::k_nu_mu_cc_with_protons, "#nu_{#mu} CC with protons" },
        { EventCategory::k_nu_mu_cc_with_pions, "#nu_{#mu} CC with pions" },
        { EventCategory::k_nu_mu_cc_with_protons_pions, "#nu_{#mu} CC with protons and pions" },
        { EventCategory::k_nu_mu_cc_other, "Other #nu_{#mu} CC" },
        { EventCategory::k_other, "Other" }
    };

    std::map<EventCategory, int> event_category_to_color_map_ = {
        { EventCategory::k_unknown, kBlack },
        { EventCategory::k_nc, kOrange },
        { EventCategory::k_nu_e_cc, kViolet },
        { EventCategory::k_nu_mu_cc_with_strange, kAzure - 2 },
        { EventCategory::k_nu_mu_cc_with_protons, kAzure - 1 },
        { EventCategory::k_nu_mu_cc_with_pions, kAzure },
        { EventCategory::k_nu_mu_cc_with_protons_pions, kAzure + 1 },
        { EventCategory::k_nu_mu_cc_other, kAzure + 2 },
        { EventCategory::k_other, kRed + 1 }
    };
};