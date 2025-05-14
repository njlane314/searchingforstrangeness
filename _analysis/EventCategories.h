#pragma once

#include <string>
#include <vector>
#include <map>
#include <stdexcept> // For std::out_of_range

// ROOT color constants can be included if preferred, e.g. #include "TColor.h"
// For simplicity, we can use integer color codes that ROOT understands.

namespace AnalysisFramework {

struct CategoryDisplayInfo {
    std::string label; // Human-readable label, can include LaTeX for plotting
    int color;         // ROOT color index (e.g., kRed, kBlue) or a custom hex
    std::string short_label; // Optional shorter label

    CategoryDisplayInfo(std::string l = "Unknown", int c = 1, std::string sl = "")
        : label(std::move(l)), color(c), short_label(std::move(sl)) {
        if (short_label.empty()) short_label = label;
    }
};

class EventCategories {
public:
    EventCategories() {
        // Populate categories based on the scheme from AddEventCategories
        // Data, EXT, Dirt
        categories_[0]   = {"On-Beam Data", kBlack, "Data"};
        categories_[1]   = {"Beam-Off (EXT)", kGray + 2, "EXT"};
        categories_[2]   = {"Dirt MC", kOrange + 7, "Dirt"};

        // Strangeness Signal Categories (NuMu CC focus)
        categories_[10]  = {R"($\nu_{\mu}$ CC Single Strange)", kRed + 2, R"($\nu_{\mu}$CC 1S)"};
        categories_[11]  = {R"($\nu_{\mu}$ CC Multiple Strange)", kRed - 3, R"($\nu_{\mu}$CC MS)"};
        categories_[12]  = {R"($\nu_e$ CC Single Strange)", kMagenta + 2, R"($\nu_e$CC 1S)"};
        categories_[13]  = {R"($\nu_e$ CC Multiple Strange)", kMagenta - 3, R"($\nu_e$CC MS)"};
        categories_[14]  = {"NC Single Strange", kPink + 1, "NC 1S"};
        categories_[15]  = {"NC Multiple Strange", kPink + 6, "NC MS"};
        categories_[19]  = {"Other True Strange", kRed - 9, "Oth.Str."};
        
        // Non-Strange Backgrounds: NuMu CC
        categories_[100] = {R"($\nu_{\mu}$ CC 0$\pi^{\pm}$ 0p)", kBlue + 2, R"($\nu_{\mu}$CC0$\pi$0p)"};
        categories_[101] = {R"($\nu_{\mu}$ CC 0$\pi^{\pm}$ 1p)", kBlue, R"($\nu_{\mu}$CC0$\pi$1p)"};
        categories_[102] = {R"($\nu_{\mu}$ CC 0$\pi^{\pm}$ Np)", kBlue - 4, R"($\nu_{\mu}$CC0$\pi$Np)"};
        categories_[103] = {R"($\nu_{\mu}$ CC 1$\pi^{\pm}$ 0p)", kCyan + 2, R"($\nu_{\mu}$CC1$\pi$0p)"};
        categories_[104] = {R"($\nu_{\mu}$ CC 1$\pi^{\pm}$ 1p)", kCyan, R"($\nu_{\mu}$CC1$\pi$1p)"};
        categories_[105] = {R"($\nu_{\mu}$ CC 1$\pi^{\pm}$ Np)", kCyan - 3, R"($\nu_{\mu}$CC1$\pi$Np)"};
        categories_[106] = {R"($\nu_{\mu}$ CC M$\pi^{\pm}$ AnyP)", kTeal + 2, R"($\nu_{\mu}$CCM$\pi$)"};

        // Non-Strange Backgrounds: NuMu NC
        categories_[110] = {R"($\nu_{\mu}$ NC 0$\pi^{\pm}$ 0p)", kGreen + 3, R"($\nu_{\mu}$NC0$\pi$0p)"};
        categories_[111] = {R"($\nu_{\mu}$ NC 0$\pi^{\pm}$ 1p)", kGreen + 1, R"($\nu_{\mu}$NC0$\pi$1p)"};
        categories_[112] = {R"($\nu_{\mu}$ NC 0$\pi^{\pm}$ Np)", kGreen - 5, R"($\nu_{\mu}$NC0$\pi$Np)"};
        categories_[113] = {R"($\nu_{\mu}$ NC 1$\pi^{\pm}$ 0p)", kSpring + 9, R"($\nu_{\mu}$NC1$\pi$0p)"};
        categories_[114] = {R"($\nu_{\mu}$ NC 1$\pi^{\pm}$ 1p)", kSpring + 5, R"($\nu_{\mu}$NC1$\pi$1p)"};
        categories_[115] = {R"($\nu_{\mu}$ NC 1$\pi^{\pm}$ Np)", kSpring - 5, R"($\nu_{\mu}$NC1$\pi$Np)"};
        categories_[116] = {R"($\nu_{\mu}$ NC M$\pi^{\pm}$ AnyP)", kYellow + 2, R"($\nu_{\mu}$NCM$\pi$)"};

        // Non-Strange Backgrounds: NuE CC
        categories_[200] = {R"($\nu_e$ CC 0$\pi^{\pm}$ 0p)", kOrange + 1, R"($\nu_e$CC0$\pi$0p)"};
        categories_[201] = {R"($\nu_e$ CC 0$\pi^{\pm}$ 1p)", kOrange - 3, R"($\nu_e$CC0$\pi$1p)"};
        categories_[202] = {R"($\nu_e$ CC 0$\pi^{\pm}$ Np)", kOrange - 9, R"($\nu_e$CC0$\pi$Np)"};
        categories_[203] = {R"($\nu_e$ CC 1$\pi^{\pm}$ 0p)", kOrange + 2, R"($\nu_e$CC1$\pi$0p)"};
        categories_[204] = {R"($\nu_e$ CC 1$\pi^{\pm}$ 1p)", kOrange - 2, R"($\nu_e$CC1$\pi$1p)"};
        categories_[205] = {R"($\nu_e$ CC 1$\pi^{\pm}$ Np)", kOrange - 8, R"($\nu_e$CC1$\pi$Np)"};
        categories_[206] = {R"($\nu_e$ CC M$\pi^{\pm}$ AnyP)", kOrange + 7, R"($\nu_e$CCM$\pi$)"};
        
        // Non-Strange Backgrounds: NuE NC
        categories_[210] = {R"($\nu_e$ NC Non-Strange)", kYellow - 3, R"($\nu_e$NC)"};
        
        // Fallback/Other
        categories_[900]  = {"Non-Strange in IS MC", kViolet - 5, "IS NonS"};
        categories_[998]  = {"Other MC (non $\nu_e$, $\nu_{\mu}$)", kGray, "Oth.MC"};
        categories_[9999] = {"Uncategorized", kBlack, "Uncat."};
    }

    const CategoryDisplayInfo& GetCategoryInfo(int category_id) const {
        try {
            return categories_.at(category_id);
        } catch (const std::out_of_range& oor) {
            // Return a default "Unknown" category if ID not found
            static const CategoryDisplayInfo unknown_info("Unknown Category", kGray + 1, "Unknown");
            return unknown_info;
        }
    }

    std::string GetLabel(int category_id) const {
        return GetCategoryInfo(category_id).label;
    }

    std::string GetShortLabel(int category_id) const {
        return GetCategoryInfo(category_id).short_label;
    }

    int GetColor(int category_id) const {
        return GetCategoryInfo(category_id).color;
    }

    std::vector<int> GetAllCategoryIds() const {
        std::vector<int> ids;
        for(const auto& pair : categories_) {
            ids.push_back(pair.first);
        }
        std::sort(ids.begin(), ids.end());
        return ids;
    }

private:
    std::map<int, CategoryDisplayInfo> categories_;
    // ROOT Color constants (examples, use actual values or TColor::GetColor(...) )
    // These are just standard ROOT color indices.
    static const int kBlack = 1;
    static const int kRed = 632; // kRed
    static const int kGreen = 416; // kGreen
    static const int kBlue = 600; // kBlue
    static const int kYellow = 400; // kYellow
    static const int kMagenta = 616; // kMagenta
    static const int kCyan = 432; // kCyan
    static const int kOrange = 800; // kOrange
    static const int kSpring = 820; // kSpring
    static const int kTeal = 840; // kTeal
    static const int kViolet = 880; // kViolet
    static const int kPink = 900; // kPink
    static const int kGray = 920; // kGray
};

} // namespace AnalysisFramework
