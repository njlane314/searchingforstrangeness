#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <vector>
#include <string>
#include "Sample.h"

namespace analysis {
    class Config {
    public:
        using SampleConfig = std::map<SampleType, std::vector<std::string>>;

        static void AddConfig(const std::string& name, const SampleConfig& config) {
            configs_[name] = config;
        }

        static SampleConfig GetConfig(const std::string& name) {
            auto it = configs_.find(name);
            if (it != configs_.end()) {
                return it->second;
            }
            return {};
        }

        static std::vector<std::string> GetFiles(const std::string& configName, SampleType type) {
            auto config = GetConfig(configName);
            auto it = config.find(type);
            if (it != config.end()) {
                return it->second;
            }
            return {};
        }

    private:
        static inline std::map<std::string, SampleConfig> configs_ = {
            {"default", {
                //{SampleType::SIGNAL, {"/exp/uboone/data/users/nlane/analysis/prod_strange_resample_fhc_run2_fhc_reco2_reco2_topologytrainingfilter_1_new_analysis.root"}}
                {SampleType::SIGNAL, {"/exp/uboone/data/users/nlane/analysis/prod_strange_resample_fhc_run2_fhc_reco2_reco2_imagetrainingfilter_10_new_analysis.root"}},
                //{SampleType::GENIE, {"/exp/uboone/data/users/nlane/analysis/Run2_FHC_New_Flux_Nu_Overlay_Pandora_reprocess_pandora_reco2_reco2_selectionfilter_20_new_analysis.root"}}
            }}
        };
    };
}

#endif