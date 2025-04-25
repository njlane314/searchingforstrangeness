#pragma once

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

enum class NtupleFileType {
    kOnBNB,
    kExtBNB,
    kNumuMC,
    kIntrinsicNueMC,
    kDirtMC,
    kDetVarMCCV,
    kDetVarMCLYatten,
    kDetVarMCLYdown,
    kDetVarMCLYrayl,
    kDetVarMCRecomb2,
    kDetVarMCSCE,
    kDetVarMCWMAngleXZ,
    kDetVarMCWMAngleYZ,
    kDetVarMCWMdEdx,
    kDetVarMCWMX,
    kDetVarMCWMYZ,
    kDetVarMCCVExtra,
    kAltCVMC,
    kUnknown,
};

bool ntuple_type_is_detVar(const NtupleFileType& type) {
    constexpr std::array<NtupleFileType, 12> detVar_types = {
        NtupleFileType::kDetVarMCCV, NtupleFileType::kDetVarMCLYatten,
        NtupleFileType::kDetVarMCLYdown, NtupleFileType::kDetVarMCLYrayl,
        NtupleFileType::kDetVarMCRecomb2, NtupleFileType::kDetVarMCSCE,
        NtupleFileType::kDetVarMCWMAngleXZ, NtupleFileType::kDetVarMCWMAngleYZ,
        NtupleFileType::kDetVarMCWMdEdx, NtupleFileType::kDetVarMCWMX,
        NtupleFileType::kDetVarMCWMYZ, NtupleFileType::kDetVarMCCVExtra
    };
    const auto begin = detVar_types.cbegin();
    const auto end = detVar_types.cend();
    const auto iter = std::find(begin, end, type);
    if (iter != end) return true;
    return false;
}

bool ntuple_type_is_altCV(const NtupleFileType& type) {
    if (type == NtupleFileType::kAltCVMC) return true;
    return false;
}

bool ntuple_type_is_mc(const NtupleFileType& type) {
    if (type != NtupleFileType::kOnBNB && type != NtupleFileType::kExtBNB) return true;
    return false;
}

bool ntuple_type_is_reweightable_mc(const NtupleFileType& type) {
    if (type == NtupleFileType::kNumuMC || type == NtupleFileType::kIntrinsicNueMC || type == NtupleFileType::kDirtMC) return true;
    return false;
}

class FilePropertiesManager {
public:
    FilePropertiesManager(const FilePropertiesManager&) = delete;
    FilePropertiesManager(FilePropertiesManager&&) = delete;
    FilePropertiesManager& operator=(const FilePropertiesManager&) = delete;
    FilePropertiesManager& operator=(FilePropertiesManager&&) = delete;

    static FilePropertiesManager& Instance() {
        static std::unique_ptr<FilePropertiesManager> the_instance(new FilePropertiesManager());
        return *the_instance;
    }

    struct TriggersAndPOT {
        TriggersAndPOT() {}
        TriggersAndPOT(int trig, double pot) : trigger_count_(trig), pot_(pot) {}
        int trigger_count_ = 0;
        double pot_ = 0.;
    };

    const std::map<std::string, TriggersAndPOT>& data_norm_map() const {
        return data_norm_map_;
    }

    const std::map<int, std::map<NtupleFileType, std::set<std::string>>>& ntuple_file_map() const {
        return ntuple_file_map_;
    }

    NtupleFileType get_ntuple_file_type(const std::string& file_name) const {
        for (const auto& run_pair : ntuple_file_map_) {
            const auto& type_map = run_pair.second;
            for (const auto& type_pair : type_map) {
                const auto& file_set = type_pair.second;
                for (const auto& name : file_set) {
                    if (file_name == name) return type_pair.first;
                }
            }
        }
        throw std::runtime_error("ntuple file not found");
        return NtupleFileType::kOnBNB;
    }

    std::string ntuple_type_to_string(const NtupleFileType& type) const {
        std::string result;
        for (const auto& pair : string_to_file_type_map_) {
            if (pair.second == type) {
                result = pair.first;
                break;
            }
        }
        return result;
    }

    NtupleFileType string_to_ntuple_type(const std::string& str) const {
        auto end = string_to_file_type_map_.cend();
        auto iter = string_to_file_type_map_.find(str);
        if (iter != end) return iter->second;
        return NtupleFileType::kUnknown;
    }

    const std::string& analysis_path() const {
        return analysis_path_;
    }

    void load_file_properties(const std::string& input_table_file_name = "") {
        ntuple_file_map_.clear();
        data_norm_map_.clear();

        const char* path = std::getenv("STV_ANALYSIS_DIR");
        if (path == nullptr) throw std::runtime_error("The environment variable STV_ANALYSIS_DIR is not set. Please set it and try again.");

        analysis_path_ = path;

        std::string in_file_name(input_table_file_name);
        if (in_file_name.empty()) {
            in_file_name = analysis_path_ + "/file_properties.txt";
        }

        config_file_name_ = in_file_name;

        std::ifstream in_file(in_file_name);
        if (!in_file) {
            throw std::runtime_error("The file properties configuration file \"" + in_file_name + "\" could not be opened.");
        }

        std::string temp_line;
        while (std::getline(in_file, temp_line)) {
            if (temp_line.front() == '#') continue;

            std::string file_name;
            int run;
            std::string type_str;
            std::istringstream temp_ss(temp_line);
            temp_ss >> file_name >> run >> type_str;

            NtupleFileType type = string_to_file_type_map_.at(type_str);

            if (!ntuple_file_map_.count(run)) {
                ntuple_file_map_[run] = std::map<NtupleFileType, std::set<std::string>>();
            }
            auto& run_map = ntuple_file_map_.at(run);

            if (!run_map.count(type)) {
                run_map[type] = std::set<std::string>();
            }
            auto& file_set = run_map.at(type);

            file_set.insert(file_name);

            if (type == NtupleFileType::kOnBNB || type == NtupleFileType::kExtBNB) {
                int trigs;
                double pot;
                temp_ss >> trigs >> pot;
                data_norm_map_[file_name] = TriggersAndPOT(trigs, pot);
            }
        }
    }

    const std::string& config_file_name() const {
        return config_file_name_;
    }

private:
    FilePropertiesManager() {
        this->load_file_properties();
    }

    std::map<int, std::map<NtupleFileType, std::set<std::string>>> ntuple_file_map_;
    std::map<std::string, TriggersAndPOT> data_norm_map_;
    std::map<std::string, NtupleFileType> string_to_file_type_map_ = {
        { "onBNB", NtupleFileType::kOnBNB },
        { "extBNB", NtupleFileType::kExtBNB },
        { "numuMC", NtupleFileType::kNumuMC },
        { "nueMC", NtupleFileType::kIntrinsicNueMC },
        { "dirtMC", NtupleFileType::kDirtMC },
        { "detVarCV", NtupleFileType::kDetVarMCCV },
        { "detVarLYatten", NtupleFileType::kDetVarMCLYatten },
        { "detVarLYdown", NtupleFileType::kDetVarMCLYdown },
        { "detVarLYrayl", NtupleFileType::kDetVarMCLYrayl },
        { "detVarRecomb2", NtupleFileType::kDetVarMCRecomb2 },
        { "detVarSCE", NtupleFileType::kDetVarMCSCE },
        { "detVarWMAngleXZ", NtupleFileType::kDetVarMCWMAngleXZ },
        { "detVarWMAngleYZ", NtupleFileType::kDetVarMCWMAngleYZ },
        { "detVarWMdEdx", NtupleFileType::kDetVarMCWMdEdx },
        { "detVarWMX", NtupleFileType::kDetVarMCWMX },
        { "detVarWMYZ", NtupleFileType::kDetVarMCWMYZ },
        { "detVarCVExtra", NtupleFileType::kDetVarMCCVExtra },
        { "altCVMC", NtupleFileType::kAltCVMC },
    };
    std::string analysis_path_;
    std::string config_file_name_;
};