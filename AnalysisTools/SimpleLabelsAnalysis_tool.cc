#ifndef ANALYSIS_SIMPLELABELS_CXX
#define ANALYSIS_SIMPLELABELS_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <vector>
#include <map>
#include <string>
#include <array>

namespace analysis
{
    enum class Label {
        MIP,
        HIP,
        shower,
        michel,
        diffuse,
        invisible,
        undefined
    };

    class SimpleLabelsAnalyser : public AnalysisToolBase {
    public:
        explicit SimpleLabelsAnalyser(const fhicl::ParameterSet& pset);
        ~SimpleLabelsAnalyser() = default;

        SimpleLabelsAnalyser(const SimpleLabelsAnalyser&) = delete;
        SimpleLabelsAnalyser(SimpleLabelsAnalyser&&) = delete;
        SimpleLabelsAnalyser& operator=(const SimpleLabelsAnalyser&) = delete;
        SimpleLabelsAnalyser& operator=(SimpleLabelsAnalyser&&) = delete;

        void configure(const fhicl::ParameterSet& pset) override;
        void analyseEvent(const art::Event& event, bool is_data) override;
        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag _particle_label;
        double _gamma_threshold;
        double _hadron_threshold;
        std::vector<std::string> _particle_labels;
        std::map<int, size_t> _id_to_index;

        std::pair<Label, Label> compute_label(const std::vector<simb::MCParticle>& particles, const simb::MCParticle& part, Label sl_from_parent);
        void process_particle(size_t idx, const std::vector<simb::MCParticle>& particles, Label sl_from_parent);
    };

    const std::array<std::string, 7> label_names = {
        "MIP",
        "HIP",
        "shower",
        "michel",
        "diffuse",
        "invisible",
        "undefined"
    };

    std::string label_to_string(Label label) {
        size_t index = static_cast<size_t>(label);
        if (index < label_names.size()) {
            return label_names[index];
        }
        return "unknown";
    }

    SimpleLabelsAnalyser::SimpleLabelsAnalyser(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void SimpleLabelsAnalyser::configure(const fhicl::ParameterSet& pset) {
        _particle_label = pset.get<art::InputTag>("particle_label", "largeant");
        _gamma_threshold = pset.get<double>("gamma_threshold", 0.02);
        _hadron_threshold = pset.get<double>("hadron_threshold", 0.2);
    }

    void SimpleLabelsAnalyser::analyseEvent(const art::Event& event, bool is_data) {
        if (is_data) {
            _particle_labels.clear();
            return;
        }

        auto const& particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(_particle_label);
        const auto& particles = *particle_handle;

        _id_to_index.clear();
        for (size_t i = 0; i < particles.size(); ++i) {
            _id_to_index[particles[i].TrackId()] = i;
        }

        _particle_labels.resize(particles.size());

        for (size_t i = 0; i < particles.size(); ++i) {
            if (particles[i].Mother() == 0) {
                process_particle(i, particles, Label::undefined);
            }
        }
    }

    void SimpleLabelsAnalyser::process_particle(size_t idx, const std::vector<simb::MCParticle>& particles, Label sl_from_parent) {
        const auto& part = particles[idx];
        auto [sl, slc] = compute_label(particles, part, sl_from_parent);
        _particle_labels[idx] = label_to_string(sl);
        for (int i = 0; i < part.NumberDaughters(); ++i) {
            int daughter_id = part.Daughter(i);
            auto it = _id_to_index.find(daughter_id);
            if (it != _id_to_index.end()) {
                process_particle(it->second, particles, slc);
            }
        }
    }

    std::pair<Label, Label> SimpleLabelsAnalyser::compute_label(const std::vector<simb::MCParticle>& particles, const simb::MCParticle& part, Label sl_from_parent) {
        if (sl_from_parent != Label::undefined) {
            return {sl_from_parent, sl_from_parent};
        }

        Label sl = Label::invisible;
        Label slc = Label::undefined;

        int pdg = part.PdgCode();
        double momentum = part.P();
        std::string start_process = part.Process();
        std::string end_process = part.EndProcess();
        int parent_id = part.Mother();
        int parent_pdg = 0;
        if (parent_id != 0) {
            auto it = _id_to_index.find(parent_id);
            if (it != _id_to_index.end()) {
                parent_pdg = particles[it->second].PdgCode();
            }
        }

        if (pdg == 211 || pdg == -211 || pdg == 13 || pdg == -13) {
            sl = Label::MIP;
        }
        else if (pdg == 321 || pdg == -321 || (std::abs(pdg) == 2212 && momentum >= _hadron_threshold)) {
            sl = Label::HIP;
        }
        else if (pdg == 11 || pdg == -11) {
            if (start_process == "primary") {
                sl = Label::shower;
                slc = Label::shower;
            }
            else if (std::abs(parent_pdg) == 13 && (start_process == "muMinusCaptureAtRest" || start_process == "muPlusCaptureAtRest" || start_process == "Decay")) {
                sl = Label::michel;
                slc = Label::michel;
            }
            else if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                if (momentum >= _gamma_threshold) {
                    sl = Label::shower;
                    slc = Label::shower;
                }
                else {
                    sl = Label::diffuse;
                }
            }
            else {
                sl = Label::diffuse;
            }
        }
        else if (pdg == 22) {
            if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                if (momentum >= _gamma_threshold) {
                    sl = Label::shower;
                    slc = Label::shower;
                }
                else {
                    sl = Label::diffuse;
                }
            }
            else {
                sl = Label::diffuse;
            }
        }
        else if (std::abs(pdg) == 2212 && momentum < _hadron_threshold) {
            sl = Label::diffuse;
        }

        return {sl, slc};
    }

    void SimpleLabelsAnalyser::analyseSlice(const art::Event&, std::vector<common::ProxyPfpElem_t>&, bool, bool) {}

    void SimpleLabelsAnalyser::setBranches(TTree* tree) {
        tree->Branch("particle_labels", &_particle_labels);
    }

    void SimpleLabelsAnalyser::resetTTree(TTree*) {
        _particle_labels.clear();
    }

    DEFINE_ART_CLASS_TOOL(SimpleLabelsAnalyser)
}

#endif