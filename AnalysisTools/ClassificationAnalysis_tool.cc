#ifndef CLASSIFICATION_ANALYSIS_H
#define CLASSIFICATION_ANALYSIS_H

#include "AnalysisToolBase.h"
#include "EventClassifier.h"
#include "art/Framework/Principal/Event.h"
#include "SignatureTools/SignatureToolBase.h"
#include "TTree.h"
#include <memory>
#include <vector>
#include <array>

namespace analysis 
{
    class ClassificationAnalysis : public AnalysisToolBase {
    public:
        explicit ClassificationAnalysis(const fhicl::ParameterSet& pset);
        ~ClassificationAnalysis() = default;

        ClassificationAnalysis(const ClassificationAnalysis&) = delete;
        ClassificationAnalysis(ClassificationAnalysis&&) = delete;
        ClassificationAnalysis& operator=(const ClassificationAnalysis&) = delete;
        ClassificationAnalysis& operator=(ClassificationAnalysis&&) = delete;

        void configure(const fhicl::ParameterSet& pset) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;
        void analyseEvent(const art::Event& event, bool is_data) override;
        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override {}

    private:
        struct ViewMetrics {
            std::vector<double> completeness_hit;
            std::vector<double> completeness_total_hits;
            std::vector<double> exclusivity_ratio;
            std::vector<double> hit_exclusivity_fraction;
            std::vector<int> particle_start_active;
            std::vector<int> particle_end_active;
        };

        struct MetricsVisitor {
            std::array<ViewMetrics, 3>& metrics;
            const signature::Signature& sig;

            MetricsVisitor(std::array<ViewMetrics, 3>& m, const signature::Signature& s) 
                : metrics(m), sig(s) {}

            void visit(const std::map<std::string, std::map<common::PandoraView, std::unique_ptr<signature::ClarityMetrics>>>& tool_metrics) {
                for (int view = common::TPC_VIEW_U; view <= common::TPC_VIEW_W; ++view) {
                    int idx = view - common::TPC_VIEW_U;
                    auto pandora_view = static_cast<common::PandoraView>(view);

                    if (auto tool_it = tool_metrics.find("PatternCompleteness"); tool_it != tool_metrics.end()) {
                        if (auto view_it = tool_it->second.find(pandora_view); view_it != tool_it->second.end()) {
                            if (auto cm = dynamic_cast<signature::CompletenessMetrics*>(view_it->second.get())) {
                                metrics[idx].completeness_hit.push_back(cm->hit_completeness);
                                metrics[idx].completeness_total_hits.push_back(cm->total_hits);
                            }
                        }
                    }

                    if (auto tool_it = tool_metrics.find("HitExclusivity"); tool_it != tool_metrics.end()) {
                        if (auto view_it = tool_it->second.find(pandora_view); view_it != tool_it->second.end()) {
                            if (auto em = dynamic_cast<signature::ExclusivityMetrics*>(view_it->second.get())) {
                                metrics[idx].exclusivity_ratio.push_back(em->exclusivity_ratio);
                                metrics[idx].hit_exclusivity_fraction.push_back(em->hit_exclusivity_fraction);
                            }
                        }
                    }

                    if (auto tool_it = tool_metrics.find("SignatureIntegrity"); tool_it != tool_metrics.end()) {
                        if (auto view_it = tool_it->second.find(pandora_view); view_it != tool_it->second.end()) {
                            if (auto im = dynamic_cast<signature::IntegrityMetrics*>(view_it->second.get())) {
                                for (size_t p = 0; p < sig.size(); ++p) {
                                    metrics[idx].particle_start_active.push_back(im->start_active[p] ? 1 : 0);
                                    metrics[idx].particle_end_active.push_back(im->end_active[p] ? 1 : 0);
                                }
                            }
                        }
                    }
                }
            }
        };

        std::unique_ptr<signature::EventClassifier> _classifier;
        int event_type_;
        int is_signal_;
        std::vector<int> signature_types_;
        std::vector<int> clarity_U_, clarity_V_, clarity_W_;
        int pass_clarity_;
        std::vector<int> particle_pdg_;
        std::vector<double> particle_energy_;
        std::vector<double> particle_px_, particle_py_, particle_pz_;
        std::array<ViewMetrics, 3> view_metrics_;
    };

    inline ClassificationAnalysis::ClassificationAnalysis(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    inline void ClassificationAnalysis::configure(const fhicl::ParameterSet& pset) {
        _classifier = std::make_unique<signature::EventClassifier>(pset);
    }

    inline void ClassificationAnalysis::setBranches(TTree* tree) {
        tree->Branch("event_type", &event_type_);
        tree->Branch("is_signal", &is_signal_);
        tree->Branch("signature_types", &signature_types_);
        tree->Branch("clarity_U", &clarity_U_);
        tree->Branch("clarity_V", &clarity_V_);
        tree->Branch("clarity_W", &clarity_W_);
        tree->Branch("pass_clarity", &pass_clarity_);
        tree->Branch("completeness_hit_U", &view_metrics_[0].completeness_hit);
        tree->Branch("completeness_hit_V", &view_metrics_[1].completeness_hit);
        tree->Branch("completeness_hit_W", &view_metrics_[2].completeness_hit);
        tree->Branch("completeness_total_hits_U", &view_metrics_[0].completeness_total_hits);
        tree->Branch("completeness_total_hits_V", &view_metrics_[1].completeness_total_hits);
        tree->Branch("completeness_total_hits_W", &view_metrics_[2].completeness_total_hits);
        tree->Branch("exclusivity_ratio_U", &view_metrics_[0].exclusivity_ratio);
        tree->Branch("exclusivity_ratio_V", &view_metrics_[1].exclusivity_ratio);
        tree->Branch("exclusivity_ratio_W", &view_metrics_[2].exclusivity_ratio);
        tree->Branch("hit_exclusivity_fraction_U", &view_metrics_[0].hit_exclusivity_fraction);
        tree->Branch("hit_exclusivity_fraction_V", &view_metrics_[1].hit_exclusivity_fraction);
        tree->Branch("hit_exclusivity_fraction_W", &view_metrics_[2].hit_exclusivity_fraction);
        tree->Branch("particle_pdg", &particle_pdg_);
        tree->Branch("particle_energy", &particle_energy_);
        tree->Branch("particle_px", &particle_px_);
        tree->Branch("particle_py", &particle_py_);
        tree->Branch("particle_pz", &particle_pz_);
        tree->Branch("particle_start_active_U", &view_metrics_[0].particle_start_active);
        tree->Branch("particle_start_active_V", &view_metrics_[1].particle_start_active);
        tree->Branch("particle_start_active_W", &view_metrics_[2].particle_start_active);
        tree->Branch("particle_end_active_U", &view_metrics_[0].particle_end_active);
        tree->Branch("particle_end_active_V", &view_metrics_[1].particle_end_active);
        tree->Branch("particle_end_active_W", &view_metrics_[2].particle_end_active);
    }

    inline void ClassificationAnalysis::resetTTree(TTree* tree) {
        event_type_ = -1;
        is_signal_ = 0;
        signature_types_.clear();
        clarity_U_.clear();
        clarity_V_.clear();
        clarity_W_.clear();
        pass_clarity_ = 0;
        particle_pdg_.clear();
        particle_energy_.clear();
        particle_px_.clear();
        particle_py_.clear();
        particle_pz_.clear();
        for (auto& vm : view_metrics_) {
            vm.completeness_hit.clear();
            vm.completeness_total_hits.clear();
            vm.exclusivity_ratio.clear();
            vm.hit_exclusivity_fraction.clear();
            vm.particle_start_active.clear();
            vm.particle_end_active.clear();
        }
    }

    inline void ClassificationAnalysis::analyseEvent(const art::Event& e, bool is_data) {
        if (is_data) return;

        event_type_ = static_cast<int>(_classifier->classifyEvent(e));
        is_signal_ = _classifier->isSignal(e) ? 1 : 0;

        const auto& pattern = _classifier->getPattern(e);
        const auto& clarity_results = _classifier->getClarityResults();

        if (pattern.size() != clarity_results.size()) return;

        for (size_t i = 0; i < pattern.size(); ++i) {
            const auto& [type, sig] = pattern[i];
            const auto& result = clarity_results[i];

            signature_types_.push_back(static_cast<int>(type));
            clarity_U_.push_back(result.passes(common::TPC_VIEW_U) ? 1 : 0);
            clarity_V_.push_back(result.passes(common::TPC_VIEW_V) ? 1 : 0);
            clarity_W_.push_back(result.passes(common::TPC_VIEW_W) ? 1 : 0);

            MetricsVisitor visitor(view_metrics_, sig);
            visitor.visit(result.tool_metrics);

            for (const auto& particle : sig) {
                particle_pdg_.push_back(particle->PdgCode());
                particle_energy_.push_back(particle->E());
                particle_px_.push_back(particle->Px());
                particle_py_.push_back(particle->Py());
                particle_pz_.push_back(particle->Pz());
            }
        }

        pass_clarity_ = _classifier->passClarity() ? 1 : 0;
    }

    DEFINE_ART_CLASS_TOOL(ClassificationAnalysis)
}

#endif