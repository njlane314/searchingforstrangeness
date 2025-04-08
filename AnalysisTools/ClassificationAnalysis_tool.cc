#ifndef CLASSIFICATION_ANALYSIS_H
#define CLASSIFICATION_ANALYSIS_H

#include "AnalysisToolBase.h"
#include "EventClassifier.h"
#include "art/Framework/Principal/Event.h"
#include "SignatureTools/SignatureToolBase.h"
#include "TTree.h"
#include <memory>
#include <vector>
#include <stdexcept>

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
        struct MetricsVisitor {
            std::vector<double>& completeness_hit_U, &completeness_hit_V, &completeness_hit_W;
            std::vector<double>& completeness_total_hits_U, &completeness_total_hits_V, &completeness_total_hits_W;
            std::vector<double>& exclusivity_ratio_U, &exclusivity_ratio_V, &exclusivity_ratio_W;
            std::vector<double>& hit_exclusivity_fraction_U, &hit_exclusivity_fraction_V, &hit_exclusivity_fraction_W;
            std::vector<int>& particle_start_active_U, &particle_start_active_V, &particle_start_active_W;
            std::vector<int>& particle_end_active_U, &particle_end_active_V, &particle_end_active_W;
            const signature::Signature& sig;

            MetricsVisitor(const signature::Signature& s,
                            std::vector<double>& chu, std::vector<double>& chv, std::vector<double>& chw,
                            std::vector<double>& ctu, std::vector<double>& ctv, std::vector<double>& ctw,
                            std::vector<double>& eru, std::vector<double>& erv, std::vector<double>& erw,
                            std::vector<double>& hefu, std::vector<double>& hefv, std::vector<double>& hefw,
                            std::vector<int>& sau, std::vector<int>& sav, std::vector<int>& saw,
                            std::vector<int>& eau, std::vector<int>& eav, std::vector<int>& eaw)
                : sig(s),
                    completeness_hit_U(chu), completeness_hit_V(chv), completeness_hit_W(chw),
                    completeness_total_hits_U(ctu), completeness_total_hits_V(ctv), completeness_total_hits_W(ctw),
                    exclusivity_ratio_U(eru), exclusivity_ratio_V(erv), exclusivity_ratio_W(erw),
                    hit_exclusivity_fraction_U(hefu), hit_exclusivity_fraction_V(hefv), hit_exclusivity_fraction_W(hefw),
                    particle_start_active_U(sau), particle_start_active_V(sav), particle_start_active_W(saw),
                    particle_end_active_U(eau), particle_end_active_V(eav), particle_end_active_W(eaw) {}

            void visit(const std::map<std::string, std::map<common::PandoraView, std::unique_ptr<signature::ClarityMetrics>>>& tool_metrics) {
                for (int view = common::TPC_VIEW_U; view <= common::TPC_VIEW_W; ++view) {
                    auto pandora_view = static_cast<common::PandoraView>(view);

                    auto tool_it = tool_metrics.find("PatternCompleteness");
                    if (tool_it != tool_metrics.end()) {
                        auto view_it = tool_it->second.find(pandora_view);
                        if (view_it != tool_it->second.end()) {
                            if (auto cm = dynamic_cast<signature::CompletenessMetrics*>(view_it->second.get())) {
                                getCompletenessHitVector(view).push_back(cm->hit_completeness);
                                getCompletenessTotalHitsVector(view).push_back(cm->total_hits);
                            }
                        }
                    }

                    tool_it = tool_metrics.find("HitExclusivity");
                    if (tool_it != tool_metrics.end()) {
                        auto view_it = tool_it->second.find(pandora_view);
                        if (view_it != tool_it->second.end()) {
                            if (auto em = dynamic_cast<signature::ExclusivityMetrics*>(view_it->second.get())) {
                                getExclusivityRatioVector(view).push_back(em->exclusivity_ratio);
                                getHitExclusivityFractionVector(view).push_back(em->hit_exclusivity_fraction);
                            }
                        }
                    }

                    tool_it = tool_metrics.find("SignatureIntegrity");
                    if (tool_it != tool_metrics.end()) {
                        auto view_it = tool_it->second.find(pandora_view);
                        if (view_it != tool_it->second.end()) {
                            if (auto im = dynamic_cast<signature::IntegrityMetrics*>(view_it->second.get())) {
                                auto& start_active = getParticleStartActiveVector(view);
                                auto& end_active = getParticleEndActiveVector(view);
                                for (size_t p = 0; p < sig.size(); ++p) {
                                    start_active.push_back(im->start_active[p] ? 1 : 0);
                                    end_active.push_back(im->end_active[p] ? 1 : 0);
                                }
                            }
                        }
                    }
                }
            }

        private:
            std::vector<double>& getCompletenessHitVector(int view) {
                switch (view) {
                    case common::TPC_VIEW_U: return completeness_hit_U;
                    case common::TPC_VIEW_V: return completeness_hit_V;
                    case common::TPC_VIEW_W: return completeness_hit_W;
                    default: throw std::invalid_argument("Invalid view");
                }
            }

            std::vector<double>& getCompletenessTotalHitsVector(int view) {
                switch (view) {
                    case common::TPC_VIEW_U: return completeness_total_hits_U;
                    case common::TPC_VIEW_V: return completeness_total_hits_V;
                    case common::TPC_VIEW_W: return completeness_total_hits_W;
                    default: throw std::invalid_argument("Invalid view");
                }
            }

            std::vector<double>& getExclusivityRatioVector(int view) {
                switch (view) {
                    case common::TPC_VIEW_U: return exclusivity_ratio_U;
                    case common::TPC_VIEW_V: return exclusivity_ratio_V;
                    case common::TPC_VIEW_W: return exclusivity_ratio_W;
                    default: throw std::invalid_argument("Invalid view");
                }
            }

            std::vector<double>& getHitExclusivityFractionVector(int view) {
                switch (view) {
                    case common::TPC_VIEW_U: return hit_exclusivity_fraction_U;
                    case common::TPC_VIEW_V: return hit_exclusivity_fraction_V;
                    case common::TPC_VIEW_W: return hit_exclusivity_fraction_W;
                    default: throw std::invalid_argument("Invalid view");
                }
            }

            std::vector<int>& getParticleStartActiveVector(int view) {
                switch (view) {
                    case common::TPC_VIEW_U: return particle_start_active_U;
                    case common::TPC_VIEW_V: return particle_start_active_V;
                    case common::TPC_VIEW_W: return particle_start_active_W;
                    default: throw std::invalid_argument("Invalid view");
                }
            }

            std::vector<int>& getParticleEndActiveVector(int view) {
                switch (view) {
                    case common::TPC_VIEW_U: return particle_end_active_U;
                    case common::TPC_VIEW_V: return particle_end_active_V;
                    case common::TPC_VIEW_W: return particle_end_active_W;
                    default: throw std::invalid_argument("Invalid view");
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
        std::vector<double> completeness_hit_U_, completeness_hit_V_, completeness_hit_W_;
        std::vector<double> completeness_total_hits_U_, completeness_total_hits_V_, completeness_total_hits_W_;
        std::vector<double> exclusivity_ratio_U_, exclusivity_ratio_V_, exclusivity_ratio_W_;
        std::vector<double> hit_exclusivity_fraction_U_, hit_exclusivity_fraction_V_, hit_exclusivity_fraction_W_;
        std::vector<int> particle_start_active_U_, particle_start_active_V_, particle_start_active_W_;
        std::vector<int> particle_end_active_U_, particle_end_active_V_, particle_end_active_W_;
    };

    inline ClassificationAnalysis::ClassificationAnalysis(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    inline void ClassificationAnalysis::configure(const fhicl::ParameterSet& pset) {
        _classifier = std::make_unique<signature::EventClassifier>(pset);
    }

    inline void ClassificationAnalysis::setBranches(TTree* tree) {
        tree->Branch("event_type", &event_type_, "event_type/I");
        tree->Branch("is_signal", &is_signal_, "is_signal/I");
        tree->Branch("signature_types", &signature_types_);
        tree->Branch("clarity_U", &clarity_U_);
        tree->Branch("clarity_V", &clarity_V_);
        tree->Branch("clarity_W", &clarity_W_);
        tree->Branch("pass_clarity", &pass_clarity_, "pass_clarity/I");
        tree->Branch("completeness_hit_U", &completeness_hit_U_);
        tree->Branch("completeness_hit_V", &completeness_hit_V_);
        tree->Branch("completeness_hit_W", &completeness_hit_W_);
        tree->Branch("completeness_total_hits_U", &completeness_total_hits_U_);
        tree->Branch("completeness_total_hits_V", &completeness_total_hits_V_);
        tree->Branch("completeness_total_hits_W", &completeness_total_hits_W_);
        tree->Branch("exclusivity_ratio_U", &exclusivity_ratio_U_);
        tree->Branch("exclusivity_ratio_V", &exclusivity_ratio_V_);
        tree->Branch("exclusivity_ratio_W", &exclusivity_ratio_W_);
        tree->Branch("hit_exclusivity_fraction_U", &hit_exclusivity_fraction_U_);
        tree->Branch("hit_exclusivity_fraction_V", &hit_exclusivity_fraction_V_);
        tree->Branch("hit_exclusivity_fraction_W", &hit_exclusivity_fraction_W_);
        tree->Branch("particle_pdg", &particle_pdg_);
        tree->Branch("particle_energy", &particle_energy_);
        tree->Branch("particle_px", &particle_px_);
        tree->Branch("particle_py", &particle_py_);
        tree->Branch("particle_pz", &particle_pz_);
        tree->Branch("particle_start_active_U", &particle_start_active_U_);
        tree->Branch("particle_start_active_V", &particle_start_active_V_);
        tree->Branch("particle_start_active_W", &particle_start_active_W_);
        tree->Branch("particle_end_active_U", &particle_end_active_U_);
        tree->Branch("particle_end_active_V", &particle_end_active_V_);
        tree->Branch("particle_end_active_W", &particle_end_active_W_);
    }

    inline void ClassificationAnalysis::resetTTree(TTree* tree) {
        event_type_ = -1;
        is_signal_ = 0;
        signature_types_.clear();
        clarity_U_.clear();
        clarity_V_.clear();
        clarity_W_.clear();
        pass_clarity_ = 0;
        completeness_hit_U_.clear();
        completeness_hit_V_.clear();
        completeness_hit_W_.clear();
        completeness_total_hits_U_.clear();
        completeness_total_hits_V_.clear();
        completeness_total_hits_W_.clear();
        exclusivity_ratio_U_.clear();
        exclusivity_ratio_V_.clear();
        exclusivity_ratio_W_.clear();
        hit_exclusivity_fraction_U_.clear();
        hit_exclusivity_fraction_V_.clear();
        hit_exclusivity_fraction_W_.clear();
        particle_pdg_.clear();
        particle_energy_.clear();
        particle_px_.clear();
        particle_py_.clear();
        particle_pz_.clear();
        particle_start_active_U_.clear();
        particle_start_active_V_.clear();
        particle_start_active_W_.clear();
        particle_end_active_U_.clear();
        particle_end_active_V_.clear();
        particle_end_active_W_.clear();
    }

    inline void ClassificationAnalysis::analyseEvent(const art::Event& e, bool is_data) {
        if (is_data) return;

        event_type_ = static_cast<int>(_classifier->classifyEvent(e));
        is_signal_ = _classifier->isSignal(e) ? 1 : 0;

        const auto& pattern = _classifier->getPattern(e);
        const auto& clarity_results = _classifier->getClarityResults();

        if (pattern.size() != clarity_results.size()) {
            throw std::runtime_error("Mismatch between pattern and clarity results sizes.");
        }

        for (size_t i = 0; i < pattern.size(); ++i) {
            const auto& [type, sig] = pattern[i];
            const auto& result = clarity_results[i];

            signature_types_.push_back(static_cast<int>(type));
            clarity_U_.push_back(result.passes(common::TPC_VIEW_U) ? 1 : 0);
            clarity_V_.push_back(result.passes(common::TPC_VIEW_V) ? 1 : 0);
            clarity_W_.push_back(result.passes(common::TPC_VIEW_W) ? 1 : 0);

            MetricsVisitor visitor(sig, completeness_hit_U_, completeness_hit_V_, completeness_hit_W_,
                                    completeness_total_hits_U_, completeness_total_hits_V_, completeness_total_hits_W_,
                                    exclusivity_ratio_U_, exclusivity_ratio_V_, exclusivity_ratio_W_,
                                    hit_exclusivity_fraction_U_, hit_exclusivity_fraction_V_, hit_exclusivity_fraction_W_,
                                    particle_start_active_U_, particle_start_active_V_, particle_start_active_W_,
                                    particle_end_active_U_, particle_end_active_V_, particle_end_active_W_);
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