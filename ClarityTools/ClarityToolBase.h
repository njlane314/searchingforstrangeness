#ifndef CLARITY_TOOLBASE_H
#define CLARITY_TOOLBASE_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "../CommonFunctions/Pandora.h"
#include "../SignatureTools/SignatureToolBase.h"

namespace signature 
{
    struct ClarityMetrics {
        virtual ~ClarityMetrics() = default;
        virtual std::unique_ptr<ClarityMetrics> clone() const = 0;
    };

    struct CompletenessMetrics : ClarityMetrics {
        double hit_completeness;
        double total_hits;
        CompletenessMetrics(double hc = 0.0, double th = 0.0) : hit_completeness(hc), total_hits(th) {}
        std::unique_ptr<ClarityMetrics> clone() const override { return std::make_unique<CompletenessMetrics>(*this); }
    };

    struct ExclusivityMetrics : ClarityMetrics {
        double exclusivity_ratio;
        double hit_exclusivity_fraction;
        ExclusivityMetrics(double er = 0.0, double hef = 0.0) 
            : exclusivity_ratio(er), hit_exclusivity_fraction(hef) {}
        std::unique_ptr<ClarityMetrics> clone() const override { 
            return std::make_unique<ExclusivityMetrics>(*this); 
        }
    };

    struct IntegrityMetrics : ClarityMetrics {
        std::vector<bool> start_active;
        std::vector<bool> end_active;
        IntegrityMetrics() = default;
        std::unique_ptr<ClarityMetrics> clone() const override { return std::make_unique<IntegrityMetrics>(*this); }
    };
        
    class ClarityToolBase {
    public:
        ClarityToolBase(fhicl::ParameterSet const& pset) :
        _HITproducer{pset.get<art::InputTag>("HitProducer", "gaushit")}, 
        _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")},
        _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}, 
        _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")},
        _DeadChannelTag{pset.get<art::InputTag>("DeadChannelTag", "nfbadchannels:badchannels:OverlayDetsim")}, 
        _geoService{art::ServiceHandle<geo::Geometry>()->provider()} {}

        virtual ~ClarityToolBase() noexcept = default;

        std::map<common::PandoraView, std::vector<bool>> filterThreePlane(
            const art::Event& e, const Pattern& pattern) const {
            std::map<common::PandoraView, std::vector<bool>> results;
            for (int v = common::TPC_VIEW_U; v <static_cast<int>(common::N_VIEWS); ++v) {
                auto view = static_cast<common::PandoraView>(v);
                results[view] = filter(e, pattern, view);
            }
            return results;
        }

        std::vector<bool> filter(const art::Event& e, const Pattern& pattern,
                                common::PandoraView view) const {
            std::vector<bool> results;
            results.reserve(pattern.size());
            for (const auto& [type, sig] : pattern) {
                results.push_back(filter(e, sig, type, view));
            }
            return results;
        }

        virtual bool filter(const art::Event& e, const Signature& sig,
                        SignatureType type, common::PandoraView view) const = 0;

        virtual std::string getToolName() const = 0;

        virtual std::unique_ptr<ClarityMetrics> getMetrics() const = 0;

    protected:
        const art::InputTag _HITproducer;
        const art::InputTag _MCPproducer;
        const art::InputTag _MCTproducer;
        const art::InputTag _BacktrackTag;
        const art::InputTag _DeadChannelTag;
        const geo::GeometryCore* _geoService;
    };

    class ClarityEvaluator {
    public:
        class Result {
        public:
            Result(SignatureType t, const Signature& s) : signature(s), type(t) {
                pass_flags.fill(true);
            }

            Result(const Result&) = delete;
            Result& operator=(const Result&) = delete;

            Result(Result&&) = default;
            Result& operator=(Result&&) = default;

            bool passes(common::PandoraView view) const {
                return pass_flags[static_cast<size_t>(view)];
            }
            bool allPass() const {
                return std::all_of(pass_flags.begin(), pass_flags.end(), [](bool p){ return p; });
            }
            void addMetrics(const std::string& tool_name, common::PandoraView view,
                            std::unique_ptr<ClarityMetrics> metrics) {
                tool_metrics[tool_name][view] = std::move(metrics);
            }
            Signature signature;
            SignatureType type;
            std::array<bool, common::N_VIEWS> pass_flags;
            std::map<std::string, std::map<common::PandoraView, std::unique_ptr<ClarityMetrics>>> tool_metrics;
        };

        explicit ClarityEvaluator(const std::vector<std::unique_ptr<ClarityToolBase>>& clarityTools)
            : _clarityToolsVec(clarityTools) {}

        void evaluate(const art::Event& e, SignatureType type, const Signature& sig) const {
            _results.clear();
            Result result(type, sig);
            for (int view = common::TPC_VIEW_U; view < static_cast<int>(common::N_VIEWS); ++view) {
                bool viewPass = true;
                auto pandoraView = static_cast<common::PandoraView>(view);
                for (const auto& tool : _clarityToolsVec) {
                    bool toolPass = tool->filter(e, sig, type, pandoraView);
                    if (!toolPass) {
                        viewPass = false;
                    }
                    result.addMetrics(tool->getToolName(), pandoraView, tool->getMetrics());
                }
                result.pass_flags[static_cast<size_t>(view)] = viewPass;
            }
            _results.push_back(std::move(result));
        }

        const std::vector<Result>& getResults() const {
            return _results;
        }

    private:
        const std::vector<std::unique_ptr<ClarityToolBase>>& _clarityToolsVec;
        mutable std::vector<Result> _results;
    };
} 

#endif