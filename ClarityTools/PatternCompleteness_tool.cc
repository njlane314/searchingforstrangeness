#ifndef PATTERNCOMPLETENESS_H
#define PATTERNCOMPLETENESS_H

#include "ClarityToolBase.h"
#include <unordered_map>

namespace signature 
{
    class PatternCompleteness : public ClarityToolBase {
    public:
        explicit PatternCompleteness(const fhicl::ParameterSet& pset)
            : ClarityToolBase(pset),
            _pattHitCompThresh{pset.get<double>("PattHitCompThresh")},
            _pattHitThresh{pset.get<double>("PattHitThresh")},
            _sigHitCompThresh{pset.get<double>("SigHitCompThresh")} {}

        bool filter(const art::Event& e, const Signature& sig,
                    SignatureType /*type*/, common::PandoraView view) const override {
            if (!loadEventHandles(e, view)) {
                _metrics.hit_completeness = 0.0;
                _metrics.total_hits = 0.0;
                return false;
            }

            std::vector<art::Ptr<recob::Hit>> sigHits;
            _mcpHitMap.clear();
            _metrics.total_hits = 0;
            for (const auto& mcp_s : sig) {
                int sigHitCount = 0;
                for (const auto& hit : _mc_hits) {
                    const auto& assParticles = _mcp_bkth_assoc->at(hit.key());
                    const auto& assData = _mcp_bkth_assoc->data(hit.key());
                    for (size_t i = 0; i < assParticles.size(); ++i) {
                        if (assParticles[i]->TrackId() == mcp_s->TrackId() && assData[i]->isMaxIDEN == 1) {
                            sigHits.push_back(hit);
                            sigHitCount += 1;
                        }
                    }
                }
                _mcpHitMap[mcp_s->TrackId()] = sigHitCount;
                _metrics.total_hits += sigHitCount;
            }
            _metrics.hit_completeness = _mc_hits.empty() ? 0.0 : static_cast<double>(sigHits.size()) / _mc_hits.size();
            if (_mc_hits.empty() || _metrics.hit_completeness < _pattHitCompThresh || _metrics.total_hits < _pattHitThresh) {
                return false;
            }
            for (const auto& [trackId, numHits] : _mcpHitMap) {
                if (static_cast<double>(numHits) / _metrics.total_hits < _sigHitCompThresh) {
                    return false;
                }
            }
            return true;
        }

        std::string getToolName() const override {
            return "PatternCompleteness";
        }

        std::unique_ptr<ClarityMetrics> getMetrics() const override {
            return std::make_unique<CompletenessMetrics>(_metrics);
        }

    private:
        double _pattHitCompThresh;
        double _pattHitThresh;
        double _sigHitCompThresh;
        mutable std::vector<art::Ptr<recob::Hit>> _mc_hits; 
        mutable std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc; 
        mutable std::unordered_map<int, int> _mcpHitMap;
        mutable CompletenessMetrics _metrics; 

        bool loadEventHandles(const art::Event& e, common::PandoraView view) const {
            art::Handle<std::vector<recob::Hit>> hitHandle;
            if (!e.getByLabel(_hitProducer, hitHandle) || !hitHandle.isValid()) {
                return false;
            }

            _mc_hits.clear();
            for (size_t i = 0; i < hitHandle->size(); ++i) {
                if ((*hitHandle)[i].View() == static_cast<int>(view)) {
                    _mc_hits.emplace_back(hitHandle, i);
                }
            }

            _mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
                hitHandle, e, _backtrackTag);
            return _mcp_bkth_assoc->isValid();
        }
    };

    DEFINE_ART_CLASS_TOOL(PatternCompleteness)
} 

#endif