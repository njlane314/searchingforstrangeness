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
              _patt_hit_comp_thresh{pset.get<double>("PattHitCompThresh", 0.50)},
              _patt_hit_thresh{pset.get<double>("PattHitThresh", 100)},
              _sig_hit_comp_thresh{pset.get<double>("SigHitCompThresh", 0.10)} {}

        bool filter(const art::Event& e, const Signature& sig,
                    SignatureType /*type*/, common::PandoraView view) const override {
            art::Handle<std::vector<recob::Hit>> hit_handle;
            if (!e.getByLabel(_HITproducer, hit_handle) || !hit_handle.isValid()) 
                return false;

            std::vector<art::Ptr<recob::Hit>> mc_hits;
            for (size_t i = 0; i < hit_handle->size(); ++i) {
                if ((*hit_handle)[i].View() == static_cast<int>(view)) {
                    mc_hits.emplace_back(hit_handle, i);
                }
            }

            auto mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
                hit_handle, e, _BacktrackTag);
            if (!mcp_bkth_assoc->isValid()) 
                return false;

            std::vector<art::Ptr<recob::Hit>> sig_hits;
            _mcp_hit_map.clear();
            _metrics.total_hits = 0;
            for (const auto& mcp_s : sig) {
                int sig_hit_count = 0;
                for (const auto& hit : mc_hits) {
                    const auto& ass_part = mcp_bkth_assoc->at(hit.key());
                    const auto& ass_data = mcp_bkth_assoc->data(hit.key());
                    for (size_t i = 0; i < ass_part.size(); ++i) {
                        if (ass_part[i]->TrackId() == mcp_s->TrackId() && ass_data[i]->isMaxIDEN == 1) {
                            sig_hits.push_back(hit);
                            sig_hit_count += 1;
                        }
                    }
                }
                _mcp_hit_map[mcp_s->TrackId()] = sig_hit_count;
                _metrics.total_hits += sig_hit_count;
            }

            _metrics.hit_completeness = mc_hits.empty() ? 0.0 : static_cast<double>(sig_hits.size()) / mc_hits.size();
            if (mc_hits.empty() || _metrics.hit_completeness < _patt_hit_comp_thresh || _metrics.total_hits < _patt_hit_thresh) 
                return false;

            for (const auto& [track_id, num_hits] : _mcp_hit_map) {
                if (static_cast<double>(num_hits) / _metrics.total_hits < _sig_hit_comp_thresh) {
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
        double _patt_hit_comp_thresh;
        double _patt_hit_thresh;
        double _sig_hit_comp_thresh;
        mutable std::unordered_map<int, int> _mcp_hit_map;
        mutable CompletenessMetrics _metrics;
    };

    DEFINE_ART_CLASS_TOOL(PatternCompleteness)
}

#endif