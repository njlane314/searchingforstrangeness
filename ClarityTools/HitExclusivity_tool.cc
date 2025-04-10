#ifndef HIT_EXCLUSIVITY_H
#define HIT_EXCLUSIVITY_H

#include "ClarityToolBase.h" 
#include <unordered_set>

namespace signature 
{    
    class HitExclusivity : public ClarityToolBase {
    public:
        explicit HitExclusivity(const fhicl::ParameterSet& pset) :
            ClarityToolBase{pset}, 
            _hit_exclusivity_thresh{pset.get<double>("HitExclusivityThreshold", 0.50)}, 
            _sig_exclusivity_thresh{pset.get<double>("SignatureExclusivityThreshold", 0.70)} {} 

        bool filter(const art::Event& e, const Signature& sig,
                    SignatureType /*type*/, common::PandoraView view) const override { 
            art::Handle<std::vector<recob::Hit>> hit_handle;
            e.getByLabel(_HITproducer, hit_handle);
            if (!hit_handle.isValid()) 
                return false;

            art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(
                hit_handle, e, _BacktrackTag);
            if (!mcp_bkth_assoc.isValid()) 
                return false;

            std::unordered_set<int> sig_tracks;
            for (const auto& particle : sig) {
                sig_tracks.insert(particle->TrackId());
            }

            double exclusive_q = 0.0;
            double inclusive_q = 0.0;
            std::unordered_set<size_t> ass_hits;
            std::unordered_set<size_t> exc_hits;

            for (size_t i = 0; i < hit_handle->size(); ++i) {
                const recob::Hit& hit = (*hit_handle)[i];
                if (hit.View() != static_cast<int>(view)) 
                    continue;

                size_t hit_key = i;
                const auto& ass_part = mcp_bkth_assoc.at(hit_key);
                const auto& ass_data = mcp_bkth_assoc.data(hit_key); 
                bool is_assoc = false;
                bool is_excus = false;

                for (size_t j = 0; j < ass_part.size(); ++j) {
                    if (sig_tracks.count(ass_part[j]->TrackId())) {
                        is_assoc = true;
                        double charge = ass_data[j]->numElectrons * ass_data[j]->ideNFraction;
                        inclusive_q += charge;
                        if (ass_data[j]->ideNFraction > _hit_exclusivity_thresh) {
                            exclusive_q += charge;
                            is_excus = true;
                            break;
                        }
                    }
                }
                if (is_assoc) {
                    ass_hits.insert(hit_key);
                    if (is_excus) {
                        exc_hits.insert(hit_key);
                    }
                }
            }

            _metrics.exclusivity_ratio = inclusive_q > 0 ? exclusive_q / inclusive_q : 0.0;
            _metrics.hit_exclusivity_fraction = ass_hits.empty() ? 0.0 :
                static_cast<double>(exc_hits.size()) / ass_hits.size();

            return _metrics.exclusivity_ratio > _sig_exclusivity_thresh;
        }

        std::string getToolName() const override {
            return "HitExclusivity";
        }

        std::unique_ptr<ClarityMetrics> getMetrics() const override {
            return std::make_unique<ExclusivityMetrics>(_metrics);
        }

    private:
        const double _hit_exclusivity_thresh;
        const double _sig_exclusivity_thresh;
        mutable ExclusivityMetrics _metrics;
    };

    DEFINE_ART_CLASS_TOOL(HitExclusivity)
}

#endif