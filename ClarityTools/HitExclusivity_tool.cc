#ifndef HIT_EXCLUSIVITY_H
#define HIT_EXCLUSIVITY_H

#include "ClarityToolBase.h" 
#include <unordered_set>

namespace signature {
    
class HitExclusivity : public ClarityToolBase {
public:
    explicit HitExclusivity(const fhicl::ParameterSet& pset) :
        ClarityToolBase{pset}, 
        _hitExclusivityThreshold{pset.get<double>("HitExclusivityThreshold", 0.50)}, 
        _sigExclusivityThreshold{pset.get<double>("SignatureExclusivityThreshold", 0.70)} {} 

    bool filter(const art::Event& e, const Signature& sig,
                SignatureType /*type*/, common::PandoraView view) const override { 
        art::Handle<std::vector<recob::Hit>> hitHandle;
        e.getByLabel(_hitProducer, hitHandle);
        if (!hitHandle.isValid()) return false;

        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcpAssoc(
            hitHandle, e, _backtrackTag);
        if (!mcpAssoc.isValid()) return false;

        std::unordered_set<int> sigTrackIds;
        for (const auto& particle : sig) {
            sigTrackIds.insert(particle->TrackId());
        }

        double exclusiveCharge = 0.0;
        double inclusiveCharge = 0.0;
        std::unordered_set<size_t> associatedHits;
        std::unordered_set<size_t> exclusiveHits;

        for (size_t i = 0; i < hitHandle->size(); ++i) {
            const recob::Hit& hit = (*hitHandle)[i];
            if (hit.View() != static_cast<int>(view)) continue;

            size_t hitKey = i;
            const auto& assParticles = mcpAssoc.at(hitKey);
            const auto& assData = mcpAssoc.data(hitKey); 
            bool isAssociated = false;
            bool isExclusive = false;

            for (size_t j = 0; j < assParticles.size(); ++j) {
                if (sigTrackIds.count(assParticles[j]->TrackId())) {
                    isAssociated = true;
                    double charge = assData[j]->numElectrons * assData[j]->ideNFraction;
                    inclusiveCharge += charge;
                    if (assData[j]->ideNFraction > _hitExclusivityThreshold) {
                        exclusiveCharge += charge;
                        isExclusive = true;
                        break;
                    }
                }
            }
            if (isAssociated) {
                associatedHits.insert(hitKey);
                if (isExclusive) {
                    exclusiveHits.insert(hitKey);
                }
            }
        }

        _metrics.exclusivity_ratio = inclusiveCharge > 0 ? exclusiveCharge / inclusiveCharge : 0.0;
        _metrics.hit_exclusivity_fraction = associatedHits.empty() ? 0.0 :
            static_cast<double>(exclusiveHits.size()) / associatedHits.size();

        return _metrics.exclusivity_ratio > _sigExclusivityThreshold;
    }

    std::string getToolName() const override {
        return "HitExclusivity";
    }

    std::unique_ptr<ClarityMetrics> getMetrics() const override {
        return std::make_unique<ExclusivityMetrics>(_metrics);
    }

private:
    const double _hitExclusivityThreshold;
    const double _sigExclusivityThreshold;
    mutable ExclusivityMetrics _metrics;
};

DEFINE_ART_CLASS_TOOL(HitExclusivity)

}

#endif