#ifndef HIT_EXCLUSIVITY_H
#define HIT_EXCLUSIVITY_H

#include "ClarityToolBase.h" 

namespace claritytools {

class HitExclusivity : ClarityToolBase {

public:
    explicit HitExclusivity(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)}
    , _hit_exclus_thresh{pset.get<double>("HitExclusivityThreshold", 0.5)}
    , _sig_exclus_thresh{pset.get<double>("SignatureExclusivityThreshold", 0.8)}
    {
        configure(pset);
    }

    ~HitExclusivity() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    //bool filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filter(const art::Event &e, const signature::Pattern& patt);

private:

   const double _hit_exclus_thresh;
   const double _sig_exclus_thresh;

};

bool HitExclusivity::filter(const art::Event &e, const signature::Pattern& patt)
//bool HitExclusivity::filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
{
    //std::cout << "Testing HitExclusivity" << std::endl;
    this->loadEventHandles(e);
    for (const auto& sig : patt) {
        double sig_q_inclusive = 0.0;
        double sig_q_exclusive = 0.0;
        for (const auto& mcp_s : sig) {
            for (const auto& hit : _mc_hits) {
                auto assmcp = _mcp_bkth_assoc->at(hit.key());
                auto assmdt = _mcp_bkth_assoc->data(hit.key());

                for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
                    auto amd = assmdt[ia];
                    if (assmcp[ia]->TrackId() == mcp_s->TrackId()) {
                        sig_q_inclusive += amd->numElectrons * amd->ideNFraction;
                        if (amd->ideNFraction > _hit_exclus_thresh) 
                            sig_q_exclusive += amd->numElectrons * amd->ideNFraction;
                    }
                }
            }
        }

        if (sig_q_exclusive / sig_q_inclusive < _sig_exclus_thresh)
            return false;
    }

    return true;
}

DEFINE_ART_CLASS_TOOL(HitExclusivity)

}

#endif
