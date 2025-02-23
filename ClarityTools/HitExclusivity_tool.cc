#ifndef HIT_EXCLUSIVITY_H
#define HIT_EXCLUSIVITY_H

#include "ClarityToolBase.h" 

namespace signature {

class HitExclusivity : ClarityToolBase {
public:
    explicit HitExclusivity(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)}
    , _hit_exclus_thresh{pset.get<double>("HitExclusivityThreshold", 0.80)}
    , _sig_exclus_thresh{pset.get<double>("SignatureExclusivityThreshold", 0.70)} {
        configure(pset);
    }
    ~HitExclusivity() override = default;
    void configure(fhicl::ParameterSet const& pset) override {
        ClarityToolBase::configure(pset);
    }
    bool filter(const art::Event &e, const Signature& sig, const SignatureType& type, common::PandoraView view);

private:
    const double _hit_exclus_thresh;
    const double _sig_exclus_thresh;
};

bool HitExclusivity::filter(const art::Event &e, const Signature& sig, const SignatureType& type, common::PandoraView view) {
    if(!this->loadEventHandles(e,view)) 
        return false;
    for (const auto& mcp_s : sig) {
        double sig_q_inclusive = 0.0;
        double sig_q_exclusive = 0.0;
        for (const auto& hit : _mc_hits) {
            auto assmcp = _mcp_bkth_assoc->at(hit.key());
            auto assmdt = _mcp_bkth_assoc->data(hit.key());
            for (unsigned int ia = 0; ia < assmcp.size(); ++ia) {
                auto amd = assmdt[ia];
                if (assmcp[ia]->TrackId() == mcp_s->TrackId()) {
                    sig_q_inclusive += amd->numElectrons * amd->ideNFraction;
                    double e_frac = 0.0;
                    for (unsigned int ia2 = 0; ia2 < assmcp.size(); ++ia2){
                        if(abs(assmcp[ia2]->PdgCode()) == 11 && assmcp[ia2]->Process() != "primary") e_frac += assmdt[ia2]->ideNFraction;
                    }          
                    if (amd->ideNFraction/(1.0 - e_frac) > _hit_exclus_thresh){
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