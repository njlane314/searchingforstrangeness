#ifndef HIT_EXCLUSIVITY_H
#define HIT_EXCLUSIVITY_H

#include "ClarityToolBase.h" 

namespace claritytools {

class HitExclusivity : ClarityToolBase {

public:
    explicit HitExclusivity(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)}
    , _hit_exclus_thresh{pset.get<double>("HitExclusivityThreshold", 0.80)}
    , _sig_exclus_thresh{pset.get<double>("SignatureExclusivityThreshold", 0.70)}
    {
        configure(pset);
    }

    ~HitExclusivity() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view);

private:

   const double _hit_exclus_thresh;
   const double _sig_exclus_thresh;

};

bool HitExclusivity::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{

  if(_verbose)
    std::cout << "Checking HitExclusivity in view " << view << " for signature " << signature::GetSignatureName(sig) << std::endl;

  if(!this->loadEventHandles(e,view)) return false;

  for (const auto& mcp_s : sig.second) {

    if(_verbose)
      std::cout << "Checking HitExclusivity for particle pdg=" << mcp_s->PdgCode() << " trackid=" << mcp_s->TrackId() << std::endl;

    double sig_q_inclusive = 0.0;
    double sig_q_exclusive = 0.0;

    for (const auto& hit : _mc_hits) {
      auto assmcp = _mcp_bkth_assoc->at(hit.key());
      auto assmdt = _mcp_bkth_assoc->data(hit.key());

      for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
        auto amd = assmdt[ia];
        if (assmcp[ia]->TrackId() == mcp_s->TrackId()) {
          sig_q_inclusive += amd->numElectrons * amd->ideNFraction;

          // Calculate the fraction of energy due to delta rays
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
  
    if(_verbose)
      std::cout << "sig_q_exclusive/sig_q_inclusive = " << sig_q_exclusive << "/" << sig_q_inclusive << "=" << sig_q_exclusive / sig_q_inclusive << std::endl;

    if (sig_q_exclusive / sig_q_inclusive < _sig_exclus_thresh){
      if(_verbose)
        std::cout << "Signature " << signature::GetSignatureName(sig) << " failed HitExclusivity in plane " << view << std::endl;
      return false;
    }

  }

  if(_verbose)
    std::cout << "Signature " << signature::GetSignatureName(sig) << " passed HitExclusivity in plane " << view << std::endl;

  return true;

}

DEFINE_ART_CLASS_TOOL(HitExclusivity)

}

#endif
