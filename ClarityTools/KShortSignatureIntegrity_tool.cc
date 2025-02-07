#ifndef KSHORT_SIGNATURE_INTEGRITY_H
#define KSHORT_SIGNATURE_INTEGRITY_H

#include "SignatureIntegrity.h" 

namespace claritytools {

class KShortSignatureIntegrity : public SignatureIntegrity {

public:
    explicit KShortSignatureIntegrity(const fhicl::ParameterSet& pset) :
      SignatureIntegrity{(pset)} 
    {
        configure(pset);
    }

    ~KShortSignatureIntegrity() override = default;
   
    /* 
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }
    */

    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view) override;

};

bool KShortSignatureIntegrity::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{
  // Only check mcps that are children of a KShort
  if(sig.first != signature::SignatureKaonShort) return true;

  std::cout << "Checking KShort Signature Integrity" << std::endl;

  this->loadEventHandles(e,view);

  // Only check the first two particles - don't care about secondaries
  std::vector<int> trackids;
  for (const auto& mcp_s : sig.second) trackids.push_back(mcp_s->TrackId());

  for (const auto& mcp_s : sig.second){
    if(std::find(trackids.begin(),trackids.end(),mcp_s->Mother()) != trackids.end()) continue;
    if(!checkStart(mcp_s,view)) return false;
  }

  return true;

}

DEFINE_ART_CLASS_TOOL(KShortSignatureIntegrity)

}

#endif
