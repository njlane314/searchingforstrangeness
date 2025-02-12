#ifndef KPLUS_SIGNATURE_INTEGRITY_H
#define KPLUS_SIGNATURE_INTEGRITY_H

#include "SignatureIntegrity.h" 

namespace claritytools {

class KPlusSignatureIntegrity : public SignatureIntegrity {

public:
    explicit KPlusSignatureIntegrity(const fhicl::ParameterSet& pset) :
      SignatureIntegrity{(pset)} 
    {
        configure(pset);
    }

    ~KPlusSignatureIntegrity() override = default;
   
    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view) override;

};

bool KPlusSignatureIntegrity::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{
  if(sig.first != signature::SignatureChargedKaon) return true;

  this->loadEventHandles(e,view);
 
  // Check the start and end of the kaon, and only the start of the muon/pion it decays to 
  for (const auto& mcp_s : sig.second){
    if(!checkDeadChannelFrac(mcp_s,view)) return false;
    if(abs(mcp_s->PdgCode()) == 321 && (!checkStart2(mcp_s,view) || !checkEnd2(mcp_s,view))) return false;
    else if((abs(mcp_s->PdgCode()) == 13 || abs(mcp_s->PdgCode()) == 211) && !checkStart2(mcp_s,view)) return false;
  }

  return true;

}

DEFINE_ART_CLASS_TOOL(KPlusSignatureIntegrity)

}

#endif
