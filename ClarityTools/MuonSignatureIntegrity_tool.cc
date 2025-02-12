#ifndef MUON_SIGNATURE_INTEGRITY_CC
#define MUON_SIGNATURE_INTEGRITY_CC

#include "SignatureIntegrity.h" 

namespace claritytools {

class MuonSignatureIntegrity : public SignatureIntegrity {

public:
    explicit MuonSignatureIntegrity(const fhicl::ParameterSet& pset) :
      SignatureIntegrity{(pset)} 
    {
        configure(pset);
    }

    ~MuonSignatureIntegrity() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
    {   
      if(sig.first != signature::SignaturePrimaryMuon) return true;
      this->loadEventHandles(e,view);

      for (const auto& mcp_s : sig.second) {
        if(!checkStart(mcp_s,view)) return false;
      }
      return true;
    }

};

DEFINE_ART_CLASS_TOOL(MuonSignatureIntegrity)

}

#endif
