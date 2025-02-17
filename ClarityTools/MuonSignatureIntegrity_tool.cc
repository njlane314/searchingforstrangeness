#ifndef MUON_SIGNATURE_INTEGRITY_CC
#define MUON_SIGNATURE_INTEGRITY_CC

#include "SignatureIntegrity.h" 

namespace signature {

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

    bool filter(const art::Event &e, const Signature& sig, const Type& type, common::PandoraView view)
    {   
      if(type != signature::SignaturePrimaryMuon) return true;
      this->loadEventHandles(e,view);

      for (const auto& mcp_s : sig) {
        if(!checkStart2(mcp_s,view)) return false;
      }
      return true;
    }

};

DEFINE_ART_CLASS_TOOL(MuonSignatureIntegrity)

}

#endif