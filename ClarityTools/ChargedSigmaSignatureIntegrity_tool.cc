#ifndef KPLUS_SIGNATURE_INTEGRITY_H
#define KPLUS_SIGNATURE_INTEGRITY_H

#include "SignatureIntegrity.h" 

namespace signature {

class ChargedSigmaSignatureIntegrity : public SignatureIntegrity {
public:
    explicit ChargedSigmaSignatureIntegrity(const fhicl::ParameterSet& pset) :
      SignatureIntegrity{(pset)} {
        configure(pset);
    }
    ~ChargedSigmaSignatureIntegrity() override = default;
   
    bool filter(const art::Event &e, const Signature& sig, const SignatureType& type, common::PandoraView view) override;
};

bool ChargedSigmaSignatureIntegrity::filter(const art::Event &e, const Signature& sig, const SignatureType& type, common::PandoraView view) {
    if(type != kChargedSigmaSignature) 
        return true;
    this->loadEventHandles(e,view);
    for (const auto& mcp_s : sig){
        if(!checkDeadChannelFrac(mcp_s,view)) 
            return false;
        if(abs(mcp_s->PdgCode()) == 3112 && (!checkStart2(mcp_s,view) || !checkEnd2(mcp_s,view))) 
            return false;
        else if((abs(mcp_s->PdgCode()) == 211) && !checkStart2(mcp_s,view)) 
            return false;
    }
    return true;
}

DEFINE_ART_CLASS_TOOL(ChargedSigmaSignatureIntegrity)
}

#endif