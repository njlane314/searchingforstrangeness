#ifndef LAMBDA_SIGNATURE_INTEGRITY_H
#define LAMBDA_SIGNATURE_INTEGRITY_H

#include "SignatureIntegrity.h" 

namespace signature {

class LambdaSignatureIntegrity : public SignatureIntegrity {
public:
    explicit LambdaSignatureIntegrity(const fhicl::ParameterSet& pset) :
      SignatureIntegrity{(pset)} {
        configure(pset);
    }
    ~LambdaSignatureIntegrity() override = default;
    bool filter(const art::Event &e, const Signature& sig, const SignatureType& type, common::PandoraView view) override;
};

bool LambdaSignatureIntegrity::filter(const art::Event &e, const Signature& sig, const SignatureType& type, common::PandoraView view) {
    if(type != kLambdaSignature) 
        return true;
    this->loadEventHandles(e,view);
    std::vector<int> trackids;
    for (const auto& mcp_s : sig) 
        trackids.push_back(mcp_s->TrackId());
    for (const auto& mcp_s : sig){
        if(std::find(trackids.begin(), trackids.end(), mcp_s->Mother()) != trackids.end()) 
            continue;
        if(!checkStart2(mcp_s,view)) 
            return false;
    }
    return true;
}

DEFINE_ART_CLASS_TOOL(LambdaSignatureIntegrity)

}

#endif