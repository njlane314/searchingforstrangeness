#ifndef LAMBDA_SIGNATURE_INTEGRITY_H
#define LAMBDA_SIGNATURE_INTEGRITY_H

#include "SignatureIntegrity.h" 

namespace claritytools {

class LambdaSignatureIntegrity : public SignatureIntegrity {

public:
    explicit LambdaSignatureIntegrity(const fhicl::ParameterSet& pset) :
      SignatureIntegrity{(pset)} 
    {
        configure(pset);
    }

    ~LambdaSignatureIntegrity() override = default;
   
    /* 
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }
    */

    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view) override;

};

bool LambdaSignatureIntegrity::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{
  std::cout << "Checking Signature Integrity" << std::endl;

  // Only check mcps that are children of a Lambda
  if(sig.first != signature::SignatureLambda) return true;

  this->loadEventHandles(e,view);

  for (const auto& mcp_s : sig.second){
    if(!checkStart(mcp_s,view)) return false;
  }

  return true;

}

DEFINE_ART_CLASS_TOOL(LambdaSignatureIntegrity)

}

#endif
