#ifndef SIGNATURE_MUON_CXX
#define SIGNATURE_MUON_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature {

class MuonSignature : public SignatureToolBase 
{
public:
    explicit MuonSignature(const fhicl::ParameterSet& pset) 
        : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    {
        configure(pset); 
    }

    ~MuonSignature() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override 
    {
        SignatureToolBase::configure(pset);
    }

    Type getSignatureType() const override 
    {
        return SignaturePrimaryMuon;
    }

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _MCPproducer;
};

void MuonSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found)
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::vector<art::Ptr<simb::MCParticle>> mcp_vec;
    art::fill_ptr_vector(mcp_vec, mcp_h);

    for (const auto& mcp : mcp_vec) 
    {
        if (std::abs(mcp->PdgCode()) == 13 && mcp->Process() == "primary") 
        {
            signature_found = true;
            this->fillSignature(mcp, signature);
        }
    }
}

DEFINE_ART_CLASS_TOOL(MuonSignature)

}

#endif
