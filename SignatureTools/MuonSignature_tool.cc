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
        , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    {
        configure(pset); 
    }

    ~MuonSignature() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override 
    {
        SignatureToolBase::configure(pset);
    }

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& found_signature) override;

private:
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
};

void MuonSignature::findSignature(art::Event const& evt, Signature& signature, bool& found_signature)
{
    found_signature = false;
    auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    for (size_t i = 0; i < mcp_h->size(); ++i) 
    {
        const simb::MCParticle& mcp = mcp_h->at(i);
        if (std::abs(mcp.PdgCode()) == 13 && mcp.Process() == "primary" && this->aboveThreshold(mcp)) 
        {
            this->fillSignature(art::Ptr<simb::MCParticle>(mcp_h, i), signature);
            found_signature = true;
            break;
        }
    }
}

DEFINE_ART_CLASS_TOOL(MuonSignature)

}

#endif
