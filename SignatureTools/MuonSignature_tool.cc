#ifndef SIGNATURE_MUON_CXX
#define SIGNATURE_MUON_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature {

class MuonSignature : public SignatureToolBase {
public:
    explicit MuonSignature(const fhicl::ParameterSet& pset)
        : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")} {
        configure(pset); 
    }
    ~MuonSignature() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override {
        SignatureToolBase::configure(pset);
        _mu_thresh_mom = pset.get<float>("MuonThreshMomentum", 0.1);
    }

    SignatureType getSignatureType() const override {
        return kPrimaryMuonSignature;
    }

    bool isDetectable(art::Event const& evt, Signature const& signature) const override;

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _MCPproducer;
    float _mu_thresh_mom;
};

void MuonSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found) {
    auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::vector<art::Ptr<simb::MCParticle>> mcp_vec;
    art::fill_ptr_vector(mcp_vec, mcp_h);
    for (const auto& mcp : mcp_vec) {
        if (std::abs(mcp->PdgCode()) == 13 && mcp->Process() == "primary") {
            signature_found = true;
            fillSignature(mcp, signature);
        }
    }
}

bool MuonSignature::isDetectable(art::Event const& evt, Signature const& signature) const {
    for (const auto& mcp : signature) {
        int pdg_code = std::abs(mcp->PdgCode());
        float momentum = mcp->P();
        std::string process = mcp->Process();
        if (pdg_code == 13 && process == "primary" && momentum < _mu_thresh_mom)
            return false;
    }
    return true;
}

DEFINE_ART_CLASS_TOOL(MuonSignature)

} 

#endif