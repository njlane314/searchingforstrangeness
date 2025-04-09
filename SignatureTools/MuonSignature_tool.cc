#ifndef SIGNATURE_MUON_CXX
#define SIGNATURE_MUON_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature 
{
    class MuonSignature : public SignatureToolBase {
    public:
        explicit MuonSignature(const fhicl::ParameterSet& pset) {
            configure(pset);
        }
        
        ~MuonSignature() override = default;
        
        void configure(fhicl::ParameterSet const& pset) override {
            SignatureToolBase::configure(pset);
            _muThreshMom = pset.get<float>("MuonThreshMomentum", 0.1);
        }
        
        SignatureType getSignatureType() const override {
            return kPrimaryMuonSignature;
        }
        
        bool isDetectable(art::Event const& evt, Signature const& signature) const override;
        
    protected:
        void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;
        
    private:
        float _muThreshMom;
    };

    void MuonSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found) {
        auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
        
        for (size_t i = 0; i < mcp_h->size(); ++i) {
            const art::Ptr<simb::MCParticle> mcp(mcp_h, i);
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
            
            if (pdg_code == 13 && process == "primary" && momentum < _muThreshMom)
                return false;
        }
        return true;
    }

    DEFINE_ART_CLASS_TOOL(MuonSignature)
}

#endif