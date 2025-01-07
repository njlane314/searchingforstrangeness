#ifndef SIGNATURE_CHARGEDSIGMA_CXX
#define SIGNATURE_CHARGEDSIGMA_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature {

class ChargedSigmaSignature : public SignatureToolBase 
{
    
public:
    explicit ChargedSigmaSignature(const fhicl::ParameterSet& pset)
    : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    {
        configure(pset); 
    }
    ~ChargedSigmaSignature() {}

    void configure(fhicl::ParameterSet const& pset) override
    {
        SignatureToolBase::configure(pset);
    }

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _MCPproducer;  
};

void ChargedSigmaSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found)
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
        mcp_map[mcp->TrackId()] = mcp;
    }

    auto addDaughterInteractions = [this, &signature, &mcp_map](const art::Ptr<simb::MCParticle>& particle, auto& self) -> void {
        auto daughters = common::GetDaughters(mcp_map.at(particle->TrackId()), mcp_map);
        for (const auto& daugh : daughters) {
            if (daugh->PdgCode() == particle->PdgCode()) {
                this->fillSignature(daugh, signature); 
                self(daugh, self); 
            }
        }
    };

    for (const auto &mcp : *mcp_h)
    {
        if ((abs(mcp.PdgCode()) == 3112 || abs(mcp.PdgCode()) == 3222) && mcp.Process() == "primary" && !signature_found)
        {
            auto particle_iter = mcp_map.find(mcp.TrackId());
            if (particle_iter == mcp_map.end()) 
                continue;

            auto sigma = particle_iter->second; 
            while (sigma->EndProcess() != "Decay")
            {
                if (!this->assessParticle(*sigma))
                    break;

                this->fillSignature(sigma, signature);
                auto daughters = common::GetDaughters(mcp_map.at(sigma->TrackId()), mcp_map);
                if (daughters.empty())
                    break;

                for (const auto& daugh : daughters) {
                    if (daugh->PdgCode() == sigma->PdgCode()) {
                        this->fillSignature(daugh, signature); 
                        sigma = daugh;
                    }
                }
            }

            if (sigma->EndProcess() == "Decay")
            {
                auto decay = common::GetDaughters(mcp_map.at(sigma->TrackId()), mcp_map);
                std::vector<std::vector<int>> decay_modes;
                if (abs(sigma->PdgCode()) == 3222) { // Sigma+
                    decay_modes = {
                        {2212, -211}, // p + pi-
                        {2112, 211}  // n + pi+
                    };
                } else if (abs(sigma->PdgCode()) == 3112) { // Sigma-
                    decay_modes = {
                        {2112, -211} // n + pi-
                    };
                }

                std::vector<int> found_dtrs;
                for (const auto &elem : decay) 
                    found_dtrs.push_back(elem->PdgCode());

                std::sort(found_dtrs.begin(), found_dtrs.end());

                bool valid_decay = false;
                for (const auto& mode : decay_modes) {
                    std::vector<int> sorted_mode = mode;
                    std::sort(sorted_mode.begin(), sorted_mode.end());
                    if (found_dtrs == sorted_mode) {
                        valid_decay = true;
                        break;
                    }
                }

                if (valid_decay) 
                {   
                    bool all_pass = std::all_of(decay.begin(), decay.end(), [&](const auto& daugh) {
                        return this->assessParticle(*daugh);
                    });

                    if (all_pass) 
                    {
                        signature_found = true;

                        this->fillSignature(mcp_map[mcp.TrackId()], signature);
                        for (const auto &elem : decay) 
                        {
                            const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(elem->PdgCode());
                            if (info->Charge() != 0.0) 
                            {
                                this->fillSignature(elem, signature);
                                addDaughterInteractions(elem, addDaughterInteractions);
                            }
                        }

                        break;
                    }
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(ChargedSigmaSignature)

} 

#endif
