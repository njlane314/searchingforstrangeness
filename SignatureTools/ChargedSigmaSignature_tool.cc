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

    auto addDaughterInteractions = [this, &signature, &mcp_map](const art::Ptr<simb::MCParticle>& particle, auto& self) -> art::Ptr<simb::MCParticle> {
        auto daughters = common::GetDaughters(mcp_map.at(particle->TrackId()), mcp_map);
        art::Ptr<simb::MCParticle> end_particle = particle;
        for (const auto& daugh : daughters) {
            if (daugh->PdgCode() == particle->PdgCode() && this->assessParticle(*daugh)) {
                this->fillSignature(daugh, signature); 
                end_particle = self(daugh, self); 
            }
        }
        
        return end_particle;
    };

    for (const auto &mcp : *mcp_h)
    {
        if ((abs(mcp.PdgCode()) == 3112 || abs(mcp.PdgCode()) == 3222) && mcp.Process() == "primary" && !signature_found)
        {
            auto particle_iter = mcp_map.find(mcp.TrackId());
            if (particle_iter == mcp_map.end()) 
                continue;

            auto sigma = particle_iter->second; 
            if (!this->assessParticle(*sigma))
                break;

            this->fillSignature(sigma, signature);
            auto end_sigma = addDaughterInteractions(sigma, addDaughterInteractions);

            std::cout << end_sigma->EndProcess() << std::endl;
            if (end_sigma->EndProcess() == "Decay")
            {
                auto decay = common::GetDaughters(mcp_map.at(sigma->TrackId()), mcp_map);
                std::vector<int> found_dtrs;
                std::vector<art::Ptr<simb::MCParticle>> clean_decay;

                for (const auto& elem : decay)
                {
                    std::cout << elem->PdgCode() << std::endl;
                    if (std::abs(elem->PdgCode()) == 11 || std::abs(elem->PdgCode()) == 22 || std::abs(elem->PdgCode()) >= 1000000000) 
                        continue;

                    found_dtrs.push_back(elem->PdgCode());
                    clean_decay.push_back(elem);
                }
                
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

                std::cout << valid_decay << std::endl;

                if (valid_decay) 
                {
                    bool all_pass = std::all_of(clean_decay.begin(), clean_decay.end(), [&](const auto& elem) {
                        return this->assessParticle(*elem);
                    });

                    std::cout << all_pass << std::endl;

                    if (all_pass) 
                    {
                        signature_found = true;

                        this->fillSignature(mcp_map[mcp.TrackId()], signature);
                        for (const auto& elem : clean_decay) 
                        {
                            const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(elem->PdgCode());
                            if (info->Charge() != 0.0) 
                            {
                                std::cout << "Filling signature" << std::endl;
                                this->fillSignature(elem, signature);
                                auto end_particle = addDaughterInteractions(elem, addDaughterInteractions);
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
