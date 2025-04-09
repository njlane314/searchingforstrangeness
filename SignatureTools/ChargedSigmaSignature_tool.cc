#ifndef SIGNATURE_CHARGEDSIGMA_CXX
#define SIGNATURE_CHARGEDSIGMA_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature 
{
    class ChargedSigmaSignature : public SignatureToolBase {
    public:
        explicit ChargedSigmaSignature(const fhicl::ParameterSet& pset) {
            configure(pset);
        }
        
        ~ChargedSigmaSignature() override = default;
        
        void configure(fhicl::ParameterSet const& pset) override {
            SignatureToolBase::configure(pset);
        }
        
        SignatureType getSignatureType() const override {
            return kChargedSigmaSignature;
        }
        
    protected:
        void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;
    };

    void ChargedSigmaSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found) {
        auto mcp_map = buildParticleMap(evt);
        auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
        
        for (const auto& mcp : *mcp_h) {
            if ((abs(mcp.PdgCode()) == 3112 || abs(mcp.PdgCode()) == 3222) && 
                mcp.Process() == "primary" && !signature_found) {
                
                auto particle_iter = mcp_map.find(mcp.TrackId());
                if (particle_iter == mcp_map.end()) 
                    continue;
                    
                auto sigma = particle_iter->second;
                fillSignature(sigma, signature);
        
                art::Ptr<simb::MCParticle> end_sigma = sigma;
                bool found_end = false;
                
                traverseDecayChain(sigma, mcp_map, [&end_sigma, &found_end, sigma_pdg = sigma->PdgCode()]
                                 (const art::Ptr<simb::MCParticle>& p) {
                    if (p->PdgCode() == sigma_pdg && p->EndProcess() == "Decay" && !found_end) {
                        end_sigma = p;
                        found_end = true;
                    }
                });
                
                if (end_sigma->EndProcess() == "Decay") {
                    auto decay = common::GetDaughters(end_sigma, mcp_map);
                    std::vector<int> found_dtrs;
                    std::vector<art::Ptr<simb::MCParticle>> clean_decay;
                    
                    for (const auto& elem : decay) {
                        if (std::abs(elem->PdgCode()) == 11 || 
                            std::abs(elem->PdgCode()) == 22 || 
                            std::abs(elem->PdgCode()) >= 1000000000) 
                            continue;
                        found_dtrs.push_back(elem->PdgCode());
                        clean_decay.push_back(elem);
                    }
                    
                    std::vector<std::vector<int>> decay_modes;
                    if (abs(sigma->PdgCode()) == 3222) { // Sigma+
                        decay_modes = {
                            {2212, 111}, // p + pi0
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
                    
                    if (valid_decay) {
                        signature_found = true;
                        for (const auto& elem : clean_decay) {
                            const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(elem->PdgCode());
                            if (info->Charge() == 0.0) 
                                continue;
                                
                            fillSignature(elem, signature);
                            traverseDecayChain(elem, mcp_map, [this, &signature, elem_pdg = elem->PdgCode()]
                                             (const art::Ptr<simb::MCParticle>& p) {
                                if (p->PdgCode() == elem_pdg) {
                                    const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(p->PdgCode());
                                    if (info && info->Charge() != 0.0) {
                                        fillSignature(p, signature);
                                    }
                                }
                            });
                        }
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ChargedSigmaSignature)
}

#endif