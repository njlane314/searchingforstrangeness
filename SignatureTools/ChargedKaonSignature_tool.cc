#ifndef SIGNATURE_CHARGEDKAON_CXX
#define SIGNATURE_CHARGEDKAON_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature 
{
    class ChargedKaonSignature : public SignatureToolBase {
    public:
        explicit ChargedKaonSignature(const fhicl::ParameterSet& pset) {
            configure(pset);
        }
        
        ~ChargedKaonSignature() override = default;
        
        void configure(fhicl::ParameterSet const& pset) override {
            SignatureToolBase::configure(pset);
        }
        
        SignatureType getSignatureType() const override {
            return kChargedKaonSignature;
        }
        
    protected:
        void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;
    };

    void ChargedKaonSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found) {
        auto mcp_map = buildParticleMap(evt);
        auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
        
        for (const auto& mcp : *mcp_h) {
            if (std::abs(mcp.PdgCode()) == 321 && mcp.Process() == "primary" && !signature_found) {
                auto particle_iter = mcp_map.find(mcp.TrackId());
                if (particle_iter == mcp_map.end()) 
                    continue;
                    
                auto kaon = particle_iter->second;
                fillSignature(kaon, signature);
                
                art::Ptr<simb::MCParticle> end_kaon = kaon;
                bool found_end = false;
                int kaon_pdg = kaon->PdgCode();
                
                traverseDecayChain(kaon, mcp_map, [&end_kaon, &found_end, kaon_pdg]
                                 (const art::Ptr<simb::MCParticle>& p) {
                    if (p->PdgCode() == kaon_pdg && 
                        (p->EndProcess() == "Decay" || p->EndProcess() == "FastScintillation") && 
                        !found_end) {
                        end_kaon = p;
                        found_end = true;
                    }
                });
                
                if (end_kaon->EndProcess() == "Decay" || end_kaon->EndProcess() == "FastScintillation") {
                    std::vector<art::Ptr<simb::MCParticle>> decay = common::GetDaughters(end_kaon, mcp_map);
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
                    if (kaon->PdgCode() == 321) {  // K+
                        decay_modes = {
                            {-13, 14},   // mu+ + neutrino
                            {211, 111}    // pi+ + pi0
                        };
                    } else if (kaon->PdgCode() == -321) {  // K-
                        decay_modes = {
                            {13, -14},   // mu- + antineutrino
                            {-211, 111}   // pi- + pi0
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
                            int elem_pdg = elem->PdgCode();
                            
                            traverseDecayChain(elem, mcp_map, [this, &signature, elem_pdg]
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

    DEFINE_ART_CLASS_TOOL(ChargedKaonSignature)
}

#endif