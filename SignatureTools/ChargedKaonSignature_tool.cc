#ifndef SIGNATURE_CHARGEDKAON_CXX
#define SIGNATURE_CHARGEDKAON_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature {

class ChargedKaonSignature : public SignatureToolBase {
public:
    explicit ChargedKaonSignature(const fhicl::ParameterSet& pset)
    : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")} {
        configure(pset); 
    }
    ~ChargedKaonSignature() {}

    void configure(fhicl::ParameterSet const& pset) override {
        SignatureToolBase::configure(pset);
    }

    SignatureType getSignatureType() const override {
        return kChargedKaonSignature;
    }

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _MCPproducer;  
};

void ChargedKaonSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found) {
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }
    auto addDaughterInteractions = [this, &signature, &mcp_map](const art::Ptr<simb::MCParticle>& particle, auto& self) -> art::Ptr<simb::MCParticle> {
        auto daughters = common::GetDaughters(mcp_map.at(particle->TrackId()), mcp_map);
        art::Ptr<simb::MCParticle> end_particle = particle;
        for (const auto& daugh : daughters) {
            if (daugh->PdgCode() == particle->PdgCode()) {
                this->fillSignature(daugh, signature); 
                end_particle = self(daugh, self); 
            }
        }
        return end_particle;
    };
    for (const auto &mcp : *mcp_h) {
        if (std::abs(mcp.PdgCode()) == 321 && mcp.Process() == "primary"  && !signature_found) {
            auto particle_iter = mcp_map.find(mcp.TrackId());
            if (particle_iter == mcp_map.end()) 
                continue;
            auto kaon = particle_iter->second;
            this->fillSignature(kaon, signature);
            auto end_kaon = addDaughterInteractions(kaon, addDaughterInteractions);
            if (end_kaon->EndProcess() == "Decay" || end_kaon->EndProcess() == "FastScintillation") {
                auto decay = common::GetDaughters(mcp_map.at(kaon->TrackId()), mcp_map);
                std::vector<int> found_dtrs;
                std::vector<art::Ptr<simb::MCParticle>> clean_decay;
                for (const auto& elem : decay) {
                    if (std::abs(elem->PdgCode()) == 11 || std::abs(elem->PdgCode()) == 22 || std::abs(elem->PdgCode()) >= 1000000000) 
                        continue;
                    found_dtrs.push_back(elem->PdgCode());
                    clean_decay.push_back(elem);
                }
                std::vector<std::vector<int>> decay_modes;
                if (kaon->PdgCode() == 321) {  // K+
                    decay_modes = {
                        {-13, +14},   // mu+ + neutrino
                        {211, 111}    // pi+ + pi0
                    };
                } else if (kaon->PdgCode() == -321) {  // K-
                    decay_modes = {
                        {+13, -14},   // mu- + antineutrino
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
                    this->fillSignature(mcp_map[mcp.TrackId()], signature);
                    for (const auto& elem : clean_decay) {
                        const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(elem->PdgCode());
                        if (info->Charge() == 0.0) 
                            continue;
                        this->fillSignature(elem, signature);
                        auto end_particle = addDaughterInteractions(elem, addDaughterInteractions);
                    }
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(ChargedKaonSignature)

} 

#endif
