#ifndef SIGNATURE_KAONSHORT_CXX
#define SIGNATURE_KAONSHORT_CXX

#include <iostream>
#include <optional>
#include "SignatureToolBase.h"
#include "VertexToolBase.h"

namespace signature
{

class KaonShortSignature : public SignatureToolBase, public VertexToolBase
{
public:
    explicit KaonShortSignature(const fhicl::ParameterSet& pset)
        : _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
        , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
        , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
        , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    {
        configure(pset);
    }

    ~KaonShortSignature() override = default;

    void configure(fhicl::ParameterSet const& pset) override
    {
        SignatureToolBase::configure(pset);
    }

    Type getSignatureType() const override 
    {
        return SignatureKaonShort;
    }

    TVector3 findVertex(art::Event const& evt) const override;

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _HitProducer;
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
    art::InputTag _BacktrackTag;
};

void KaonShortSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found)
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
        if (abs(mcp.PdgCode()) == 311 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" && mcp.NumberDaughters() == 1 && !signature_found) 
        {
            auto dtrs = common::GetDaughters(mcp_map.at(mcp.TrackId()), mcp_map);
            if (dtrs.size() != 1) continue; 

            auto dtr = dtrs.at(0);
            if (dtr->PdgCode() == 310 && dtr->Process() == "Decay" && dtr->EndProcess() == "Decay" && dtr->NumberDaughters() == 2 && !signature_found)
            {
                auto decay = common::GetDaughters(mcp_map.at(dtr->TrackId()), mcp_map);
                if (decay.size() != 2) continue;
                
                std::vector<int> exp_decay = {-211, 211};
                std::vector<int> fnd_decay;
                
                for (const auto &elem : decay) 
                    fnd_decay.push_back(elem->PdgCode());

                std::sort(exp_decay.begin(), exp_decay.end());
                std::sort(fnd_decay.begin(), fnd_decay.end());

                if (fnd_decay == exp_decay) 
                {
                    signature_found = true;
                    for (const auto &elem : decay) 
                    {
                        const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(elem->PdgCode());
                        if (info->Charge() == 0.0) 
                            continue;

                        this->fillSignature(elem, signature);
                        addDaughterInteractions(elem, addDaughterInteractions);
                    }
                }
            }
        }
    }
}

TVector3 KaonShortSignature::findVertex(art::Event const& evt) const
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) 
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &mcp : *mcp_h) 
    {
        if (abs(mcp.PdgCode()) == 311 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" && mcp.NumberDaughters() == 1) {
            auto dtrs = common::GetDaughters(mcp_map.at(mcp.TrackId()), mcp_map);
            if (dtrs.size() != 1) continue; 

            auto dtr = dtrs.at(0);
            if (dtr->PdgCode() == 310 && dtr->Process() == "Decay" && dtr->EndProcess() == "Decay" && dtr->NumberDaughters() == 2)
            {
                const TLorentzVector& end_position = dtr->EndPosition();
                
                return TVector3(end_position.X(), end_position.Y(), end_position.Z());
            }
        }
    }

    return TVector3();
}

DEFINE_ART_CLASS_TOOL(KaonShortSignature)
} 

#endif