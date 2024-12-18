#ifndef SIGNATURE_LAMBDA_CXX
#define SIGNATURE_LAMBDA_CXX

#include <iostream>
#include "SignatureToolBase.h"
#include "VertexToolBase.h"

namespace signature
{

class LambdaSignature : public SignatureToolBase, public VertexToolBase
{
public:
    explicit LambdaSignature(const fhicl::ParameterSet& pset)
        : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
        , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    {
        configure(pset);
    }

    ~LambdaSignature() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        SignatureToolBase::configure(pset);
    }

    TVector3 findVertex(art::Event const& evt) const override;

protected:
    void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;

private:
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
};

void LambdaSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found)
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &mcp : *mcp_h) 
    {
        if (abs(mcp.PdgCode()) == 3122 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" && mcp.NumberDaughters() == 2 && !signature_found) 
        {
            auto dtrs = common::GetDaughters(mcp_map.at(mcp.TrackId()), mcp_map);
            if (dtrs.size() == 2) 
            {
                std::vector<int> exp_decay = {-211, 2212};
                std::vector<int> fnd_decay;
                
                for (const auto &dtr : dtrs) 
                    fnd_decay.push_back(dtr->PdgCode());

                std::sort(exp_decay.begin(), exp_decay.end());
                std::sort(fnd_decay.begin(), fnd_decay.end());

                if (fnd_decay == exp_decay) 
                {
                    bool all_pass = std::all_of(dtrs.begin(), dtrs.end(), [&](const auto& dtr) {
                        return this->assessParticle(*dtr);
                    });

                    if (all_pass) 
                    {
                        signature_found = true;
                        for (const auto &dtr : dtrs) 
                        {
                            const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(dtr->PdgCode());
                            if (info->Charge() != 0.0) 
                                this->fillSignature(dtr, signature);
                        }

                        break;
                    }
                }
            }
        }
    }
}

TVector3 LambdaSignature::findVertex(art::Event const& evt) const
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) 
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &mcp : *mcp_h) {
        if (abs(mcp.PdgCode()) == 3122 && mcp.Process() == "primary" && mcp.EndProcess() == "Decay" && mcp.NumberDaughters() == 2) 
        {
            const TLorentzVector& end_position = mcp.EndPosition();
            return TVector3(end_position.X(), end_position.Y(), end_position.Z());
        }
    }

    return TVector3();
}

DEFINE_ART_CLASS_TOOL(LambdaSignature)
} 

#endif
