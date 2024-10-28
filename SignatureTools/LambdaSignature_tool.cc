#ifndef SIGNATURE_LAMBDA_CXX
#define SIGNATURE_LAMBDA_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature
{

class LambdaSignature : public SignatureToolBase 
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

protected:
    void findSignature(art::Event const& evt, SignatureCollection& signature_coll, bool& found_signature) override;

private:
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
};

void LambdaSignature::findSignature(art::Event const& evt, SignatureCollection& signature_coll, bool& found_signature)
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &t_part : *mcp_h) 
    {
        if (abs(t_part.PdgCode()) == 3122 && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 2 && !found_signature) 
        {
            auto daughters = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (daughters.size() == 2) 
            {
                std::vector<int> exp_dtrs = {-211, 2212};
                std::vector<int> fnd_dtrs;
                
                for (const auto &dtr : daughters) 
                {
                    std::cout << dtr->PdgCode() << std::endl;
                    fnd_dtrs.push_back(dtr->PdgCode());
                }

                std::sort(exp_dtrs.begin(), exp_dtrs.end());
                std::sort(fnd_dtrs.begin(), fnd_dtrs.end());

                if (fnd_dtrs == exp_dtrs) 
                {
                    bool all_above_threshold = std::all_of(daughters.begin(), daughters.end(), [&](const auto& dtr) {
                        return this->aboveThreshold(*dtr);
                    });

                    if (all_above_threshold) 
                    {
                        found_signature = true;
                        for (const auto &dtr : daughters) 
                            this->fillSignature(dtr, signature_coll);

                        break;
                    }
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(LambdaSignature)
} 

#endif
