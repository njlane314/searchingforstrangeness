#ifndef CHARGEDKAONSIGNATURE_CXX
#define CHARGEDKAONSIGNATURE_CXX

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
    void findSignature(art::Event const& evt, SignatureCollection& signature_coll, bool& found_signature) override;

private:
    art::InputTag _MCPproducer;  
};

void ChargedSigmaSignature::findSignature(art::Event const& evt, SignatureCollection& signature_coll, bool& found_signature)
{
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &mc_particle : *mcp_h) 
    {
        int pdg_code = std::abs(mc_particle.PdgCode());

        if (pdg_code == 3112)
        {
            std::cout << "End process: " << mc_particle.EndProcess() << std::endl;
        }

        if (pdg_code == 3112 && mc_particle.Process() == "primary" && (mc_particle.EndProcess() == "Decay" || mc_particle.EndProcess() == "FastScintillation")  && !found_signature) 
        {
            auto daughters = common::GetDaughters(mcp_map.at(mc_particle.TrackId()), mcp_map);
            daughters.erase(std::remove_if(daughters.begin(), daughters.end(), [](const auto& dtr) {
                return dtr->Process() != "Decay";
            }), daughters.end());

            for (const auto& dtr : daughters)
            {
                std::cout << dtr->PdgCode() << std::endl;
                std::cout << dtr->Process() << std::endl;
            }

            std::vector<int> expected_dtrs = std::vector<int>{2112, -211}; // Sigma- -> Neutron + Pi-

            std::vector<int> found_dtrs;
            for (const auto &dtr : daughters) 
                found_dtrs.push_back(dtr->PdgCode());

            std::sort(expected_dtrs.begin(), expected_dtrs.end());
            std::sort(found_dtrs.begin(), found_dtrs.end());

            if (found_dtrs == expected_dtrs) 
            {   
                bool all_above_threshold = std::all_of(daughters.begin(), daughters.end(), [&](const auto& dtr) {
                    return this->aboveThreshold(*dtr);
                });

                if (all_above_threshold) 
                {
                    found_signature = true;

                    this->fillSignature(mcp_map[mc_particle.TrackId()], signature_coll);

                    for (const auto &dtr : daughters) 
                        this->fillSignature(dtr, signature_coll);

                    break;
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(ChargedSigmaSignature)

} 

#endif
