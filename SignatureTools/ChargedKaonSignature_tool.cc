#ifndef CHARGEDKAONSIGNATURE_CXX
#define CHARGEDKAONSIGNATURE_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature {

class ChargedKaonSignature : public SignatureToolBase 
{
    
public:
    ChargedKaonSignature(const fhicl::ParameterSet& pset);
    ~ChargedKaonSignature() {}

    void configure(fhicl::ParameterSet const& pset) override;

protected:
    void findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature) override;

private:
    art::InputTag _MCPproducer;  
    std::string _decay_mode;
};

ChargedKaonSignature::ChargedKaonSignature(const fhicl::ParameterSet& pset)
    : _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _decay_mode{pset.get<std::string>("DecayMode", "muonic")} 
{
    configure(pset); 
}

void ChargedKaonSignature::configure(fhicl::ParameterSet const& pset)
{
    SignatureToolBase::configure(pset);
}

void ChargedKaonSignature::findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature)
{
    std::cout << "Looking for charged-kaon signature..." << std::endl;
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    for (const auto &mc_particle : *mcp_h) 
    {
        int pdg_code = std::abs(mc_particle.PdgCode());

        if (pdg_code == 321 && mc_particle.Process() == "primary" && (mc_particle.EndProcess() == "Decay" || mc_particle.EndProcess() == "FastScintillation")  && !found_signature) 
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

            std::vector<int> expected_dtrs;
            if (_decay_mode == "muonic")  
            {
                expected_dtrs = (mc_particle.PdgCode() == 321) ? std::vector<int>{-13, +14}  // K+ -> Muon+ + Neutrino
                                                    : std::vector<int>{+13, -14}; // K- -> Muon- + Antineutrino
            }
            else if (_decay_mode == "pionic") 
            {
                expected_dtrs = (mc_particle.PdgCode() == 321) ? std::vector<int>{211, 111}  // K+ -> Pi+ + Pi0
                                                    : std::vector<int>{-211, 111}; // K- -> Pi- + Pi0
            }

            std::vector<int> found_dtrs;
            for (const auto &dtr : daughters) 
                found_dtrs.push_back(dtr->PdgCode());

            std::sort(expected_dtrs.begin(), expected_dtrs.end());
            std::sort(found_dtrs.begin(), found_dtrs.end());

            if (found_dtrs == expected_dtrs) 
            {   
                std::cout << "Found dtrs = exp" << std::endl;
                bool all_above_threshold = std::all_of(daughters.begin(), daughters.end(), [&](const auto& dtr) {
                    return this->aboveThreshold(*dtr);
                });

                if (all_above_threshold) 
                {
                    std::cout << "Found kaon signature" << std::endl;
                    found_signature = true;
                    for (const auto &dtr : daughters) 
                        this->fillTrace(dtr, trace_coll);

                    break;
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(ChargedKaonSignature)

} 

#endif
