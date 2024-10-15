#ifndef SIGNATURE_NEUTRALKAON_CXX
#define SIGNATURE_NEUTRALKAON_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature
{

class NeutralKaonSignature : public SignatureToolBase 
{

public:
    NeutralKaonSignature(const fhicl::ParameterSet& pset);
    ~NeutralKaonSignature(){};
    
    void configure(fhicl::ParameterSet const & pset) override;

protected:
    void findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature) override;

private:
    art::InputTag _HitProducer;
    art::InputTag _MCPproducer;
    art::InputTag _MCTproducer;
    art::InputTag _BacktrackTag;
};

NeutralKaonSignature::NeutralKaonSignature(const fhicl::ParameterSet& pset)
    : _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
{
    configure(pset); 
}

void NeutralKaonSignature::configure(fhicl::ParameterSet const & pset)
{
    SignatureToolBase::configure(pset);
}

void NeutralKaonSignature::findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature)
{
    auto const &mct_h = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    const simb::MCTruth& mct = mct_h->at(0);
    if (!mct.NeutrinoSet()) return;

    auto const &neutrino = mct.GetNeutrino();
    common::True2RecoMappingXYZ(neutrino.Nu().T(), neutrino.Nu().Vx(), neutrino.Nu().Vy(), neutrino.Nu().Vz(), nu_vtx.data());

    if (neutrino.CCNC() != simb::kCC) return;

    for (const auto &t_part : *mcp_h) {
        if (abs(t_part.PdgCode()) == 13 && t_part.Process() == "primary" && this->aboveThreshold(t_part)) {
            this->fillTrace(t_part, trace_coll);
        }

        if (abs(t_part.PdgCode()) == 310 && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1 && !found_signature) {
            auto daughters = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (daughters.size() == 2) {
                std::vector<int> exp_dtrs = {-211, 211};
                std::vector<int> fnd_dtrs;
                
                for (const auto &dtr : daughters) 
                    fnd_dtrs.push_back(dtr->PdgCode());

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
                            this->fillTrace(dtr, trace_coll);

                        break;
                    }
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(NeutralKaonSignature)
} 

#endif