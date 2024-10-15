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
    
    void configure(fhicl::ParameterSet const & pset);
    void identifySignalParticles(art::Event const& evt, 
                                         std::vector<Particle>& particle_elements, 
                                         bool& found_signature, 
                                         std::array<float, 3>& nu_vtx);
    
private:
    
};

NeutralKaonSignature::NeutralKaonSignature(const fhicl::ParameterSet& pset)
    : _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _PionThreshold{pset.get<float>("PionThreshold", 0.1)}
    , _MuonThreshold{pset.get<float>("MuonThreshold", 0.1)}
    , _ProtonThreshold{pset.get<float>("ProtonThreshold", 0.1)}
{
}

void NeutralKaonSignature::configure(fhicl::ParameterSet const & pset)
{
}

void NeutralKaonSignature::identifySignalParticles(art::Event const& evt, 
                                         std::vector<Particle>& particle_elements, 
                                         bool& found_signature, 
                                         std::array<float, 3>& nu_vtx)
{
    particle_elements.clear();
    found_signature = false;
    std::fill(nu_vtx.begin(), nu_vtx.end(), 0.0f);

    auto const &mct_h = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++)
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    const simb::MCTruth& mct = mct_h->at(0);
    if (!mct.NeutrinoSet()) return;

    auto const &neutrino = mct.GetNeutrino();
    auto const &nu = neutrino.Nu();
    common::True2RecoMappingXYZ(nu.T(), nu.Vx(), nu.Vy(), nu.Vz(), nu_vtx.data());

    if (neutrino.CCNC() != simb::kCC) return;
    const simb::MCParticle& lepton = neutrino.Lepton();

    if (abs(lepton.PdgCode()) != muon->PdgCode() || lepton.Momentum().Vect().Mag() < _MuonThreshold) return;

    for (const auto &t_part : *mcp_h) 
    {
        if (abs(t_part.PdgCode()) == muon->PdgCode() && t_part.Process() == "primary" && t_part.Momentum().Vect().Mag() > _MuonThreshold)
        {
            Particle muon_particle;
            muon_particle.pdg_code = t_part.PdgCode();
            muon_particle.track_id = t_part.TrackId();
            muon_particle.momentum = {static_cast<float>(t_part.Px()), static_cast<float>(t_part.Py()), static_cast<float>(t_part.Pz())};
            particle_elements.push_back(muon_particle);
        }

        if (abs(t_part.PdgCode()) == neutral_kaon->PdgCode() && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1 && !found_signature) 
        {
            std::vector<art::Ptr<simb::MCParticle>> dtrs = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (dtrs.size() != 1) continue; 

            auto g_part = dtrs.at(0);
            if (g_part->PdgCode() == kaon_short->PdgCode() && g_part->Process() == "Decay" && g_part->EndProcess() == "Decay" && g_part->NumberDaughters() == 2 && !found_signature)
            {
                auto daughters = common::GetDaughters(mcp_map.at(g_part->TrackId()), mcp_map);
                if (daughters.size() != 2) continue;

                std::vector<int> exp_dtrs = {-211, 211};
                std::vector<int> fnd_dtrs;

                for (const auto &dtr : daughters) 
                    fnd_dtrs.push_back(dtr->PdgCode());
                    
                std::sort(exp_dtrs.begin(), exp_dtrs.end());
                std::sort(fnd_dtrs.begin(), fnd_dtrs.end());
                if (fnd_dtrs == exp_dtrs) 
                {
                    for (const auto &dtr : daughters) 
                    {
                        Particle particle;
                        particle.pdg_code = dtr->PdgCode();
                        particle.track_id = dtr->TrackId();
                        particle.momentum = {static_cast<float>(dtr->Px()), static_cast<float>(dtr->Py()), static_cast<float>(dtr->Pz())};
                        
                        particle_elements.push_back(particle);
                    }

                    found_signature = std::all_of(daughters.begin(), daughters.end(), [&](const auto& dtr) {
                        return dtr->Momentum().Vect().Mag() >= _PionThreshold;
                    });

                    if (found_signature) break;
                }
            }
        }
    }
}

DEFINE_ART_CLASS_TOOL(NeutralKaonSignature)
} 

#endif