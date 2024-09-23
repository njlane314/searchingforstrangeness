#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "AnalysisTools/AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "TVector3.h"

#include "CommonFunctions/Backtracking.h"
#include "CommonFunctions/Scores.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Geometry.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

class SignalTruthFilter : public art::EDFilter
{
public:
    explicit SignalTruthFilter(fhicl::ParameterSet const &p);

    SignalTruthFilter(SignalTruthFilter const &) = delete;
    SignalTruthFilter(SignalTruthFilter &&) = delete;
    SignalTruthFilter &operator=(SignalTruthFilter const &) = delete;
    SignalTruthFilter &operator=(SignalTruthFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _MCTproducer;
    art::InputTag _MCPproducer;

    float _PionThreshold;
    float _MuonThreshold;

    TParticlePDG *neutral_kaon = TDatabasePDG::Instance()->GetParticle(311);
    TParticlePDG *kaon_short = TDatabasePDG::Instance()->GetParticle(310);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
};

SignalTruthFilter::SignalTruthFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
{
    _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
    _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
    _PionThreshold = pset.get<float>("PionThreshold", 0.1); 
    _MuonThreshold = pset.get<float>("MuonThreshold", 0.1);
}

bool SignalTruthFilter::filter(art::Event &e)
{
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++) 
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    auto neutrino = mct_h->at(0).GetNeutrino();
    if (neutrino.CCNC() != simb::kCC) return false;
    auto lepton = neutrino.Lepton();
    if (abs(lepton.PdgCode()) != muon->PdgCode() ||  TVector3(lepton.Px(), lepton.Py(), lepton.Pz()).Mag() < _MuonThreshold)
        return false;

    bool found_signature = false;
    for (const auto &t_part : *mcp_h) 
    {
        if (abs(t_part.PdgCode()) != neutral_kaon->PdgCode() || t_part.Process() != "primary" ||
            t_part.EndProcess() != "Decay" || t_part.NumberDaughters() != 1) 
            continue;

        auto daughters = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
        if (daughters.size() != 1) continue;

        auto g_part = daughters.at(0);
        if (g_part->PdgCode() != kaon_short->PdgCode() || g_part->Process() != "Decay" ||
            g_part->EndProcess() != "Decay" || g_part->NumberDaughters() != 2) 
            continue;

        auto decay_daughters = common::GetDaughters(mcp_map.at(g_part->TrackId()), mcp_map);
        if (decay_daughters.size() != 2) continue;

        std::vector<int> expected_pdg = {-211, 211};  
        std::vector<int> found_pdg;

        for (const auto &dtr : decay_daughters) 
        {
            found_pdg.push_back(dtr->PdgCode());
        }

        std::sort(expected_pdg.begin(), expected_pdg.end());
        std::sort(found_pdg.begin(), found_pdg.end());

        if (found_pdg != expected_pdg) continue; 

        found_signature = std::all_of(decay_daughters.begin(), decay_daughters.end(), [&](const auto &dtr) 
        {
            return TVector3(dtr->Px(), dtr->Py(), dtr->Pz()).Mag() >= _PionThreshold;
        });

        if (found_signature) break;
    }

    return found_signature;
}

DEFINE_ART_MODULE(SignalTruthFilter)