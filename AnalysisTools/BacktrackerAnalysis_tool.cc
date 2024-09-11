#ifndef ANALYSIS_BACKTRACKER_CXX
#define ANALYSIS_BACKTRACKER_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "TVector3.h"

#include "../CommonFunctions/BacktrackingFuncs.h"
#include "../CommonFunctions/TrackShowerScoreFuncs.h"
#include "../CommonFunctions/SpaceChargeCorrections.h"
#include "../CommonFunctions/Scatters.h"
#include "../CommonFunctions/Geometry.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace analysis
{

class BacktrackerAnalysis : public AnalysisToolBase
{

public:
  
    BacktrackerAnalysis(const fhicl::ParameterSet &pset);

    ~BacktrackerAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);

    void analyzeEvent(art::Event const &e, bool is_data) override;

    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;

    void setBranches(TTree *_tree) override;

    void resetTTree(TTree *_tree) override;

private:

    art::InputTag _MCTproducer;
    art::InputTag _MCRproducer;
    art::InputTag _MCPproducer;
    art::InputTag _Hproducer;
    art::InputTag _BacktrackTag;
    art::InputTag _PFPproducer;
    art::InputTag _CLSproducer; 
    art::InputTag _SLCproducer;
    art::InputTag _TRKproducer;
    art::InputTag _VTXproducer;
    art::InputTag _PCAproducer;
    art::InputTag _SHRproducer;

    TParticlePDG *neutral_kaon = TDatabasePDG::Instance()->GetParticle(311);
    TParticlePDG *kaon_short = TDatabasePDG::Instance()->GetParticle(310);
    TParticlePDG *kaon_long = TDatabasePDG::Instance()->GetParticle(130);
    TParticlePDG *lambda = TDatabasePDG::Instance()->GetParticle(3122);
    TParticlePDG *sigma_plus = TDatabasePDG::Instance()->GetParticle(3222); 
    TParticlePDG *sigma_minus = TDatabasePDG::Instance()->GetParticle(3112);
    TParticlePDG *sigma_zero = TDatabasePDG::Instance()->GetParticle(3212);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
    TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);

    std::vector<int> _bt_pdg_v;
    std::vector<float> _bt_px_v;
    std::vector<float> _bt_py_v;
    std::vector<float> _bt_pz_v;
    std::vector<float> _bt_energy_v;
    std::vector<std::vector<unsigned int>> _bt_tids_v;
};

BacktrackerAnalysis::BacktrackerAnalysis(const fhicl::ParameterSet &pset)
{
    _PFPproducer = pset.get<art::InputTag>("PFPproducer");
    _CLSproducer = pset.get<art::InputTag>("CLSproducer");
    _SLCproducer = pset.get<art::InputTag>("SLCproducer");
    _TRKproducer = pset.get<art::InputTag>("TRKproducer");
    _VTXproducer = pset.get<art::InputTag>("VTXproducer");
    _PCAproducer = pset.get<art::InputTag>("PCAproducer");
    _SHRproducer = pset.get<art::InputTag>("SHRproducer");
    _MCTproducer = pset.get<art::InputTag>("MCTproducer", "");
    _MCRproducer = pset.get<art::InputTag>("MCRproducer", "");
    _MCPproducer = pset.get<art::InputTag>("MCPproducer", "");
    _Hproducer = pset.get<art::InputTag>("Hproducer", "");
    _BacktrackTag = pset.get<art::InputTag>("BacktrackTag", ""); 
}

void BacktrackerAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void BacktrackerAnalysis::analyzeEvent(art::Event const &e, bool is_data)
{
}

void BacktrackerAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected)
{
    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, _PFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(_PFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(_CLSproducer),
                                                        proxy::withAssociated<recob::Slice>(_SLCproducer),
                                                        proxy::withAssociated<recob::Track>(_TRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(_VTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(_PCAproducer),
                                                        proxy::withAssociated<recob::Shower>(_SHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(_PFPproducer));

    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, _CLSproducer,
                                                                                            proxy::withAssociated<recob::Hit>(_CLSproducer));

    std::vector<common::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

    if (!is_data)
    {
        const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(_MCRproducer));
        const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(_MCRproducer));
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(_Hproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, _BacktrackTag));
        btparts_v = common::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);

        std::cout << "[BacktrackerAnalysis::analyzeSlice] The size of the backtracker particle vector is " << btparts_v.size() << std::endl;
        for (const auto &btpart : btparts_v)
        {
            _bt_pdg_v.push_back(btpart.pdg);
            _bt_px_v.push_back(btpart.px);
            _bt_py_v.push_back(btpart.py);
            _bt_pz_v.push_back(btpart.pz);
            _bt_energy_v.push_back(btpart.e);
            _bt_tids_v.push_back(btpart.tids);
        }
    }
}

void BacktrackerAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("bt_pdg", &_bt_pdg_v);
    _tree->Branch("bt_px", &_bt_px_v);
    _tree->Branch("bt_py", &_bt_py_v);
    _tree->Branch("bt_pz", &_bt_pz_v);
    _tree->Branch("bt_energy", &_bt_energy_v);
    _tree->Branch("bt_tids", &_bt_tids_v);
}

void BacktrackerAnalysis::resetTTree(TTree *_tree)
{
    _bt_pdg_v.clear();
    _bt_px_v.clear();
    _bt_py_v.clear();
    _bt_pz_v.clear();
    _bt_energy_v.clear();
    _bt_tids_v.clear();
}

DEFINE_ART_CLASS_TOOL(BacktrackerAnalysis)
} 

#endif