#ifndef ANALYSIS_NEUTRALDECAY_CXX
#define ANALYSIS_NEUTRALDECAY_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "TVector3.h"

#include "../CommonFuncs/BacktrackingFuncs.h"
#include "../CommonFuncs/TrackShowerScoreFuncs.h"
#include "../CommonFuncs/SpaceChargeCorrections.h"
#include "../CommonFuncs/Scatters.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace analysis
{

class StrangenessAnalysis : public AnalysisToolBase
{

public:
  
    StrangenessAnalysis(const fhicl::ParameterSet &pset);

    ~StrangenessAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);

    void analyzeEvent(art::Event const &e, bool fData) override;

    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    void SaveTruth(art::Event const &e);

    void setBranches(TTree *_tree) override;

    void resetTTree(TTree *_tree) override;

private:

    art::InputTag fMCTproducer;
    art::InputTag fMCRproducer;
    art::InputTag fMCPproducer;
    art::InputTag fHproducer;
    art::InputTag fBacktrackTag;

    TParticlePDG *neutral_kaon = TDatabasePDG::Instance()->GetParticle(311);
    TParticlePDG *kaon_short = TDatabasePDG::Instance()->GetParticle(310);
    TParticlePDG *kaon_long = TDatabasePDG::Instance()->GetParticle(130);
    TParticlePDG *lambda = TDatabasePDG::Instance()->GetParticle(3122);
    TParticlePDG *sigma_plus = TDatabasePDG::Instance()->GetParticle(3222); 
    TParticlePDG *sigma_minus = TDatabasePDG::Instance()->GetParticle(3112);
    TParticlePDG *sigma_zero = TDatabasePDG::Instance()->GetParticle(3212);

    unsigned int _mc_piplus_n_elas;
    unsigned int _mc_piplus_n_inelas;
    unsigned int _mc_piminus_n_elas;
    unsigned int _mc_piminus_n_inelas;

    int _mc_kshrt_piplus_tid; 
    int _mc_kshrt_piminus_tid;

    int _mc_kshrt_piplus_pdg;
    float _mc_kshrt_piplus_energy;
    float _mc_kshrt_piplus_px, _mc_kshrt_piplus_py, _mc_kshrt_piplus_pz;

    int _mc_kshrt_piminus_pdg;
    float _mc_kshrt_piminus_energy;
    float _mc_kshrt_piminus_px, _mc_kshrt_piminus_py, _mc_kshrt_piminus_pz;

    float _mc_kshrt_total_energy;

    float _mc_neutrino_vertex_x, _mc_neutrino_vertex_y, _mc_neutrino_vertex_z;
    float _mc_kaon_decay_x, _mc_kaon_decay_y, _mc_kaon_decay_z;
    float _mc_kaon_decay_distance;

    float _mc_piplus_impact_param;
    float _mc_piminus_impact_param;

    std::string _mc_piplus_endprocess;
    std::string _mc_piminus_endprocess;

    bool _mc_is_kshort_decay_pionic;

    bool _mc_has_lambda;
    bool _mc_has_sigma_minus;
    bool _mc_has_sigma_plus;
    bool _mc_has_sigma_zero;
};

StrangenessAnalysis::StrangenessAnalysis(const fhicl::ParameterSet &pset)
{
    fMCTproducer = pset.get<art::InputTag>("MCTproducer", "");
    fMCRproducer = pset.get<art::InputTag>("MCRproducer", "");
    fMCPproducer = pset.get<art::InputTag>("MCPproducer", "");
    fHproducer = pset.get<art::InputTag>("Hproducer", "");
    fBacktrackTag = pset.get<art::InputTag>("BacktrackTag", ""); 
}

void StrangenessAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

 void StrangenessAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
    if (fData) return;

    std::cout << "[StrangenessAnalysis] Analysing event..." << std::endl;
  
    // Load generator truth 
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);

    // Load transportation truth
    auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCRproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++)
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    auto mct = mct_h->at(0);
    auto neutrino = mct.GetNeutrino();
    auto nu = neutrino.Nu();

    _mc_neutrino_vertex_x = nu.Vx();
    _mc_neutrino_vertex_y = nu.Vy();
    _mc_neutrino_vertex_z = nu.Vz();

    TVector3 neutrino_vertex(_mc_neutrino_vertex_x, _mc_neutrino_vertex_y, _mc_neutrino_vertex_z);

    _mc_is_kshort_decay_pionic = false;

    size_t mct_n_part = mct.NParticles();
    for (size_t i = 0; i < mct_n_part; i++) 
    {
        auto const &t_part = mct.GetParticle(i);

        if (abs(t_part.PdgCode()) == lambda->PdgCode() && t_part.Process() == "primary")
            _mc_has_lambda = true;

        if (abs(t_part.PdgCode()) == sigma_plus->PdgCode() && t_part.Process() == "primary")
            _mc_has_sigma_plus = true;

        if (abs(t_part.PdgCode()) == sigma_zero->PdgCode() && t_part.Process() == "primary")
            _mc_has_sigma_zero = true;

        if (abs(t_part.PdgCode()) == sigma_minus->PdgCode() && t_part.Process() == "primary")
            _mc_has_sigma_minus = true;

        // Look for K0 at the generator level
        if (abs(t_part.PdgCode()) == neutral_kaon->PdgCode() && t_part.Process() == "primary" && !_mc_is_kshort_decay_pionic) 
        {
            // Find the corresponding K0-short or K0-long from GEANT4
            for (size_t j = 0; j < mcp_h->size(); j++) 
            {
                auto const &g_part = mcp_h->at(j);

                if (g_part.PdgCode() == kaon_short->PdgCode() && g_part.EndProcess() == "Decay" && !_mc_is_kshort_decay_pionic) 
                {
                    art::Ptr<simb::MCParticle> kaon_ptr(mcp_h, j);

                    auto daughters = common::GetDaughters(kaon_ptr, mcp_map);
                    if (daughters.size() == 2) 
                    {
                        std::vector<int> exp_dtrs = {-211, 211};
                        std::vector<int> fnd_dtrs;

                        for (const auto &dtr : daughters) 
                        {
                            fnd_dtrs.push_back(dtr->PdgCode());
                        }

                        std::sort(exp_dtrs.begin(), exp_dtrs.end());
                        std::sort(fnd_dtrs.begin(), fnd_dtrs.end());

                        if (fnd_dtrs == exp_dtrs) 
                        {
                            _mc_kshrt_total_energy = kaon_ptr->E();

                            _mc_kaon_decay_x = kaon_ptr->EndX();
                            _mc_kaon_decay_y = kaon_ptr->EndY();
                            _mc_kaon_decay_z = kaon_ptr->EndZ();

                            TVector3 kaon_decay(_mc_kaon_decay_x, _mc_kaon_decay_y, _mc_kaon_decay_z);
                            float decay_length = (kaon_decay - neutrino_vertex).Mag();
                            float phi_h = std::atan2(kaon_decay.Y(), kaon_decay.X());

                            _mc_kaon_decay_distance = decay_length;

                            for (const auto &dtr : daughters) 
                            {
                                std::cout << dtr->PdgCode() << std::endl;
                                TVector3 pion_mom(dtr->Px(), dtr->Py(), dtr->Pz());
                                float phi_i = std::atan2(pion_mom.Y(), pion_mom.X());
                                float d_0 = decay_length * std::sin(phi_i - phi_h);

                                unsigned int n_elas = 0;
                                unsigned int n_inelas = 0;
    
                                art::Ptr<simb::MCParticle> scat_part;
                                std::string scat_end_process;

                                common::GetNScatters(mcp_h, dtr, scat_part, n_elas, n_inelas);
                                scat_end_process = common::GetEndState(dtr, mcp_h);
                                std::cout << scat_end_process << std::endl;

                                if (dtr->PdgCode() == 211) // pion-plus
                                { 
                                    _mc_kshrt_piplus_tid = dtr->TrackId();
                                    _mc_kshrt_piplus_pdg = dtr->PdgCode();
                                    _mc_kshrt_piplus_energy = dtr->E();
                                    _mc_kshrt_piplus_px = dtr->Px();
                                    _mc_kshrt_piplus_py = dtr->Py();
                                    _mc_kshrt_piplus_pz = dtr->Pz();
                                    _mc_piplus_impact_param = d_0;
                                    _mc_piplus_n_elas = n_elas;
                                    _mc_piplus_n_inelas = n_inelas;
                                    _mc_piplus_endprocess = scat_end_process;
                                } 
                                else if (dtr->PdgCode() == -211) // pion-minus
                                { 
                                    _mc_kshrt_piminus_tid = dtr->TrackId();
                                    _mc_kshrt_piminus_pdg = dtr->PdgCode();
                                    _mc_kshrt_piminus_energy = dtr->E();
                                    _mc_kshrt_piminus_px = dtr->Px();
                                    _mc_kshrt_piminus_py = dtr->Py();
                                    _mc_kshrt_piminus_pz = dtr->Pz();
                                    _mc_piminus_impact_param = d_0;
                                    _mc_piminus_n_elas = n_elas;
                                    _mc_piminus_n_inelas = n_inelas;
                                    _mc_piminus_endprocess = scat_end_process;
                                }
                            }

                            _mc_is_kshort_decay_pionic = true;
                        }
                    }
                }
            }
        }

        if (_mc_is_kshort_decay_pionic) break;  
    }
}

void StrangenessAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
    /*common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer,
                                                                                            proxy::withAssociated<recob::Hit>(fCLSproducer));
    // somehow proxies don't work for the slice-hit association, so go back to old assns
    art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));*/

    // Build larpandora info:
    /*lar_pandora::LArPandoraHelper larpandora;
    lar_pandora::PFParticleVector pfparticles;
    lar_pandora::PFParticleMap particleMap;
    larpandora.CollectPFParticles(e, "pandora", pfparticles);
    larpandora.BuildPFParticleMap(pfparticles, particleMap);

    // load backtrack information
    std::vector<common::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!fData)
    {
        const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
        const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
        btparts_v = common::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
    }*/

    /*size_t pfpidx = 0;
    _n_pfps = 0;
    for (auto pfp : slice_pfp_v)
    {
        if (pfp->IsPrimary()) 
            continue;
        
        auto trkshrscore = common::GetTrackShowerScore(pfp_pxy);
        if (trkshrscore > 0.5)
        {
	        unsigned int trk_hits = 0;
            for (const auto &trk : pfp_pxy.get<recob::Track>())
            {
                
            }
        }
    }*/

    return;
}


void StrangenessAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("mc_piplus_tid", &_mc_kshrt_piplus_tid, "mc_piplus_tid/i");
    _tree->Branch("mc_piminus_tid", &_mc_kshrt_piminus_tid, "mc_piminus_tid/i");

    _tree->Branch("mc_piplus_n_elas", &_mc_piplus_n_elas, "mc_piplus_n_elas/i");
    _tree->Branch("mc_piplus_n_inelas", &_mc_piplus_n_inelas, "mc_piplus_n_inelas/i");
    _tree->Branch("mc_piminus_n_elas", &_mc_piminus_n_elas, "mc_piminus_n_elas/i");
    _tree->Branch("mc_piminus_n_inelas", &_mc_piminus_n_inelas, "mc_piminus_n_inelas/i");

    _tree->Branch("mc_kshrt_piplus_pdg", &_mc_kshrt_piplus_pdg, "mc_kshrt_piplus_pdg/I");
    _tree->Branch("mc_kshrt_piplus_energy", &_mc_kshrt_piplus_energy, "mc_kshrt_piplus_energy/F");
    _tree->Branch("mc_kshrt_piplus_px", &_mc_kshrt_piplus_px, "mc_kshrt_piplus_px/F");
    _tree->Branch("mc_kshrt_piplus_py", &_mc_kshrt_piplus_py, "mc_kshrt_piplus_py/F");
    _tree->Branch("mc_kshrt_piplus_pz", &_mc_kshrt_piplus_pz, "mc_kshrt_piplus_pz/F");

    _tree->Branch("mc_kshrt_piminus_pdg", &_mc_kshrt_piminus_pdg, "mc_kshrt_piminus_pdg/I");
    _tree->Branch("mc_kshrt_piminus_energy", &_mc_kshrt_piminus_energy, "mc_kshrt_piminus_energy/F");
    _tree->Branch("mc_kshrt_piminus_px", &_mc_kshrt_piminus_px, "mc_kshrt_piminus_px/F");
    _tree->Branch("mc_kshrt_piminus_py", &_mc_kshrt_piminus_py, "mc_kshrt_piminus_py/F");
    _tree->Branch("mc_kshrt_piminus_pz", &_mc_kshrt_piminus_pz, "mc_kshrt_piminus_pz/F");

    _tree->Branch("mc_kshrt_total_energy", &_mc_kshrt_total_energy, "mc_kshrt_total_energy/F");

    _tree->Branch("mc_neutrino_vertex_x", &_mc_neutrino_vertex_x, "mc_neutrino_vertex_x/F");
    _tree->Branch("mc_neutrino_vertex_y", &_mc_neutrino_vertex_y, "mc_neutrino_vertex_y/F");
    _tree->Branch("mc_neutrino_vertex_z", &_mc_neutrino_vertex_z, "mc_neutrino_vertex_z/F");

    _tree->Branch("mc_kaon_decay_x", &_mc_kaon_decay_x, "mc_kaon_decay_x/F");
    _tree->Branch("mc_kaon_decay_y", &_mc_kaon_decay_y, "mc_kaon_decay_y/F");
    _tree->Branch("mc_kaon_decay_z", &_mc_kaon_decay_z, "mc_kaon_decay_z/F");
    _tree->Branch("mc_kaon_decay_distance", &_mc_kaon_decay_distance, "mc_kaon_decay_distance/F");

    _tree->Branch("mc_piplus_impact_param", &_mc_piplus_impact_param, "mc_piplus_impact_param/F");
    _tree->Branch("mc_piminus_impact_param", &_mc_piminus_impact_param, "mc_piminus_impact_param/F");

    _tree->Branch("mc_piplus_endprocess", &_mc_piplus_endprocess); 
    _tree->Branch("mc_piminus_endprocess", &_mc_piminus_endprocess); 

    _tree->Branch("mc_is_pionic_kaon", &_mc_is_kshort_decay_pionic);

    _tree->Branch("mc_has_lambda", &_mc_has_lambda);
    _tree->Branch("mc_has_sigma_plus", &_mc_has_sigma_plus);
    _tree->Branch("mc_has_sigma_minus", &_mc_has_sigma_minus);
    _tree->Branch("mc_has_sigma_zero", &_mc_has_sigma_zero);
}

void StrangenessAnalysis::resetTTree(TTree *_tree)
{
    _mc_piplus_n_elas = 0;
    _mc_piplus_n_inelas = 0;
    _mc_piminus_n_elas = 0;
    _mc_piminus_n_inelas = 0;

    _mc_kshrt_piplus_pdg = 0;
    _mc_kshrt_piplus_energy = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_px = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_py = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_pz = std::numeric_limits<float>::lowest();

    _mc_kshrt_piminus_pdg = 0;
    _mc_kshrt_piminus_energy = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_px = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_py = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_pz = std::numeric_limits<float>::lowest();

    _mc_kshrt_total_energy = std::numeric_limits<float>::lowest();

    _mc_neutrino_vertex_x = std::numeric_limits<float>::lowest();
    _mc_neutrino_vertex_y = std::numeric_limits<float>::lowest();
    _mc_neutrino_vertex_z = std::numeric_limits<float>::lowest();

    _mc_kaon_decay_x = std::numeric_limits<float>::lowest();
    _mc_kaon_decay_y = std::numeric_limits<float>::lowest();
    _mc_kaon_decay_z = std::numeric_limits<float>::lowest();
    _mc_kaon_decay_distance = std::numeric_limits<float>::lowest();

    _mc_piplus_impact_param = std::numeric_limits<float>::lowest();
    _mc_piminus_impact_param = std::numeric_limits<float>::lowest();

    _mc_piplus_endprocess = ""; 
    _mc_piminus_endprocess = "";

    _mc_kshrt_piminus_tid = -1;
    _mc_kshrt_piplus_tid = -1;

    _mc_is_kshort_decay_pionic = false;

    _mc_has_lambda = false;
    _mc_has_sigma_plus = false;
    _mc_has_sigma_minus = false;
    _mc_has_sigma_zero = false;
}

DEFINE_ART_CLASS_TOOL(StrangenessAnalysis)
} 

#endif