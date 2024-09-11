#ifndef ANALYSIS_SIGNAL_CXX
#define ANALYSIS_SIGNAL_CXX

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

class SignalAnalysis : public AnalysisToolBase
{

public:
  
    SignalAnalysis(const fhicl::ParameterSet &pset);

    ~SignalAnalysis(){};

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

    bool _mc_has_muon;
    int _mc_muon_tid;
    int _mc_muon_pdg;
    float _mc_muon_energy;
    float _mc_muon_px; 
    float _mc_muon_py;
    float _mc_muon_pz;
    float _mc_muon_startx;
    float _mc_muon_starty;
    float _mc_muon_startz;
    float _mc_muon_endx;
    float _mc_muon_endy;
    float _mc_muon_endz;

    unsigned int _mc_piplus_n_elas;
    unsigned int _mc_piplus_n_inelas;
    unsigned int _mc_piminus_n_elas;
    unsigned int _mc_piminus_n_inelas;

    int _mc_kshrt_piplus_tid; 
    int _mc_kshrt_piminus_tid;

    int _mc_kshrt_piplus_pdg;
    float _mc_kshrt_piplus_energy;
    float _mc_kshrt_piplus_px, _mc_kshrt_piplus_py, _mc_kshrt_piplus_pz;
    float _mc_kshrt_piplus_startx, _mc_kshrt_piplus_starty, _mc_kshrt_piplus_startz;
    float _mc_kshrt_piplus_endx, _mc_kshrt_piplus_endy, _mc_kshrt_piplus_endz;

    int _mc_kshrt_piminus_pdg;
    float _mc_kshrt_piminus_energy;
    float _mc_kshrt_piminus_px, _mc_kshrt_piminus_py, _mc_kshrt_piminus_pz;
    float _mc_kshrt_piminus_startx, _mc_kshrt_piminus_starty, _mc_kshrt_piminus_startz;
    float _mc_kshrt_piminus_endx, _mc_kshrt_piminus_endy, _mc_kshrt_piminus_endz;

    float _mc_kshrt_total_energy;

    float _mc_neutrino_vertex_x, _mc_neutrino_vertex_y, _mc_neutrino_vertex_z;
    float _mc_kaon_decay_x, _mc_kaon_decay_y, _mc_kaon_decay_z;
    float _mc_kaon_decay_distance;

    float _mc_piplus_impact_param;
    float _mc_piminus_impact_param;

    float _mc_piplus_phi;
    float _mc_piminus_phi;

    std::vector<float> _pfp_muon_purity_v;
    std::vector<float> _pfp_muon_completeness_v;
    std::vector<float> _pfp_piplus_purity_v;
    std::vector<float> _pfp_piplus_completeness_v;
    std::vector<float> _pfp_piminus_purity_v;
    std::vector<float> _pfp_piminus_completeness_v;

    std::string _mc_piplus_endprocess;
    std::string _mc_piminus_endprocess;

    bool _mc_is_kshort_decay_pionic;

    bool _mc_has_lambda;
    bool _mc_has_sigma_minus;
    bool _mc_has_sigma_plus;
    bool _mc_has_sigma_zero;

    std::vector<float> _pfp_trk_sep_v;
    std::vector<float> _pfp_trk_phi_v;
    std::vector<float> _pfp_trk_d_v;
};

SignalAnalysis::SignalAnalysis(const fhicl::ParameterSet &pset)
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

void SignalAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void SignalAnalysis::analyzeEvent(art::Event const &e, bool is_data)
{
    if (is_data) 
        return;
  
    // Load generator truth 
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);

    // Load transportation truth
    auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

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
    for (size_t i = 0; i < mcp_h->size(); i++)
    {
        auto const &t_part = mcp_h->at(i);

        if (abs(t_part.PdgCode()) == lambda->PdgCode() && t_part.Process() == "primary")
            _mc_has_lambda = true;

        if (abs(t_part.PdgCode()) == sigma_plus->PdgCode() && t_part.Process() == "primary")
            _mc_has_sigma_plus = true;

        if (abs(t_part.PdgCode()) == sigma_zero->PdgCode() && t_part.Process() == "primary")
            _mc_has_sigma_zero = true;

        if (abs(t_part.PdgCode()) == sigma_minus->PdgCode() && t_part.Process() == "primary")
            _mc_has_sigma_minus = true;

        // Look for primary muon at the generator level
        if (abs(t_part.PdgCode()) == muon->PdgCode() && t_part.Process() == "primary")
        {
            _mc_has_muon = true;
            _mc_muon_tid = t_part.TrackId();
            _mc_muon_pdg = t_part.PdgCode();
            _mc_muon_energy = t_part.E();
            _mc_muon_px = t_part.Px(); 
            _mc_muon_py = t_part.Py();
            _mc_muon_pz = t_part.Pz();
            _mc_muon_startx = t_part.Vx();
            _mc_muon_starty = t_part.Vy();
            _mc_muon_startz = t_part.Vz();
            _mc_muon_endx = t_part.EndX();
            _mc_muon_endy = t_part.EndY();
            _mc_muon_endz = t_part.EndZ();
        }

        // Look for K0 at the generator level
        if (abs(t_part.PdgCode()) == neutral_kaon->PdgCode() && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1 && !_mc_is_kshort_decay_pionic) 
        {
            std::vector<art::Ptr<simb::MCParticle>> dtrs = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (dtrs.size() != 1) 
                continue; 

            auto g_part = dtrs.at(0);
            if (g_part->PdgCode() == kaon_short->PdgCode() && g_part->Process() == "Decay" && g_part->EndProcess() == "Decay" && g_part->NumberDaughters() == 2 && !_mc_is_kshort_decay_pionic)
            {
                std::cout << "[SignalAnalysis::analyzeEvent] Found kaon-short..." << std::endl;
                auto daughters = common::GetDaughters(mcp_map.at(g_part->TrackId()), mcp_map);
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
                        TVector3 kaon_decay(g_part->EndX(), g_part->EndY(), g_part->EndZ());
                        float decay_length = (kaon_decay - neutrino_vertex).Mag();

                        _mc_kshrt_total_energy = g_part->E();

                        _mc_kaon_decay_x = g_part->EndX();
                        _mc_kaon_decay_y = g_part->EndY();
                        _mc_kaon_decay_z = g_part->EndZ();
                        
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

                            if (dtr->PdgCode() == 211) // pion-plus
                            { 
                                std::cout << "[SignalAnalysis::analyzeEvent] Found decay pion-minus with track identifier " << dtr->TrackId() << "..." << std::endl;
                                _mc_kshrt_piplus_tid = dtr->TrackId();
                                _mc_kshrt_piplus_pdg = dtr->PdgCode();
                                _mc_kshrt_piplus_energy = dtr->E();
                                _mc_kshrt_piplus_px = dtr->Px();
                                _mc_kshrt_piplus_py = dtr->Py();
                                _mc_kshrt_piplus_pz = dtr->Pz();
                                _mc_kshrt_piplus_startx = dtr->Vx();
                                _mc_kshrt_piplus_starty = dtr->Vy();
                                _mc_kshrt_piplus_startz = dtr->Vz();
                                _mc_kshrt_piplus_endx = dtr->EndX();
                                _mc_kshrt_piplus_endy = dtr->EndY();
                                _mc_kshrt_piplus_endz = dtr->EndZ();
                                _mc_piplus_phi = phi_i;
                                _mc_piplus_impact_param = d_0;
                                _mc_piplus_n_elas = n_elas;
                                _mc_piplus_n_inelas = n_inelas;
                                _mc_piplus_endprocess = scat_end_process;
                            } 
                            else if (dtr->PdgCode() == -211) // pion-minus
                            { 
                                std::cout << "[SignalAnalysis::analyzeEvent] Found decay pion-minus with track identifier " << dtr->TrackId() << "..." << std::endl;
                                _mc_kshrt_piminus_tid = dtr->TrackId();
                                _mc_kshrt_piminus_pdg = dtr->PdgCode();
                                _mc_kshrt_piminus_energy = dtr->E();
                                _mc_kshrt_piminus_px = dtr->Px();
                                _mc_kshrt_piminus_py = dtr->Py();
                                _mc_kshrt_piminus_pz = dtr->Pz();
                                _mc_kshrt_piminus_startx = dtr->Vx();
                                _mc_kshrt_piminus_starty = dtr->Vy();
                                _mc_kshrt_piminus_startz = dtr->Vz();
                                _mc_kshrt_piminus_endx = dtr->EndX();
                                _mc_kshrt_piminus_endy = dtr->EndY();
                                _mc_kshrt_piminus_endz = dtr->EndZ();
                                _mc_piminus_phi = phi_i;
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

        if (_mc_is_kshort_decay_pionic)
            break;
    }
}

void SignalAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected)
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

    if (_mc_is_kshort_decay_pionic)
    {
        std::cout << slice_pfp_v.size() << std::endl;
        for (const common::ProxyPfpElem_t &pfp_pxy : slice_pfp_v)
        {
            double reco_nu_vtx_sce_x; 
            double reco_nu_vtx_sce_y;
            double reco_nu_vtx_sce_z;

            if (pfp_pxy->IsPrimary()) 
            {
                double xyz[3] = {};

                auto vtx = pfp_pxy.get<recob::Vertex>();
                if (vtx.size() == 1)
                {
                    vtx.at(0)->XYZ(xyz);
                    auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

                    double reco_nu_vtx_x = nuvtx.X();
                    double reco_nu_vtx_y = nuvtx.Y();
                    double reco_nu_vtx_z = nuvtx.Z();

                    float reco_nu_vtx_sce[3];
                    common::ApplySCECorrectionXYZ(reco_nu_vtx_x, reco_nu_vtx_y, reco_nu_vtx_z, reco_nu_vtx_sce);
                    reco_nu_vtx_sce_x = reco_nu_vtx_sce[0];
                    reco_nu_vtx_sce_y = reco_nu_vtx_sce[1];
                    reco_nu_vtx_sce_z = reco_nu_vtx_sce[2];
                }
            }
            
            auto trk_v = pfp_pxy.get<recob::Track>();

            if (trk_v.size() == 1) 
            {
                auto trk = trk_v.at(0);

                auto trk_strt = trk->Start();
                float trk_strt_sce[3];
                common::ApplySCECorrectionXYZ(trk_strt.X(), trk_strt.Y(), trk_strt.Z(), trk_strt_sce);
                float phi_h = std::atan2(trk_strt_sce[0], trk_strt_sce[1]);

                auto trk_end = trk->End();
                float trk_end_sce[3];
                common::ApplySCECorrectionXYZ(trk_end.X(), trk_end.Y(), trk_end.Z(), trk_end_sce);

                float trk_sep = common::distance3d(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z, trk_end_sce[0], trk_end_sce[1], trk_end_sce[2]);
                float trk_phi = std::atan2(trk_end_sce[0], trk_end_sce[1]);
                float trk_d = trk_sep * sin(trk_phi - phi_h);

                _pfp_trk_sep_v.push_back(trk_sep); 
                _pfp_trk_phi_v.push_back(trk_phi);
                _pfp_trk_d_v.push_back(trk_d);
            }
            else 
            {
                _pfp_trk_sep_v.push_back(std::numeric_limits<float>::lowest()); 
                _pfp_trk_phi_v.push_back(std::numeric_limits<float>::lowest());
                _pfp_trk_d_v.push_back(std::numeric_limits<float>::lowest());
            }
        }

        art::ValidHandle<std::vector<recob::Hit>> in_hits = e.getValidHandle<std::vector<recob::Hit>>(_Hproducer);
        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> mcp_bkth_assoc = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(in_hits, e, _BacktrackTag));

        int mcp_mu_hits = 0;
        int mcp_piplus_hits = 0;
        int mcp_piminus_hits = 0;

        for (unsigned int ih = 0; ih < in_hits->size(); ih++)
        {
            auto assmcp = mcp_bkth_assoc->at(ih);
            auto assmdt = mcp_bkth_assoc->data(ih);
            for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
            {
                auto mcp = assmcp[ia];
                auto amd = assmdt[ia];
                if (amd->isMaxIDE != 1)
                    continue;

                if (mcp->TrackId() == _mc_muon_tid)
                    mcp_mu_hits++;
                else if (mcp->TrackId() == _mc_kshrt_piplus_tid)
                    mcp_piplus_hits++;
                else if (mcp->TrackId() == _mc_kshrt_piminus_tid)
                    mcp_piminus_hits++;
            }
        }

        for (const common::ProxyPfpElem_t &pfp_pxy : slice_pfp_v)
        {
            std::vector<art::Ptr<recob::Hit>> pfp_hits;
            auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();

            for (auto ass_clus : clus_pxy_v)
            {
                const auto &clus = clus_proxy[ass_clus.key()];
                auto clus_hit_v = clus.get<recob::Hit>();

                for (const auto &hit : clus_hit_v)
                {
                    pfp_hits.push_back(hit);
                }
            } 

            std::vector<simb::MCParticle> assmcp;
            std::vector<anab::BackTrackerHitMatchingData> amd;

            int pfp_mu_hits = 0; 
            int pfp_piplus_hits = 0;
            int pfp_piminus_hits = 0;

            for (auto hit : pfp_hits)
            {
                auto assmcp = mcp_bkth_assoc->at(hit.key());       
                auto assmdt = mcp_bkth_assoc->data(hit.key());    
                for (size_t i = 0; i < assmcp.size(); i++)
                {
                    if (assmdt[i]->isMaxIDE != 1) 
                        continue;

                    const int track_id = assmcp[i]->TrackId();
                    if (track_id == _mc_muon_tid)
                        pfp_mu_hits++;
                    else if (track_id == _mc_kshrt_piplus_tid)
                        pfp_piplus_hits++;
                    else if (track_id == _mc_kshrt_piminus_tid)
                        pfp_piminus_hits++;
                }
            }

            float muon_purity = (pfp_hits.size() > 0) ? static_cast<float>(pfp_mu_hits) / pfp_hits.size() : 0.0f;
            float muon_completeness = (mcp_mu_hits > 0) ? static_cast<float>(pfp_mu_hits) / mcp_mu_hits : 0.0f;

            float piplus_purity = (pfp_hits.size() > 0) ? static_cast<float>(pfp_piplus_hits) / pfp_hits.size() : 0.0f;
            float piplus_completeness = (mcp_piplus_hits > 0) ? static_cast<float>(pfp_piplus_hits) / mcp_piplus_hits : 0.0f;

            float piminus_purity = (pfp_hits.size() > 0) ? static_cast<float>(pfp_piminus_hits) / pfp_hits.size() : 0.0f;
            float piminus_completeness = (mcp_piminus_hits > 0) ? static_cast<float>(pfp_piminus_hits) / mcp_piminus_hits : 0.0f;

            _pfp_muon_purity_v.push_back(muon_purity);
            _pfp_muon_completeness_v.push_back(muon_completeness);
            _pfp_piplus_purity_v.push_back(piplus_purity);
            _pfp_piplus_completeness_v.push_back(piplus_completeness);
            _pfp_piminus_purity_v.push_back(piminus_purity);
            _pfp_piminus_completeness_v.push_back(piminus_completeness);
            
            std::cout << "Muon Purity: " << muon_purity << ", Muon Completeness: " << muon_completeness << std::endl;
            std::cout << "Pion+ Purity: " << piplus_purity << ", Pion+ Completeness: " << piplus_completeness << std::endl;
            std::cout << "Pion- Purity: " << piminus_purity << ", Pion- Completeness: " << piminus_completeness << std::endl;
        }
    }

    return;
}

void SignalAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("mc_muon_tid", &_mc_muon_tid, "mc_muon_tid/I");
    _tree->Branch("mc_muon_pdg", &_mc_muon_pdg, "mc_muon_pdg/I");
    _tree->Branch("mc_muon_energy", &_mc_muon_energy, "mc_muon_energy/F");
    _tree->Branch("mc_muon_px", &_mc_muon_px, "mc_muon_px/F");
    _tree->Branch("mc_muon_py", &_mc_muon_py, "mc_muon_py/F");
    _tree->Branch("mc_muon_pz", &_mc_muon_pz, "mc_muon_pz/F");
    _tree->Branch("mc_muon_startx", &_mc_muon_startx, "mc_muon_startx/F");
    _tree->Branch("mc_muon_starty", &_mc_muon_starty, "mc_muon_starty/F");
    _tree->Branch("mc_muon_startz", &_mc_muon_startz, "mc_muon_startz/F");
    _tree->Branch("mc_muon_endx", &_mc_muon_endx, "mc_muon_endx/F");
    _tree->Branch("mc_muon_endy", &_mc_muon_endy, "mc_muon_endy/F");
    _tree->Branch("mc_muon_endz", &_mc_muon_endz, "mc_muon_endz/F");

    _tree->Branch("mc_piplus_tid", &_mc_kshrt_piplus_tid, "mc_piplus_tid/i");
    _tree->Branch("mc_piminus_tid", &_mc_kshrt_piminus_tid, "mc_piminus_tid/i");

    _tree->Branch("mc_kshrt_piplus_startx", &_mc_kshrt_piplus_startx, "mc_kshrt_piplus_startx/F");
    _tree->Branch("mc_kshrt_piplus_starty", &_mc_kshrt_piplus_starty, "mc_kshrt_piplus_starty/F");
    _tree->Branch("mc_kshrt_piplus_startz", &_mc_kshrt_piplus_startz, "mc_kshrt_piplus_startz/F");
    _tree->Branch("mc_kshrt_piminus_startx", &_mc_kshrt_piminus_startx, "mc_kshrt_piminus_startx/F");
    _tree->Branch("mc_kshrt_piminus_starty", &_mc_kshrt_piminus_starty, "mc_kshrt_piminus_starty/F");
    _tree->Branch("mc_kshrt_piminus_startz", &_mc_kshrt_piminus_startz, "mc_kshrt_piminus_startz/F");

    _tree->Branch("mc_kshrt_piplus_endx", &_mc_kshrt_piplus_endx, "mc_kshrt_piplus_endx/F");
    _tree->Branch("mc_kshrt_piplus_endy", &_mc_kshrt_piplus_endy, "mc_kshrt_piplus_endy/F");
    _tree->Branch("mc_kshrt_piplus_endz", &_mc_kshrt_piplus_endz, "mc_kshrt_piplus_endz/F");
    _tree->Branch("mc_kshrt_piminus_endx", &_mc_kshrt_piminus_endx, "mc_kshrt_piminus_endx/F");
    _tree->Branch("mc_kshrt_piminus_endy", &_mc_kshrt_piminus_endy, "mc_kshrt_piminus_endy/F");
    _tree->Branch("mc_kshrt_piminus_endz", &_mc_kshrt_piminus_endz, "mc_kshrt_piminus_endz/F");

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

    _tree->Branch("pfp_muon_purity", &_pfp_muon_purity_v);
    _tree->Branch("pfp_muon_completeness", &_pfp_muon_completeness_v); 
    _tree->Branch("pfp_piplus_purity", &_pfp_piplus_purity_v);
    _tree->Branch("pfp_piplus_completeness", &_pfp_piplus_completeness_v);
    _tree->Branch("pfp_piminus_purity", &_pfp_piminus_purity_v);
    _tree->Branch("pfp_piminus_completeness", &_pfp_piminus_completeness_v);

    _tree->Branch("mc_piplus_phi", &_mc_piplus_phi, "mc_piplus_phi/F");
    _tree->Branch("mc_piminus_phi", &_mc_piminus_phi, "mc_piminus_phi/F");

    _tree->Branch("mc_piplus_endprocess", &_mc_piplus_endprocess); 
    _tree->Branch("mc_piminus_endprocess", &_mc_piminus_endprocess); 

    _tree->Branch("mc_is_kshort_decay_pionic", &_mc_is_kshort_decay_pionic);

    _tree->Branch("mc_has_muon", &_mc_has_muon);
    _tree->Branch("mc_has_lambda", &_mc_has_lambda);
    _tree->Branch("mc_has_sigma_plus", &_mc_has_sigma_plus);
    _tree->Branch("mc_has_sigma_minus", &_mc_has_sigma_minus);
    _tree->Branch("mc_has_sigma_zero", &_mc_has_sigma_zero);

    _tree->Branch("all_pfp_trk_sep", &_pfp_trk_sep_v);
    _tree->Branch("all_pfp_trk_d", &_pfp_trk_d_v);
    _tree->Branch("all_pfp_trk_phi", &_pfp_trk_phi_v);
}

void SignalAnalysis::resetTTree(TTree *_tree)
{
    _mc_muon_tid = -1;
    _mc_muon_pdg = 0;
    _mc_muon_energy = std::numeric_limits<float>::lowest();
    _mc_muon_px = std::numeric_limits<float>::lowest();
    _mc_muon_py = std::numeric_limits<float>::lowest();
    _mc_muon_pz = std::numeric_limits<float>::lowest();
    _mc_muon_startx = std::numeric_limits<float>::lowest();
    _mc_muon_starty = std::numeric_limits<float>::lowest();
    _mc_muon_startz = std::numeric_limits<float>::lowest();
    _mc_muon_endx = std::numeric_limits<float>::lowest();
    _mc_muon_endy = std::numeric_limits<float>::lowest();
    _mc_muon_endz = std::numeric_limits<float>::lowest();

    _mc_piplus_n_elas = 0;
    _mc_piplus_n_inelas = 0;
    _mc_piminus_n_elas = 0;
    _mc_piminus_n_inelas = 0;

    _mc_kshrt_piplus_pdg = 0;
    _mc_kshrt_piplus_energy = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_px = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_py = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_pz = std::numeric_limits<float>::lowest();

    _mc_kshrt_piplus_startx = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_starty = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_startz = std::numeric_limits<float>::lowest();

    _mc_kshrt_piplus_endx = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_endy = std::numeric_limits<float>::lowest();
    _mc_kshrt_piplus_endz = std::numeric_limits<float>::lowest();

    _mc_kshrt_piminus_pdg = 0;
    _mc_kshrt_piminus_energy = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_px = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_py = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_pz = std::numeric_limits<float>::lowest();

    _mc_kshrt_piminus_startx = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_starty = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_startz = std::numeric_limits<float>::lowest();

    _mc_kshrt_piminus_endx = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_endy = std::numeric_limits<float>::lowest();
    _mc_kshrt_piminus_endz = std::numeric_limits<float>::lowest();

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

    _mc_piplus_phi = std::numeric_limits<float>::lowest();
    _mc_piminus_phi = std::numeric_limits<float>::lowest();

    _mc_piplus_endprocess = ""; 
    _mc_piminus_endprocess = "";

    _mc_kshrt_piminus_tid = -1;
    _mc_kshrt_piplus_tid = -1;

    _mc_is_kshort_decay_pionic = false;

    _mc_has_muon = false;
    _mc_has_lambda = false;
    _mc_has_sigma_plus = false;
    _mc_has_sigma_minus = false;
    _mc_has_sigma_zero = false;

    _pfp_trk_d_v.clear();
    _pfp_trk_phi_v.clear();
    _pfp_trk_sep_v.clear();

    _pfp_muon_purity_v.clear();
    _pfp_muon_completeness_v.clear();
    _pfp_piplus_purity_v.clear();
    _pfp_piplus_completeness_v.clear();
    _pfp_piminus_purity_v.clear();
    _pfp_piminus_completeness_v.clear();
}

DEFINE_ART_CLASS_TOOL(SignalAnalysis)
} 

#endif