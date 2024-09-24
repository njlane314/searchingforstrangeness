#ifndef ANALYSIS_EVENTCATEGORY
#define ANALYSIS_EVENTCATEGORY

#include <iostream>
#include "AnalysisToolBase.h"

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
#include "CommonFunctions/Pandora.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace analysis
{

struct Particle {
    int tid = -1;
    int pdg = 0;
    float energy = std::numeric_limits<float>::lowest();
    float px = std::numeric_limits<float>::lowest(), py = std::numeric_limits<float>::lowest(), pz = std::numeric_limits<float>::lowest();
    float startx = std::numeric_limits<float>::lowest(), starty = std::numeric_limits<float>::lowest(), startz = std::numeric_limits<float>::lowest();
    float endx = std::numeric_limits<float>::lowest(), endy = std::numeric_limits<float>::lowest(), endz = std::numeric_limits<float>::lowest();
    unsigned int n_elas = 0; 
    unsigned int n_inelas = 0;
    std::string endstate = "";
};

struct Neutrino {
    int nu_pdg = -1;         
    int ccnc = -1;    
    int interaction = -1;
    float nu_e = std::numeric_limits<float>::lowest();   
    float nu_theta = std::numeric_limits<float>::lowest();   
    float nu_pt = std::numeric_limits<float>::lowest();     
    int target_nucleus = -1;
    int hit_nucleon = -1;
    float W = std::numeric_limits<float>::lowest();
    float X = std::numeric_limits<float>::lowest();
    float Y = std::numeric_limits<float>::lowest();
    float QSqr = std::numeric_limits<float>::lowest();

    float px = std::numeric_limits<float>::lowest(), py = std::numeric_limits<float>::lowest(), pz = std::numeric_limits<float>::lowest();

    float true_vtx_sce_x = std::numeric_limits<float>::lowest(), true_vtx_sce_y = std::numeric_limits<float>::lowest(), true_vtx_sce_z = std::numeric_limits<float>::lowest();
    float true_vtx_sce_u_wire = std::numeric_limits<float>::lowest(), true_vtx_sce_w_wire = std::numeric_limits<float>::lowest(), true_vtx_sce_v_wire = std::numeric_limits<float>::lowest();
};

class EventCategoryAnalysis : public AnalysisToolBase
{

public:
    EventCategoryAnalysis(const fhicl::ParameterSet &pset);
    ~EventCategoryAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);
    void analyzeEvent(art::Event const &e, bool is_data) override;
    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

    void setParticleBranches(TTree *_tree, const std::string &prefix, Particle &particle);

    void fillNeutrino(Neutrino& p, const simb::MCNeutrino& neutrino, const simb::MCParticle& nu)
    {
        p.nu_pdg = nu.PdgCode();
        p.nu_e = nu.Trajectory().E(0);
        p.nu_theta = neutrino.Theta();
        p.nu_pt = neutrino.Pt();
        p.ccnc = neutrino.CCNC();
        p.interaction = neutrino.Mode();
        p.target_nucleus = neutrino.Target();
        p.hit_nucleon = neutrino.HitNuc();
        p.W = neutrino.W();
        p.X = neutrino.X();
        p.Y = neutrino.Y();
        p.QSqr = neutrino.QSqr();

        p.px = nu.Px();
        p.py = nu.Py();
        p.pz = nu.Pz();
    
        float true_nu_vtx_sce[3];
        common::True2RecoMappingXYZ(nu.T(), nu.Vx(), nu.Vy(), nu.Vz(), true_nu_vtx_sce);

        p.true_vtx_sce_x = true_nu_vtx_sce[0];
        p.true_vtx_sce_y = true_nu_vtx_sce[1];
        p.true_vtx_sce_z = true_nu_vtx_sce[2];

        p.true_vtx_sce_u_wire = (common::ProjectToWireView(p.true_vtx_sce_x, p.true_vtx_sce_y, p.true_vtx_sce_z, common::TPC_VIEW_U)).Z(); 
        p.true_vtx_sce_v_wire = (common::ProjectToWireView(p.true_vtx_sce_x, p.true_vtx_sce_y, p.true_vtx_sce_z, common::TPC_VIEW_V)).Z(); 
        p.true_vtx_sce_w_wire = (common::ProjectToWireView(p.true_vtx_sce_x, p.true_vtx_sce_y, p.true_vtx_sce_z, common::TPC_VIEW_W)).Z(); 
    }

    void fillParticle(const art::ValidHandle<std::vector<simb::MCParticle>>& mcp_h, const art::Ptr<simb::MCParticle>& particle, Particle& p)
    {
        p.tid = particle->TrackId();
        p.pdg = particle->PdgCode();
        p.energy = particle->E();
        p.px = particle->Px();
        p.py = particle->Py();
        p.pz = particle->Pz();
        p.startx = particle->Vx();
        p.starty = particle->Vy();
        p.startz = particle->Vz();
        p.endx = particle->EndX();
        p.endy = particle->EndY();
        p.endz = particle->EndZ();

        if (abs(particle->PdgCode()) == 211)
        {
            art::Ptr<simb::MCParticle> final_scatter_particle;
            common::GetNScatters(mcp_h, particle, final_scatter_particle, p.n_elas, p.n_inelas);

            p.endstate = common::GetEndState(particle, mcp_h);
        }
    }

private:
    art::InputTag _MCTproducer;
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

    float _PionThreshold;
    float _MuonThreshold;

    TParticlePDG *neutral_kaon = TDatabasePDG::Instance()->GetParticle(311);
    TParticlePDG *kaon_short = TDatabasePDG::Instance()->GetParticle(310);
    TParticlePDG *kaon_long = TDatabasePDG::Instance()->GetParticle(130);
    TParticlePDG *lambda = TDatabasePDG::Instance()->GetParticle(3122);
    TParticlePDG *sigma_plus = TDatabasePDG::Instance()->GetParticle(3222); 
    TParticlePDG *sigma_minus = TDatabasePDG::Instance()->GetParticle(3112);
    TParticlePDG *sigma_zero = TDatabasePDG::Instance()->GetParticle(3212);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
    TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);

    bool _pass_preselection;

    Neutrino _neutrino;

    Particle _mcp_mu;
    Particle _mcp_kshort;
    Particle _mcp_piplus;
    Particle _mcp_piminus;

    std::vector<int> _final_state;

    bool _found_signature;

    float _reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z;
    float _reco_nu_vtx_sce_u_wire, _reco_nu_vtx_sce_v_wire, _reco_nu_vtx_sce_w_wire;

    int _n_tracks;
    int _n_showers;

    std::vector<int> _mcp_pfp_id_v;
    std::vector<float> _mcp_mu_purity_v;
    std::vector<float> _mcp_mu_comple_v;
    std::vector<float> _mcp_piplus_purity_v;
    std::vector<float> _mcp_piplus_comple_v;
    std::vector<float> _mcp_piminus_purity_v;
    std::vector<float> _mcp_piminus_comple_v;
};

EventCategoryAnalysis::EventCategoryAnalysis(const fhicl::ParameterSet &pset)
{
    _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
    _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
    _Hproducer = pset.get<art::InputTag>("Hproducer", "gaushit");
    _BacktrackTag = pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
    _PFPproducer = pset.get<art::InputTag>("PFPproducer", "pandora");
    _CLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
    _SLCproducer = pset.get<art::InputTag>("SLCproducer", "pandora");
    _TRKproducer = pset.get<art::InputTag>("TRKproducer", "pandora");
    _VTXproducer = pset.get<art::InputTag>("VTXproducer", "pandora");
    _PCAproducer = pset.get<art::InputTag>("PCAproducer", "pandora");
    _SHRproducer = pset.get<art::InputTag>("SHRproducer", "pandora"); 

    _PionThreshold = pset.get<float>("PionThreshold", 0.1); 
    _MuonThreshold = pset.get<float>("MuonThreshold", 0.1);
}

void EventCategoryAnalysis::configure(fhicl::ParameterSet const &pset)
{
}

void EventCategoryAnalysis::analyzeEvent(art::Event const &e, bool is_data)
{
    if (is_data) 
        return;
  
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++)
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    auto mct = mct_h->at(0);
    _found_signature = false;

    if (mct.NeutrinoSet())
    {
        auto const &neutrino = mct.GetNeutrino();
        auto const &nu = neutrino.Nu();
        fillNeutrino(_neutrino, neutrino, nu);

        if (neutrino.CCNC() != simb::kCC) return;
        auto lepton = neutrino.Lepton();
        if (abs(lepton.PdgCode()) != muon->PdgCode() ||  TVector3(lepton.Px(), lepton.Py(), lepton.Pz()).Mag() < _MuonThreshold)
            return;
        else
        {
            art::Ptr<simb::MCParticle> lepton_ptr = mcp_map[lepton.TrackId()];
            fillParticle(mcp_h, lepton_ptr, _mcp_mu);
        }

        for (const auto &t_part : *mcp_h)
        {
            if (t_part.Process() != "primary") 
                continue;

            _final_state.push_back(t_part.PdgCode());
        }

        for (const auto &t_part : *mcp_h) 
        {
            if (abs(t_part.PdgCode()) == neutral_kaon->PdgCode() && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1 && !_found_signature) 
            {
                std::vector<art::Ptr<simb::MCParticle>> dtrs = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
                if (dtrs.size() != 1) 
                    continue; 

                auto g_part = dtrs.at(0);
                if (g_part->PdgCode() == kaon_short->PdgCode() && g_part->Process() == "Decay" && g_part->EndProcess() == "Decay" && g_part->NumberDaughters() == 2 && !_found_signature)
                {
                    fillParticle(mcp_h, g_part, _mcp_kshort);
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
                            for (const auto &dtr : daughters) 
                            {
                                unsigned int n_elas = 0;
                                unsigned int n_inelas = 0;

                                art::Ptr<simb::MCParticle> scat_part;
                                std::string scat_end_process;

                                common::GetNScatters(mcp_h, dtr, scat_part, n_elas, n_inelas);
                                scat_end_process = common::GetEndState(dtr, mcp_h);

                                if (dtr->PdgCode() == 211) // pion-plus
                                    fillParticle(mcp_h, dtr, _mcp_piplus);
                                
                                else if (dtr->PdgCode() == -211) // pion-minus
                                    fillParticle(mcp_h, dtr, _mcp_piminus);
                            }

                            _found_signature = std::all_of(daughters.begin(), daughters.end(), [&](const auto &dtr) 
                            {
                                return TVector3(dtr->Px(), dtr->Py(), dtr->Pz()).Mag() >= _PionThreshold;
                            });
                        }
                    }
                }
            }
        }
    }
}

void EventCategoryAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected)
{
    if (selected)
        _pass_preselection = true;

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

    if (_found_signature)
    {
        art::ValidHandle<std::vector<recob::Hit>> in_hits = e.getValidHandle<std::vector<recob::Hit>>(_Hproducer);
        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> mcp_bkth_assoc = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(in_hits, e, _BacktrackTag));

        int mcp_mu_hits = 0, mcp_piplus_hits = 0, mcp_piminus_hits = 0;
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

                if (mcp->TrackId() == _mcp_mu.tid)
                    mcp_mu_hits++;
                else if (mcp->TrackId() == _mcp_piplus.tid)
                    mcp_piplus_hits++;
                else if (mcp->TrackId() == _mcp_piminus.tid)
                    mcp_piminus_hits++;
            }
        }

        for (const common::ProxyPfpElem_t &pfp_pxy : slice_pfp_v)
        {
            if (pfp_pxy->IsPrimary())
            {
                double xyz[3] = {};
                auto vtx = pfp_pxy.get<recob::Vertex>();
                if (vtx.size() == 1)
                {
                    vtx.at(0)->XYZ(xyz);
                    auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

                    float _reco_nu_vtx_sce[3];
                    common::ApplySCECorrectionXYZ(nuvtx.X(), nuvtx.Y(), nuvtx.Z(), _reco_nu_vtx_sce);
                    _reco_nu_vtx_sce_x = _reco_nu_vtx_sce[0];
                    _reco_nu_vtx_sce_y = _reco_nu_vtx_sce[1];
                    _reco_nu_vtx_sce_z = _reco_nu_vtx_sce[2];

                    _reco_nu_vtx_sce_u_wire = (common::ProjectToWireView(_reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z, common::TPC_VIEW_U)).Z();
                    _reco_nu_vtx_sce_v_wire = (common::ProjectToWireView(_reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z, common::TPC_VIEW_V)).Z();
                    _reco_nu_vtx_sce_w_wire = (common::ProjectToWireView(_reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z, common::TPC_VIEW_W)).Z();
                }

                continue;
            }

            float trkscore = common::GetTrackShowerScore(pfp_pxy);
            if ((trkscore >= 0) && (trkscore >= 0.5))
            {
                _n_tracks++;
            }
            else if ((trkscore >= 0) && (trkscore < 0.5))
            {
                _n_showers++;
            }

            std::vector<art::Ptr<recob::Hit>> pfp_hits;
            auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();

            for (auto ass_clus : clus_pxy_v)
            {
                const auto &clus = clus_proxy[ass_clus.key()];
                auto clus_hit_v = clus.get<recob::Hit>();

                for (const auto &hit : clus_hit_v)
                    pfp_hits.push_back(hit);
            } 

            std::vector<simb::MCParticle> assmcp;
            std::vector<anab::BackTrackerHitMatchingData> amd;

            int pfp_mu_hits = 0, pfp_piplus_hits = 0, pfp_piminus_hits = 0;
            for (auto hit : pfp_hits)
            {
                auto assmcp = mcp_bkth_assoc->at(hit.key());       
                auto assmdt = mcp_bkth_assoc->data(hit.key());    
                for (size_t i = 0; i < assmcp.size(); i++)
                {
                    if (assmdt[i]->isMaxIDE != 1) 
                        continue;

                    const int track_id = assmcp[i]->TrackId();
                    if (track_id == _mcp_mu.tid)
                        pfp_mu_hits++;
                    else if (track_id == _mcp_piplus.tid)
                        pfp_piplus_hits++;
                    else if (track_id == _mcp_piminus.tid)
                        pfp_piminus_hits++;
                }
            }

            float muon_purity = (pfp_hits.size() > 0) ? static_cast<float>(pfp_mu_hits) / pfp_hits.size() : 0.0f;
            float muon_comple = (mcp_mu_hits > 0) ? static_cast<float>(pfp_mu_hits) / mcp_mu_hits : 0.0f;

            float piplus_purity = (pfp_hits.size() > 0) ? static_cast<float>(pfp_piplus_hits) / pfp_hits.size() : 0.0f;
            float piplus_comple = (mcp_piplus_hits > 0) ? static_cast<float>(pfp_piplus_hits) / mcp_piplus_hits : 0.0f;

            float piminus_purity = (pfp_hits.size() > 0) ? static_cast<float>(pfp_piminus_hits) / pfp_hits.size() : 0.0f;
            float piminus_comple = (mcp_piminus_hits > 0) ? static_cast<float>(pfp_piminus_hits) / mcp_piminus_hits : 0.0f;

            _mcp_pfp_id_v.push_back(pfp_pxy->Self());
            _mcp_mu_purity_v.push_back(muon_purity);
            _mcp_mu_comple_v.push_back(muon_comple);
            _mcp_piplus_purity_v.push_back(piplus_purity);
            _mcp_piplus_comple_v.push_back(piplus_comple);
            _mcp_piminus_purity_v.push_back(piminus_purity);
            _mcp_piminus_comple_v.push_back(piminus_comple);
        }
    }

    return;
}

void EventCategoryAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("pass_preselection", &_pass_preselection, "pass_preselection/B");

    _tree->Branch("nu_pdg", &_neutrino.nu_pdg, "nu_pdg/I");
    _tree->Branch("ccnc", &_neutrino.ccnc, "ccnc/I");
    _tree->Branch("nu_e", &_neutrino.nu_e, "nu_e/F");
    _tree->Branch("nu_theta", &_neutrino.nu_theta, "nu_theta/F");
    _tree->Branch("nu_pt", &_neutrino.nu_pt, "nu_pt/F");
    _tree->Branch("target_nucleus", &_neutrino.target_nucleus, "target_nucleus/I");
    _tree->Branch("hit_nucleon", &_neutrino.hit_nucleon, "hit_nucleon/I");
    _tree->Branch("W", &_neutrino.W, "W/F");
    _tree->Branch("X", &_neutrino.X, "X/F");
    _tree->Branch("Y", &_neutrino.Y, "Y/F");
    _tree->Branch("QSqr", &_neutrino.QSqr, "QSqr/F");
    
    _tree->Branch("true_vtx_sce_x", &_neutrino.true_vtx_sce_x, "true_vtx_sce_x/F");
    _tree->Branch("true_vtx_sce_y", &_neutrino.true_vtx_sce_y, "true_vtx_sce_y/F");
    _tree->Branch("true_vtx_sce_z", &_neutrino.true_vtx_sce_z, "true_vtx_sce_z/F");
    _tree->Branch("true_vtx_sce_u_wire", &_neutrino.true_vtx_sce_u_wire, "true_vtx_sce_u_wire/F");
    _tree->Branch("true_vtx_sce_v_wire", &_neutrino.true_vtx_sce_v_wire, "true_vtx_sce_v_wire/F");
    _tree->Branch("true_vtx_sce_w_wire", &_neutrino.true_vtx_sce_w_wire, "true_vtx_sce_w_wire/F");

    _tree->Branch("reco_vtx_sce_x", &_reco_nu_vtx_sce_x, "reco_vtx_sce_x/F");
    _tree->Branch("reco_vtx_sce_y", &_reco_nu_vtx_sce_y, "reco_vtx_sce_y/F");
    _tree->Branch("reco_vtx_sce_z", &_reco_nu_vtx_sce_z, "reco_vtx_sce_z/F");
    _tree->Branch("reco_vtx_sce_u_wire", &_reco_nu_vtx_sce_u_wire, "reco_vtx_sce_u_wire/F");
    _tree->Branch("reco_vtx_sce_v_wire", &_reco_nu_vtx_sce_v_wire, "reco_vtx_sce_v_wire/F");
    _tree->Branch("reco_vtx_sce_w_wire", &_reco_nu_vtx_sce_w_wire, "reco_vtx_sce_w_wire/F");

    _tree->Branch("n_tracks", &_n_tracks, "n_tracks/I");
    _tree->Branch("n_showers", &_n_showers, "n_showers/I");

    setParticleBranches(_tree, "mcp_muon", _mcp_mu);
    setParticleBranches(_tree, "mcp_kshort", _mcp_kshort);
    setParticleBranches(_tree, "mcp_piplus", _mcp_piplus);
    setParticleBranches(_tree, "mcp_piminus", _mcp_piminus);

    _tree->Branch("mcp_final_state", &_final_state);

    _tree->Branch("mcp_pfp_id", &_mcp_pfp_id_v);
    _tree->Branch("mcp_muon_purity", &_mcp_mu_purity_v);
    _tree->Branch("mcp_muon_comple", &_mcp_mu_comple_v); 
    _tree->Branch("mcp_piplus_purity", &_mcp_piplus_purity_v);
    _tree->Branch("mcp_piplus_comple", &_mcp_piplus_comple_v);
    _tree->Branch("mcp_piminus_purity", &_mcp_piminus_purity_v);
    _tree->Branch("mcp_piminus_comple", &_mcp_piminus_comple_v);
}

void EventCategoryAnalysis::setParticleBranches(TTree *_tree, const std::string &prefix, Particle &particle)
{
    _tree->Branch((prefix + "_tid").c_str(), &particle.tid, (prefix + "_tid/I").c_str());
    _tree->Branch((prefix + "_pdg").c_str(), &particle.pdg, (prefix + "_pdg/I").c_str());
    _tree->Branch((prefix + "_energy").c_str(), &particle.energy, (prefix + "_energy/F").c_str());
    _tree->Branch((prefix + "_px").c_str(), &particle.px, (prefix + "_px/F").c_str());
    _tree->Branch((prefix + "_py").c_str(), &particle.py, (prefix + "_py/F").c_str());
    _tree->Branch((prefix + "_pz").c_str(), &particle.pz, (prefix + "_pz/F").c_str());
    _tree->Branch((prefix + "_startx").c_str(), &particle.startx, (prefix + "_startx/F").c_str());
    _tree->Branch((prefix + "_starty").c_str(), &particle.starty, (prefix + "_starty/F").c_str());
    _tree->Branch((prefix + "_startz").c_str(), &particle.startz, (prefix + "_startz/F").c_str());
    _tree->Branch((prefix + "_endx").c_str(), &particle.endx, (prefix + "_endx/F").c_str());
    _tree->Branch((prefix + "_endy").c_str(), &particle.endy, (prefix + "_endy/F").c_str());
    _tree->Branch((prefix + "_endz").c_str(), &particle.endz, (prefix + "_endz/F").c_str());
    _tree->Branch((prefix + "_n_elas").c_str(), &particle.n_elas, (prefix + "_n_elas/I").c_str());
    _tree->Branch((prefix + "_n_inelas").c_str(), &particle.n_inelas, (prefix + "_n_inelas/I").c_str());
    _tree->Branch((prefix + "_endstate").c_str(), &particle.endstate);
}

void EventCategoryAnalysis::resetTTree(TTree *_tree)
{
    _pass_preselection = false;

    _neutrino = Neutrino();
    _mcp_mu = Particle();
    _mcp_kshort = Particle();
    _mcp_piplus = Particle();
    _mcp_piminus = Particle();

    _final_state.clear();

    _reco_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_z = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_u_wire = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_v_wire = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_w_wire = std::numeric_limits<float>::lowest();

    _n_tracks = 0;
    _n_showers = 0;

    _mcp_pfp_id_v.clear();
    _mcp_mu_purity_v.clear();
    _mcp_mu_comple_v.clear();
    _mcp_piplus_purity_v.clear();
    _mcp_piplus_comple_v.clear();
    _mcp_piminus_purity_v.clear();
    _mcp_piminus_comple_v.clear();
}

DEFINE_ART_CLASS_TOOL(EventCategoryAnalysis)
} 

#endif