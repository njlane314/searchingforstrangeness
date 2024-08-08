#ifndef ANALYSIS_DEFAULTANALYSIS_CXX
#define ANALYSIS_DEFAULTANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "../CommonFuncs/BackTrackingFuncs.h"
#include "../CommonFuncs/Geometry.h"
#include "../CommonFuncs/SpaceChargeCorrections.h"
#include "../CommonFuncs/Containment.h"
#include "../CommonFuncs/TrackShowerScoreFuncs.h"
#include "../CommonFuncs/ProximityClustering.h"

#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "canvas/Persistency/Common/TriggerResults.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace analysis
{

class StrangenessAnalysis_tool : public AnalysisToolBase
{

public:
    StrangenessAnalysis_tool(const fhicl::ParameterSet &pset);
    ~StrangenessAnalysis_tool(){};

    void configure(fhicl::ParameterSet const &pset);
    void analyzeEvent(art::Event const &e, bool fData) override;
    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    void SaveTruth(art::Event const &e);

    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
    TParticlePDG *neutron = TDatabasePDG::Instance()->GetParticle(2112);
    TParticlePDG *electron = TDatabasePDG::Instance()->GetParticle(11);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
    TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);
    TParticlePDG *pi0 = TDatabasePDG::Instance()->GetParticle(111);
    TParticlePDG *electron_neutrino = TDatabasePDG::Instance()->GetParticle(12);
    TParticlePDG *muon_neutrino = TDatabasePDG::Instance()->GetParticle(14);
    TParticlePDG *kaon_short = TDatabasePDG::Instance()->GetParticle(310);
    TParticlePDG *kaon_long = TDatabasePDG::Instance()->GetParticle(130);
    TParticlePDG *lambda = TDatabasePDG::Instance()->GetParticle(3122);

    art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
    art::InputTag fCLSproducer;     // cluster associated to PFP
    art::InputTag fMCTproducer;     // MCTruth from neutrino generator
    art::InputTag fMCPproducer;     // MCParticle from Geant4 stage
    art::InputTag fMCFluxproducer;  // MCFlux producer
    art::InputTag fBacktrackTag;
    art::InputTag fHproducer;
    art::InputTag fMCRproducer;
    art::InputTag fSLCproducer; // slice associated to PFP
    float fTrkShrScore;         /**< Threshold on the Pandora track score (default 0.5) */

    float fFidvolXstart;
    float fFidvolXend;
    float fFidvolYstart;
    float fFidvolYend;
    float fFidvolZstart;
    float fFidvolZend;

    const int k_nu_mu_cc_kshrt_single = 1;
    const int k_nu_mu_cc_kshrt_associated = 2;
    const int k_nu_mu_cc_only_hyperon = 3;
    const int k_nu_mu_other = 4;
    const int k_nc = 5;
    const int k_cosmic = 6;
    const int k_outfv = 7;
    const int k_other = 8;
    const int k_data = 0;

    float fProtonThreshold;
    float fPionThreshold;
    float fElectronThreshold;
    float fMuonThreshold;

    int _category;

    float _true_nu_vtx_t, _true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z;
    float _true_nu_vtx_sce_x, _true_nu_vtx_sce_y, _true_nu_vtx_sce_z;
    float _reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z;
    float _reco_nu_vtx_sce_x, _reco_nu_vtx_sce_y, _reco_nu_vtx_sce_z;

    float _nu_e;  /**< neutrino energy [GeV] */
    float _nu_l;  /**< propagation length [m] */
    float _nu_pt; /**< transverse momentum of interaction [GeV/c] */
    float _theta; /**< angle between incoming and outgoing leptons, in radians */

    int _nu_pdg;           /**< neutrino PDG code */
    int _ccnc;             /**< CC or NC tag from GENIE */
    int _nu_parent_pdg;    /**< neutrino parent's PDG code [http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/] */
    int _nu_hadron_pdg;    /**< PDG code of hadron eventually producing neutrino [http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/] */
    int _nu_decay_mode;    /**< decay mode that lead to this neutrino in beam simulation [http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/] */
    int _interaction;      /**< Interaction code from GENIE */
    bool _isVtxInFiducial; /**< true if neutrino in fiducial volume */
    bool _truthFiducial;   /**< is the truth information contained. Require all track start/end point in FV and showers deposit > 60% of energy in TPC or deposit at least 100 MeV in TPC */

    int _nmuon;                         /**< is there a final-state muon from the neutrino? [1=yes 0=no] */
    float _muon_e, _muon_p, _muon_c;    /**< energy, purity, completeness. */
    int _nelec;                         /**< is there a final-state electron from the neutrino? [1=yes 0=no] */
    float _elec_e, _pi_p, _pi_c;    /**< energy, purity, completeness. */
    float _elec_vx, _elec_vy, _elec_vz; /**< electron vertex. */
    float _elec_px, _elec_py, _elec_pz; /**< electron momentum vector [normalized] */
    int _nneutron;                      /**< how many neutrons are there? */
    int _nproton;                        /**< how many protons are there? */
    int _proton;                           /**< is there a final-state proton from the neutrino? [1=yes 0=no] */
    float _proton_e, _proton_p, _proton_c; /**< energy, purity, completeness. */
    int _npion;                            /**< how many pions are there? */
    int _pion;                             /**< is there a final-state charged pion from the neutrino? [1=yes 0=no] */
    float _pion_e, _pion_p, _pion_c;       /**< energy, purity, completeness. */
    int _nkshrt;                             /**< is there a final-state kaon short from the neutrino? [1=yes 0=no] */
    int _kshrt;                             /**< is there a final-state kaon short from the neutrino*/
    float _kshrt_e;                         /**< energy */
    int _nlambda;                          /**< how many lambdas are there? */
    int _lambda;                             /**< is there a final-state lambda from the neutrino? */
    float _lambda_e;                        /**< energy */

    std::string _endkshrtprocess; /**< End muon process name */
    float _endkshrthadronic;       /**< End muon Michel electron energy */

    std::vector<int> _pfp_slice_idx; /**< index of PFP is vector of PFPs in nu slice */

    std::vector<int>   _backtracked_idx;            // index of PFP [key]
    std::vector<int>   _backtracked_pdg;            // PDG code of backtracked particle
    std::vector<float> _backtracked_e;              // energy of backtracked particle
    std::vector<int>   _backtracked_tid;            // track-id of backtracked particle
    std::vector<float> _backtracked_purity;         // purity of backtracking
    std::vector<float> _backtracked_completeness;   // completeness of backtracking
    std::vector<float> _backtracked_overlay_purity; // purity of overlay

    std::vector<float> _backtracked_px;
    std::vector<float> _backtracked_py;
    std::vector<float> _backtracked_pz;

    std::vector<float> _backtracked_start_x;
    std::vector<float> _backtracked_start_y;
    std::vector<float> _backtracked_start_z;
    std::vector<float> _backtracked_start_U;
    std::vector<float> _backtracked_start_V;
    std::vector<float> _backtracked_start_Y;
    std::vector<float> _backtracked_sce_start_x;
    std::vector<float> _backtracked_sce_start_y;
    std::vector<float> _backtracked_sce_start_z;
    std::vector<float> _backtracked_sce_start_U;
    std::vector<float> _backtracked_sce_start_V;
    std::vector<float> _backtracked_sce_start_Y;

    float _lep_e; // lepton energy (if one exists) [GeV]
    int _pass;    // does the slice pass the selection

    int evnhits;                     // number of hits in event
    int slpdg;                       // PDG code of primary pfp in slice
    int slnhits;                     // number of hits in slice
    float _topo_score;               /**< topological score of the slice */
    std::vector<int> pfpdg;          // PDG code of pfp in slice
    std::vector<int> pfnhits;        // number of hits in pfp
    std::vector<int> pfnplanehits_U; // number of hits in pfp plane U
    std::vector<int> pfnplanehits_V; // number of hits in pfp plane V
    std::vector<int> pfnplanehits_Y; // number of hits in pfp plane Y
    std::vector<int> pfpplanesubclusters_U;
    std::vector<int> pfpplanesubclusters_V;
    std::vector<int> pfpplanesubclusters_Y;
    std::vector<float> pfpplanesubhitfracmax_U;
    std::vector<float> pfpplanesubhitfracmax_V;
    std::vector<float> pfpplanesubhitfracmax_Y;
    float slclustfrac; //fraction of clustered hits in the slice

    std::vector<uint> _generation;    // generation, 1 is primary
    std::vector<uint> _shr_daughters; // number of shower daughters
    std::vector<uint> _trk_daughters; // number of track daughters

    unsigned int _n_pfps;
    std::vector<float> _trk_score_v;
    unsigned int _n_tracks;
    unsigned int _n_showers;

    unsigned int _hits_u;
    unsigned int _hits_v;
    unsigned int _hits_y;

    std::vector<int> _mc_pdg;
    std::vector<float> _mc_E;

    std::vector<float> _mc_px;
    std::vector<float> _mc_py;
    std::vector<float> _mc_pz;

    std::vector<float> _mc_vx;
    std::vector<float> _mc_vy;
    std::vector<float> _mc_vz;

    std::vector<float> _mc_endx;
    std::vector<float> _mc_endy;
    std::vector<float> _mc_endz;

    std::vector<float> _mc_completeness;
    std::vector<float> _mc_purity;

    float _true_pt;
    float _true_pt_visible;
    float _true_p;
    float _true_p_visible;
    float _true_e_visible;
    float _leeweight;
};

StrangenessAnalysis_tool::StrangenessAnalysis_tool(const fhicl::ParameterSet &p)
{
    fCRTVetoproducer = p.get<art::InputTag>("CRTVetoproducer", ""); // default is no CRT veto
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fMCFluxproducer = p.get<art::InputTag>("MCFluxproducer");
    fBacktrackTag = p.get<art::InputTag>("BacktrackTag");
    fHproducer = p.get<art::InputTag>("Hproducer");
    fMCRproducer = p.get<art::InputTag>("MCRproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fTrkShrScore = p.get<float>("TrkShrScore", 0.5);

    fProtonThreshold = p.get<float>("ProtonThreshold", 0.04);
    fMuonThreshold = p.get<float>("MuonThreshold", 0.02);
    fPionThreshold = p.get<float>("PionThreshold", 0.04);
    fElectronThreshold = p.get<float>("ElectronThreshold", 0.03);

    fFidvolXstart = p.get<double>("fidvolXstart", 10);
    fFidvolXend = p.get<double>("fidvolXend", 10);

    fFidvolYstart = p.get<double>("fidvolYstart", 10);
    fFidvolYend = p.get<double>("fidvolYend", 10);

    fFidvolZstart = p.get<double>("fidvolZstart", 10);
    fFidvolZend = p.get<double>("fidvolZend", 50);
}

void StrangenessAnalysis_tool::configure(fhicl::ParameterSet const &p)
{
}

void StrangenessAnalysis_tool::analyzeEvent(art::Event const &e, bool fData)
{
    std::cout << "[StrangenessAnalysis_tool::analyzeEvent] Run: " << e.run() << ", SubRun: " << e.subRun() << ", Event: " << e.event() << std::endl;

    if (!fData)
    {
        SaveTruth(e);
    }

    art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
    evnhits = inputHits->size();
}

void StrangenessAnalysis_tool::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));

    art::ValidHandle<std::vector<recob::Slice>> inputSlice = e.getValidHandle<std::vector<recob::Slice>>(fSLCproducer);
    auto assocSliceHit = std::unique_ptr<art::FindManyP<recob::Hit>>(new art::FindManyP<recob::Hit>(inputSlice, e, fSLCproducer));

    lar_pandora::LArPandoraHelper larpandora;
    lar_pandora::PFParticleVector pfparticles;
    lar_pandora::PFParticleMap particleMap;
    larpandora.CollectPFParticles(e, "pandora", pfparticles);
    larpandora.BuildPFParticleMap(pfparticles, particleMap);

    std::vector<common::BtPart> btparts_v;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!fData)
    {
        const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer));
        const std::vector<sim::MCTrack> &inputMCTrack = *(e.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer));
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
        btparts_v = common::initBacktrackingParticleVec(inputMCShower, inputMCTrack, *inputHits, assocMCPart);
    }

    size_t pfpidx = 0;
    _n_pfps = 0;
    for (auto pfp : slice_pfp_v)
    {
        if (pfp->IsPrimary())
        {
            slpdg = pfp->PdgCode();
            auto slice_pxy_v = pfp.get<recob::Slice>();
            if (slice_pxy_v.size() != 1)
            {
                std::cout << "WRONG!!! n slices = " << slice_pxy_v.size() << " " << __FILE__ << " " << __LINE__ << std::endl;
                return;
            }

            auto slicehits = assocSliceHit->at(slice_pxy_v[0].key());
            slnhits = slicehits.size();

            auto metadata_pxy_v = pfp.get<larpandoraobj::PFParticleMetadata>();

            if (metadata_pxy_v.size() != 0)
            {
                for (unsigned int j = 0; j < metadata_pxy_v.size(); ++j)
                {
                    const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(metadata_pxy_v.at(j));
                    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
                    if (!pfParticlePropertiesMap.empty())
                    {
                        _topo_score = pfParticlePropertiesMap.at("NuScore");
                    } 
                }
            }

            double xyz[3] = {};

            auto vtx = pfp.get<recob::Vertex>();
            if (vtx.size() == 1)
            {
                vtx.at(0)->XYZ(xyz);
                auto nuvtx = TVector3(xyz[0], xyz[1], xyz[2]);

                _reco_nu_vtx_x = nuvtx.X();
                _reco_nu_vtx_y = nuvtx.Y();
                _reco_nu_vtx_z = nuvtx.Z();

                float _reco_nu_vtx_sce[3];
                common::ApplySCECorrectionXYZ(_reco_nu_vtx_x, _reco_nu_vtx_y, _reco_nu_vtx_z, _reco_nu_vtx_sce);
                _reco_nu_vtx_sce_x = _reco_nu_vtx_sce[0];
                _reco_nu_vtx_sce_y = _reco_nu_vtx_sce[1];
                _reco_nu_vtx_sce_z = _reco_nu_vtx_sce[2];
            }
            else
            {
                std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
            }

            continue;
        } 

        _n_pfps++;
        _pfp_slice_idx.push_back(pfpidx++);
        pfpdg.push_back(pfp->PdgCode());

        _generation.push_back(larpandora.GetGeneration(particleMap, particleMap.at(pfp->Self())));
        uint this_num_trk_d = 0;
        uint this_num_shr_d = 0;
        for (size_t daughter : pfp->Daughters())
        {
            if (larpandora.IsTrack(particleMap.at(daughter)))
            {
                this_num_trk_d++; // Track daughter
            }
            else
            {
                this_num_shr_d++; // Shower daughter
            }
        }
        _shr_daughters.push_back(this_num_shr_d);
        _trk_daughters.push_back(this_num_trk_d);

        float trkscore = common::GetTrackShowerScore(pfp);
        if ((trkscore >= 0) && (trkscore >= fTrkShrScore))
        {
            _n_tracks++;
        }
        else if ((trkscore >= 0) && (trkscore < fTrkShrScore))
        {
            _n_showers++;
        }
        _trk_score_v.push_back(trkscore);

        std::vector<art::Ptr<recob::Hit>> hit_v;
        pfnhits.push_back(hit_v.size());

        if (!fData)
        {
            float purity = 0., completeness = 0., overlay_purity = 0.;
            int ibt = common::getAssocBtPart(hit_v, assocMCPart, btparts_v, purity, completeness, overlay_purity);
            if (ibt >= 0)
            {
                auto &mcp = btparts_v[ibt];
                auto PDG = mcp.pdg;
                _backtracked_idx.push_back(pfp->Self());
                _backtracked_e.push_back(mcp.e);
                _backtracked_tid.push_back(mcp.tids.at(0));
                _backtracked_pdg.push_back(PDG);
                _backtracked_purity.push_back(purity);
                _backtracked_completeness.push_back(completeness);
                _backtracked_overlay_purity.push_back(overlay_purity);

                _backtracked_px.push_back(mcp.px);
                _backtracked_py.push_back(mcp.py);
                _backtracked_pz.push_back(mcp.pz);
                _backtracked_start_x.push_back(mcp.start_x);
                _backtracked_start_y.push_back(mcp.start_y);
                _backtracked_start_z.push_back(mcp.start_z);

                for (size_t i = 0; i < _mc_E.size(); i++)
                {
                    if (_mc_E[i] == mcp.e)
                    {
                        _mc_completeness[i] = completeness;
                        _mc_purity[i] = purity;
                    }
                }

                float reco_st[3] = {mcp.start_x, mcp.start_y, mcp.start_z};

                if (PDG == 11 || PDG == 22)
                {
                    reco_st[0] += common::x_offset(mcp.start_t);
                }
                else
                {
                    common::True2RecoMappingXYZ(mcp.start_t, mcp.start_x, mcp.start_y, mcp.start_z, reco_st);
                }

                _backtracked_sce_start_x.push_back(reco_st[0]);
                _backtracked_sce_start_y.push_back(reco_st[1]);
                _backtracked_sce_start_z.push_back(reco_st[2]);

                _backtracked_sce_start_U.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 0));
                _backtracked_sce_start_V.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 1));
                _backtracked_sce_start_Y.push_back(common::YZtoPlanecoordinate(reco_st[1], reco_st[2], 2));
            }
            else
            {
                _backtracked_idx.push_back(0);
                _backtracked_tid.push_back(0);
                _backtracked_e.push_back(std::numeric_limits<float>::lowest());
                _backtracked_tid.push_back(std::numeric_limits<int>::lowest());
                _backtracked_pdg.push_back(0);
                _backtracked_purity.push_back(std::numeric_limits<float>::lowest());
                _backtracked_completeness.push_back(std::numeric_limits<float>::lowest());
                _backtracked_overlay_purity.push_back(std::numeric_limits<float>::lowest());

                _backtracked_px.push_back(std::numeric_limits<float>::lowest());
                _backtracked_py.push_back(std::numeric_limits<float>::lowest());
                _backtracked_pz.push_back(std::numeric_limits<float>::lowest());
                _backtracked_start_x.push_back(std::numeric_limits<float>::lowest());
                _backtracked_start_y.push_back(std::numeric_limits<float>::lowest());
                _backtracked_start_z.push_back(std::numeric_limits<float>::lowest());

                _backtracked_start_U.push_back(std::numeric_limits<float>::lowest());
                _backtracked_start_V.push_back(std::numeric_limits<float>::lowest());
                _backtracked_start_Y.push_back(std::numeric_limits<float>::lowest());

                _backtracked_sce_start_x.push_back(std::numeric_limits<float>::lowest());
                _backtracked_sce_start_y.push_back(std::numeric_limits<float>::lowest());
                _backtracked_sce_start_z.push_back(std::numeric_limits<float>::lowest());

                _backtracked_sce_start_U.push_back(std::numeric_limits<float>::lowest());
                _backtracked_sce_start_V.push_back(std::numeric_limits<float>::lowest());
                _backtracked_sce_start_Y.push_back(std::numeric_limits<float>::lowest());
            }
        }   
    }

    if (!fData)
    {
        bool there_is_true_mu = _nmuon > 0;
        bool there_is_true_kaon_short = _nkshrt > 0;
        bool there_is_true_lambda = _nlambda > 0;

        if (!_isVtxInFiducial)
        {
            _category = k_outfv;
        }
        else if (abs(_nu_pdg) == muon_neutrino->PdgCode())
        {
            if (there_is_true_mu)
            {
                if (there_is_true_kaon_short && !there_is_true_lambda)
                {
                    _category = k_nu_mu_cc_kshrt_single;
                }
                else if (there_is_true_kaon_short && there_is_true_lambda)
                {
                    _category = k_nu_mu_cc_kshrt_associated;
                }
                else if (!there_is_true_kaon_short && there_is_true_lambda)
                {
                    _category = k_nu_mu_cc_only_hyperon;
                }
                else
                {   
                    _category = k_nu_mu_other; 
                }
            }
            else
            {
                _category = k_nc;
            }
        }
        else
        {
            _category = k_cosmic;
        }
    }
    else
    {
        _category = k_data;
    }

    if (selected)
        _pass = 1;
}   

void StrangenessAnalysis_tool::setBranches(TTree *_tree)
{
    _tree->Branch("leeweight", &_leeweight, "leeweight/F");

    _tree->Branch("true_pt", &_true_pt, "true_pt/F");
    _tree->Branch("true_pt_visible", &_true_pt_visible, "true_pt_visible/F");
    _tree->Branch("true_p", &_true_p, "true_p/F");
    _tree->Branch("true_p_visible", &_true_p_visible, "true_p_visible/F");

    _tree->Branch("true_e_visible", &_true_e_visible, "true_e_visible/F");

    // neutrino information
    _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
    _tree->Branch("ccnc", &_ccnc, "ccnc/I");
    _tree->Branch("nu_parent_pdg", &_nu_parent_pdg, "nu_parent_pdg/I");
    _tree->Branch("nu_hadron_pdg", &_nu_hadron_pdg, "nu_hadron_pdg/I");
    _tree->Branch("nu_decay_mode", &_nu_decay_mode, "nu_decay_mode/I");
    _tree->Branch("interaction", &_interaction, "interaction/I");
    _tree->Branch("nu_e", &_nu_e, "nu_e/F");
    _tree->Branch("nu_l", &_nu_l, "nu_l/F");
    _tree->Branch("nu_pt", &_nu_pt, "nu_pt/F");
    _tree->Branch("theta", &_theta, "theta/F");
    _tree->Branch("isVtxInFiducial", &_isVtxInFiducial, "isVtxInFiducial/O");
    _tree->Branch("truthFiducial", &_truthFiducial, "truthFiducial/O");

    _tree->Branch("true_nu_vtx_t", &_true_nu_vtx_t, "true_nu_vtx_t/F");
    _tree->Branch("true_nu_vtx_x", &_true_nu_vtx_x, "true_nu_vtx_x/F");
    _tree->Branch("true_nu_vtx_y", &_true_nu_vtx_y, "true_nu_vtx_y/F");
    _tree->Branch("true_nu_vtx_z", &_true_nu_vtx_z, "true_nu_vtx_z/F");
    _tree->Branch("true_nu_vtx_sce_x", &_true_nu_vtx_sce_x, "true_nu_vtx_sce_x/F");
    _tree->Branch("true_nu_vtx_sce_y", &_true_nu_vtx_sce_y, "true_nu_vtx_sce_y/F");
    _tree->Branch("true_nu_vtx_sce_z", &_true_nu_vtx_sce_z, "true_nu_vtx_sce_z/F");
    _tree->Branch("reco_nu_vtx_x", &_reco_nu_vtx_x, "reco_nu_vtx_x/F");
    _tree->Branch("reco_nu_vtx_y", &_reco_nu_vtx_y, "reco_nu_vtx_y/F");
    _tree->Branch("reco_nu_vtx_z", &_reco_nu_vtx_z, "reco_nu_vtx_z/F");
    _tree->Branch("reco_nu_vtx_sce_x", &_reco_nu_vtx_sce_x, "reco_nu_vtx_sce_x/F");
    _tree->Branch("reco_nu_vtx_sce_y", &_reco_nu_vtx_sce_y, "reco_nu_vtx_sce_y/F");
    _tree->Branch("reco_nu_vtx_sce_z", &_reco_nu_vtx_sce_z, "reco_nu_vtx_sce_z/F");

    // individual particles in the neutrino slice
    // legend:
    // _e -> energy of particle in GeV
    // _c -> completeness from back-tracking [0,1]
    // _p -> purity from back-tracking [0,1]
    // muon
    _tree->Branch("nmuon", &_nmuon, "nmuon/I");
    _tree->Branch("muon_e", &_muon_e, "muon_e/F");
    _tree->Branch("muon_c", &_muon_c, "muon_c/F");
    _tree->Branch("muon_p", &_muon_p, "muon_p/F");
    // electron
    _tree->Branch("nelec", &_nelec, "nelec/I");
    _tree->Branch("elec_e", &_elec_e, "elec_e/F");
    _tree->Branch("elec_c", &_pi_c, "elec_c/F");
    _tree->Branch("elec_p", &_pi_p, "elec_p/F");
    _tree->Branch("elec_vx", &_elec_vx, "elec_vx/F");
    _tree->Branch("elec_vy", &_elec_vy, "elec_vy/F");
    _tree->Branch("elec_vz", &_elec_vz, "elec_vz/F");
    _tree->Branch("elec_px", &_elec_px, "elec_px/F");
    _tree->Branch("elec_py", &_elec_py, "elec_py/F");
    _tree->Branch("elec_pz", &_elec_pz, "elec_pz/F");

    _tree->Branch("nneutron", &_nneutron, "nneutron/I");

    // first [highest momentum] proton
    _tree->Branch("nproton", &_nproton, "nproton/I");
    _tree->Branch("proton_e", &_proton_e, "proton_e/F");
    _tree->Branch("proton_c", &_proton_c, "proton_c/F");
    _tree->Branch("proton_p", &_proton_p, "proton_p/F");

    // charged pions
    _tree->Branch("npion", &_npion, "npion/I");
    _tree->Branch("pion_e", &_pion_e, "pion_e/F");
    _tree->Branch("pion_c", &_pion_c, "pion_c/F");
    _tree->Branch("pion_p", &_pion_p, "pion_p/F");

    _tree->Branch("pfp_slice_idx", "std::vector<int>", &_pfp_slice_idx);

    // PFParticle backtracking
    // _tree->Branch("backtracked_idx"   ,"std::vector<int>"  ,&_backtracked_idx   );
    // _tree->Branch("backtracked_tid"   ,"std::vector<int>"  ,&_backtracked_tid   );

    _tree->Branch("category", &_category, "category/I");

    _tree->Branch("backtracked_pdg", "std::vector<int>", &_backtracked_pdg);
    _tree->Branch("backtracked_e", "std::vector<float>", &_backtracked_e);
    _tree->Branch("backtracked_tid", "std::vector<int>", &_backtracked_tid);
    _tree->Branch("backtracked_purity", "std::vector<float>", &_backtracked_purity);
    _tree->Branch("backtracked_completeness", "std::vector<float>", &_backtracked_completeness);
    _tree->Branch("backtracked_overlay_purity", "std::vector<float>", &_backtracked_overlay_purity);

    _tree->Branch("backtracked_px", "std::vector<float>", &_backtracked_px);
    _tree->Branch("backtracked_py", "std::vector<float>", &_backtracked_py);
    _tree->Branch("backtracked_pz", "std::vector<float>", &_backtracked_pz);

    _tree->Branch("backtracked_start_x", "std::vector<float>", &_backtracked_start_x);
    _tree->Branch("backtracked_start_y", "std::vector<float>", &_backtracked_start_y);
    _tree->Branch("backtracked_start_z", "std::vector<float>", &_backtracked_start_z);
    _tree->Branch("backtracked_start_U", "std::vector<float>", &_backtracked_start_U);
    _tree->Branch("backtracked_start_V", "std::vector<float>", &_backtracked_start_V);
    _tree->Branch("backtracked_start_Y", "std::vector<float>", &_backtracked_start_Y);
    _tree->Branch("backtracked_sce_start_x", "std::vector<float>", &_backtracked_sce_start_x);
    _tree->Branch("backtracked_sce_start_y", "std::vector<float>", &_backtracked_sce_start_y);
    _tree->Branch("backtracked_sce_start_z", "std::vector<float>", &_backtracked_sce_start_z);
    _tree->Branch("backtracked_sce_start_U", "std::vector<float>", &_backtracked_sce_start_U);
    _tree->Branch("backtracked_sce_start_V", "std::vector<float>", &_backtracked_sce_start_V);
    _tree->Branch("backtracked_sce_start_Y", "std::vector<float>", &_backtracked_sce_start_Y);

    _tree->Branch("lep_e", &_lep_e, "lep_e/F");
    _tree->Branch("pass", &_pass, "pass/I");

    _tree->Branch("evnhits", &evnhits, "evnhits/I");
    _tree->Branch("slpdg", &slpdg, "slpdg/I");
    _tree->Branch("slnhits", &slnhits, "slnhits/I");
    _tree->Branch("n_pfps", &_n_pfps, "n_pfps/I");
    _tree->Branch("n_tracks", &_n_tracks, "n_tracks/I");
    _tree->Branch("n_showers", &_n_showers, "n_showers/I");

    _tree->Branch("pfp_generation_v", "std::vector< uint >", &_generation);
    _tree->Branch("pfp_trk_daughters_v", "std::vector< uint >", &_trk_daughters);
    _tree->Branch("pfp_shr_daughters_v", "std::vector< uint >", &_shr_daughters);

    _tree->Branch("trk_score_v", "std::vector< float >", &_trk_score_v);

    _tree->Branch("pfpdg", "std::vector<int>", &pfpdg);
    _tree->Branch("pfnhits", "std::vector<int>", &pfnhits);
    _tree->Branch("pfnplanehits_U", "std::vector<int>", &pfnplanehits_U);
    _tree->Branch("pfnplanehits_V", "std::vector<int>", &pfnplanehits_V);
    _tree->Branch("pfnplanehits_Y", "std::vector<int>", &pfnplanehits_Y);
    _tree->Branch("pfpplanesubclusters_U", "std::vector<int>", &pfpplanesubclusters_U);
    _tree->Branch("pfpplanesubclusters_V", "std::vector<int>", &pfpplanesubclusters_V);
    _tree->Branch("pfpplanesubclusters_Y", "std::vector<int>", &pfpplanesubclusters_Y);
    _tree->Branch("pfpplanesubhitfracmax_U", "std::vector<float>", &pfpplanesubhitfracmax_U);
    _tree->Branch("pfpplanesubhitfracmax_V", "std::vector<float>", &pfpplanesubhitfracmax_V);
    _tree->Branch("pfpplanesubhitfracmax_Y", "std::vector<float>", &pfpplanesubhitfracmax_Y);

    _tree->Branch("hits_u", &_hits_u, "hits_u/i");
    _tree->Branch("hits_v", &_hits_v, "hits_v/i");
    _tree->Branch("hits_y", &_hits_y, "hits_y/i");
    _tree->Branch("topological_score", &_topo_score, "topological_score/F");
    _tree->Branch("slclustfrac", &slclustfrac, "slclustfrac/F");

    _tree->Branch("mc_pdg", "std::vector< int >", &_mc_pdg);
    _tree->Branch("mc_E", "std::vector< float >", &_mc_E);

    _tree->Branch("mc_vx", "std::vector< float >", &_mc_vx);
    _tree->Branch("mc_vy", "std::vector< float >", &_mc_vy);
    _tree->Branch("mc_vz", "std::vector< float >", &_mc_vz);

    _tree->Branch("mc_endx", "std::vector< float >", &_mc_endx);
    _tree->Branch("mc_endy", "std::vector< float >", &_mc_endy);
    _tree->Branch("mc_endz", "std::vector< float >", &_mc_endz);

    _tree->Branch("mc_px", "std::vector< float >", &_mc_px);
    _tree->Branch("mc_py", "std::vector< float >", &_mc_py);
    _tree->Branch("mc_pz", "std::vector< float >", &_mc_pz);

    _tree->Branch("mc_completeness", "std::vector< float >", &_mc_completeness);
    _tree->Branch("mc_purity", "std::vector< float >", &_mc_purity);

    _tree->Branch("endmuonprocess", &_endkshrtprocess);

    _tree->Branch("endmuonmichel", &_endkshrthadronic, "endmuonmichel/F");
}

void StrangenessAnalysis_tool::resetTTree(TTree *_tree)
{
    _leeweight = 0;
    _nu_e = std::numeric_limits<float>::lowest();
    _nu_l = std::numeric_limits<float>::lowest();
    _theta = std::numeric_limits<float>::lowest();
    _nu_pt = std::numeric_limits<float>::lowest();

    _nu_pdg = std::numeric_limits<int>::lowest();
    _ccnc = std::numeric_limits<int>::lowest();
    _nu_parent_pdg = std::numeric_limits<int>::lowest();
    _nu_hadron_pdg = std::numeric_limits<int>::lowest();
    _nu_decay_mode = std::numeric_limits<int>::lowest();
    _interaction = std::numeric_limits<int>::lowest();
    _pass = 0;

    _category = 0;

    _true_nu_vtx_t = std::numeric_limits<float>::lowest();
    _true_nu_vtx_x = std::numeric_limits<float>::lowest();
    _true_nu_vtx_y = std::numeric_limits<float>::lowest();
    _true_nu_vtx_z = std::numeric_limits<float>::lowest();
    _true_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
    _true_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
    _true_nu_vtx_sce_z = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_x = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_y = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_z = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_x = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_y = std::numeric_limits<float>::lowest();
    _reco_nu_vtx_sce_z = std::numeric_limits<float>::lowest();

    _isVtxInFiducial = false;

    _nmuon = 0;
    _muon_e = 0;
    _muon_p = 0;
    _muon_c = 0;

    _nelec = 0;
    _elec_e = 0;
    _pi_p = 0;
    _pi_c = 0;
    _elec_vx = std::numeric_limits<float>::lowest();
    _elec_vy = std::numeric_limits<float>::lowest();
    _elec_vz = std::numeric_limits<float>::lowest();
    _elec_px = std::numeric_limits<float>::lowest();
    _elec_py = std::numeric_limits<float>::lowest();

    _npion = 0;
    _pion_e = 0;
    _pion_p = 0;
    _pion_c = 0;

    _nneutron = 0;

    _nproton = 0;
    _proton_e = 0;
    _proton_p = 0;
    _proton_c = 0;

    _endkshrtprocess = "";
    _endkshrthadronic = 0;

    _n_pfps = 0;
    _n_tracks = 0;
    _n_showers = 0;
    _trk_score_v.clear();
    // _backtracked_idx.clear();
    // _backtracked_tid.clear();
    _backtracked_e.clear();
    _backtracked_tid.clear();
    _backtracked_pdg.clear();
    _backtracked_purity.clear();
    _backtracked_completeness.clear();
    _backtracked_overlay_purity.clear();

    _backtracked_px.clear();
    _backtracked_py.clear();
    _backtracked_pz.clear();

    _backtracked_start_x.clear();
    _backtracked_start_y.clear();
    _backtracked_start_z.clear();
    _backtracked_start_U.clear();
    _backtracked_start_V.clear();
    _backtracked_start_Y.clear();
    _backtracked_sce_start_x.clear();
    _backtracked_sce_start_y.clear();
    _backtracked_sce_start_z.clear();
    _backtracked_sce_start_U.clear();
    _backtracked_sce_start_V.clear();
    _backtracked_sce_start_Y.clear();

    evnhits = std::numeric_limits<int>::lowest();
    slpdg = std::numeric_limits<int>::lowest();
    _topo_score = std::numeric_limits<float>::lowest();
    slnhits = std::numeric_limits<int>::lowest();
    pfpdg.clear();
    pfnhits.clear();
    pfnplanehits_U.clear();
    pfnplanehits_V.clear();
    pfnplanehits_Y.clear();
    pfpplanesubclusters_U.clear();
    pfpplanesubclusters_V.clear();
    pfpplanesubclusters_Y.clear();
    pfpplanesubhitfracmax_U.clear();
    pfpplanesubhitfracmax_V.clear();
    pfpplanesubhitfracmax_Y.clear();
    _generation.clear();
    _shr_daughters.clear();
    _trk_daughters.clear();
    slclustfrac = std::numeric_limits<float>::lowest();

    _hits_u = 0;
    _hits_v = 0;
    _hits_y = 0;

    _mc_E.clear();
    _mc_pdg.clear();

    _mc_px.clear();
    _mc_py.clear();
    _mc_pz.clear();

    _mc_vx.clear();
    _mc_vy.clear();
    _mc_vz.clear();

    _mc_endx.clear();
    _mc_endy.clear();
    _mc_endz.clear();

    _mc_completeness.clear();
    _mc_purity.clear();

    _true_pt = 0;
    _true_pt_visible = 0;
    _true_p = 0;
    _true_p_visible = 0;

    _true_e_visible = 0;
}

void StrangenessAnalysis_tool::SaveTruth(art::Event const &e)
{
    auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    auto const &mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    auto const& mcflux_h = e.getValidHandle<std::vector<simb::MCFlux>>(fMCFluxproducer);

    // reference: http://www.hep.utexas.edu/~zarko/wwwgnumi/v19/
    /*
        Decay mode that produced neutrino:

        1  K0L -> nue pi- e+
        2  K0L -> nuebar pi+ e-
        3  K0L -> numu pi- mu+
        4  K0L -> numubar pi+ mu-
        5  K+  -> numu mu+
        6  K+  -> nue pi0 e+
        7  K+  -> numu pi0 mu+
        8  K-  -> numubar mu-
        9  K-  -> nuebar pi0 e-
        10  K-  -> numubar pi0 mu-
        11  mu+ -> numubar nue e+
        12  mu- -> numu nuebar e-
        13  pi+ -> numu mu+
        14  pi- -> numubar mu-
    */
    auto flux = mcflux_h->at(0);
    _nu_parent_pdg = flux.fptype;
    _nu_hadron_pdg = flux.ftptype;
    _nu_decay_mode = flux.fndecay;

    auto mct = mct_h->at(0);
    auto neutrino = mct.GetNeutrino();
    auto nu = neutrino.Nu();

    _ccnc = neutrino.CCNC();
    _interaction = neutrino.Mode();
    _nu_pdg = nu.PdgCode();
    _nu_e = nu.Trajectory().E(0);
    _nu_l = flux.fdk2gen + flux.fgen2vtx;
    std::cout << "total length = " << flux.fdk2gen << " + " << flux.fgen2vtx << " = " << _nu_l << std::endl;
    _lep_e = neutrino.Lepton().E();

    _true_nu_vtx_t = nu.T();
    _true_nu_vtx_x = nu.Vx();
    _true_nu_vtx_y = nu.Vy();
    _true_nu_vtx_z = nu.Vz();

    float _true_nu_vtx_sce[3];
    common::True2RecoMappingXYZ(_true_nu_vtx_t, _true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z, _true_nu_vtx_sce);

    _true_nu_vtx_sce_x = _true_nu_vtx_sce[0];
    _true_nu_vtx_sce_y = _true_nu_vtx_sce[1];
    _true_nu_vtx_sce_z = _true_nu_vtx_sce[2];

    _theta = neutrino.Theta();
    _nu_pt = neutrino.Pt();

    double vtx[3] = {_true_nu_vtx_x, _true_nu_vtx_y, _true_nu_vtx_z};
    _isVtxInFiducial = common::isFiducial(vtx, fFidvolXstart, fFidvolYstart, fFidvolZstart, fFidvolXend, fFidvolYend, fFidvolZend);

    _nelec = 0;
    _nmuon = 0;
    _nproton = 0;
    _npion = 0;
    _nneutron = 0;

    float muonMomentum = 0;

    TLorentzVector total_p;
    TLorentzVector total_p_visible;

    size_t npart = mct.NParticles();
    for (size_t i = 0; i < npart; i++)
    {
        auto const &part = mct.GetParticle(i);

        if (part.PdgCode() == kaon_short->PdgCode()) { 
            _nkshrt += 1; 
            _kshrt_e = part.Momentum(0).E(); 
        }

        if (part.PdgCode() == lambda->PdgCode()) {
            _nlambda += 1;
            _lambda_e = part.Momentum(0).E();
        }

        if (part.StatusCode() != 1)
        {
            continue;
        }

        total_p += part.Momentum(0);

        if ((std::abs(part.PdgCode()) == muon->PdgCode()) && (part.StatusCode() == 1))
        {
            muonMomentum = part.Momentum(0).E();

            if (part.Momentum(0).E() - muon->Mass() > fMuonThreshold)
            {
                _true_e_visible += part.Momentum(0).E() - muon->Mass();
                total_p_visible += part.Momentum(0);
                _nmuon += 1;
            }

            if (part.Momentum(0).E() > _muon_e)
                _muon_e = part.Momentum(0).E();
        } 

        // if electron
        else if ((std::abs(part.PdgCode()) == electron->PdgCode()) and (part.StatusCode() == 1))
        {
            if (part.Momentum(0).E() - electron->Mass() > fElectronThreshold)
            {
                _nelec += 1;
                total_p_visible += part.Momentum(0);
                _true_e_visible += part.Momentum(0).E() - electron->Mass();
            }
            if (part.Momentum(0).E() > _elec_e)
                _elec_e = part.Momentum(0).E();
        } // if electron

        // if proton
        else if ((part.PdgCode() == proton->PdgCode()) and (part.StatusCode() == 1))
        {
            if (part.Momentum(0).E() - proton->Mass() > fProtonThreshold)
            {
                total_p_visible += part.Momentum(0);
                _true_e_visible += part.Momentum(0).E() - proton->Mass();
                _nproton += 1;
            }
            if (part.Momentum(0).E() > _proton_e)
                _proton_e = part.Momentum(0).E();
        } // if proton

        // if neutron
        else if ((part.PdgCode() == neutron->PdgCode()) and (part.StatusCode() == 1))
        {
            _nneutron += 1;
        }

        // if pion
        else if ((std::abs(part.PdgCode()) == pion->PdgCode()) and (part.StatusCode() == 1))
        {
            if (part.Momentum(0).E() - pion->Mass() > fPionThreshold)
            {
                _npion += 1;
                total_p_visible += part.Momentum(0);
                _true_e_visible += part.Momentum(0).E() - pion->Mass();
            }
            if (part.Momentum(0).E() > _pion_e)
                _pion_e = part.Momentum(0).E();
        } // if pion    

        else
        {
            TParticlePDG *particle_pdg = TDatabasePDG::Instance()->GetParticle(part.PdgCode());
            if (particle_pdg != NULL) // PDG codes corresponding to ions e.g. 2000000101 are not in the database
                _true_e_visible += part.Momentum(0).E() - particle_pdg->Mass();
        }

    } 

    _true_pt = total_p.Perp();
    _true_pt_visible = total_p_visible.Perp();
    _true_p = total_p.Mag();
    _true_p_visible = total_p_visible.Mag();

    for (size_t p = 0; p < mcp_h->size(); p++)
    {
        auto mcp = mcp_h->at(p);
        if (!(mcp.Process() == "primary" && mcp.StatusCode() == 1))
        {
            continue;
        }

        _mc_E.push_back(mcp.E());

        _mc_pdg.push_back(mcp.PdgCode());

        _mc_px.push_back(mcp.Px());
        _mc_py.push_back(mcp.Py());
        _mc_pz.push_back(mcp.Pz());

        _mc_vx.push_back(mcp.Vx());
        _mc_vy.push_back(mcp.Vy());
        _mc_vz.push_back(mcp.Vz());

        _mc_endx.push_back(mcp.EndX());
        _mc_endy.push_back(mcp.EndY());
        _mc_endz.push_back(mcp.EndZ());

        _mc_completeness.push_back(std::numeric_limits<float>::lowest());
        _mc_purity.push_back(std::numeric_limits<float>::lowest());
    }

    // find kshrt -> pi+ pi-
    _endkshrtprocess = "";
    _endkshrthadronic = 0;
    bool containedkshrt = false;
    float kshrtendpointx = 0;

    int kshrt_id = -1;

    for (size_t p = 0; p < mcp_h->size(); p++)
    {
        auto mcp = mcp_h->at(p);

        if (mcp.StatusCode() != 1)
        {
            continue;
        }

        if ((mcp.Momentum(0).E() - muonMomentum) < 0.0001)
        {
            kshrt_id = mcp.TrackId();
          
            art::ServiceHandle<geo::Geometry> geo;
            geo::TPCGeo const &thisTPC = geo->TPC();
            geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
            if ((mcp.EndPosition().X() > theTpcGeo.MinX()) && (mcp.EndPosition().X() < theTpcGeo.MaxX()) &&
                (mcp.EndPosition().Y() > theTpcGeo.MinY()) && (mcp.EndPosition().Y() < theTpcGeo.MaxY()) &&
                (mcp.EndPosition().Z() > theTpcGeo.MinZ()) && (mcp.EndPosition().Z() < theTpcGeo.MaxZ()))
            {
                _endkshrtprocess = mcp.EndProcess();
                containedkshrt = true;
                kshrtendpointx = mcp.EndPosition().X();
            }

            break;
        } 
    }  

    bool foundpositivepion; 
    bool foundnegativepion;

    if (containedkshrt)
    {
        for (size_t p = 0; p < mcp_h->size(); p++)
        {
            auto mcp = mcp_h->at(p);

            if (fabs(mcp.PdgCode()) != pion->PdgCode())
                continue;

            if (mcp.Mother() == kshrt_id)
            {
                if ((mcp.Vx() - kshrtendpointx) < 0.0001)
                {
                    if (mcp.PdgCode() > 0)
                    {
                        foundpositivepion = true;
                    }
                    else
                    {
                        foundnegativepion = true;
                    }
                } 
            }   
        }     
    }       

    if (foundpositivepion && foundnegativepion)
    {

    }

    return;
}

DEFINE_ART_CLASS_TOOL(StrangenessAnalysis_tool)
} 

#endif