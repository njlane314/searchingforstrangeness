#ifndef ANALYSIS_SIGNAL_CXX
#define ANALYSIS_SIGNAL_CXX

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraSliceIdHelper.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"
#include "Common/GeometryUtils.h"
#include "Common/SpaceChargeCorrections.h"
#include "Common/TruthContainment.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

#include <cmath>
#include <unordered_map>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <array>

namespace analysis {

class SignalAnalysis : public AnalysisToolBase {
public:
    explicit SignalAnalysis(fhicl::ParameterSet const& p);
    ~SignalAnalysis() override = default;

    void configure(const fhicl::ParameterSet& pset) override;
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;
    void analyseEvent(const art::Event& event, bool is_data) override;
    void analyseSlice(const art::Event& event,
                      std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                      bool is_data, bool is_selected) override;

private:
    art::InputTag fMCTproducer;
    art::InputTag fMCPproducer;
    art::InputTag fMCRproducer;
    art::InputTag fCLSproducer;
    art::InputTag fHITproducer;
    art::InputTag fBKTproducer;

    float fFidvolXstart, fFidvolXend;
    float fFidvolYstart, fFidvolYend;
    float fFidvolZstart, fFidvolZend;

    bool  fRequireNuMu;
    bool  fRequireCC;
    bool  fAllowAntiLambda;
    bool  fRequireSpecificDecay;
    bool  fRequirePrimaryLambda;
    bool  fChooseClosestToNuVertex;
    int   fMaxAncestorSteps;
    float fProtonKEThreshold;
    float fPionKEThreshold;

    bool  fAcceptAntiNu;
    float fPRMinCompleteness;
    float fPRMinPurity;

    int   _nu_pdg;
    int   _ccnc;
    int   _interaction_mode;
    int   _interaction_type;
    float _nu_E;
    float _nu_theta;
    float _nu_pt;
    float _nu_vx, _nu_vy, _nu_vz;
    float _nu_t;
    float _nu_sce_vx, _nu_sce_vy, _nu_sce_vz;
    bool  _nu_vtx_in_fid;

    float _true_nu_dir_x, _true_nu_dir_y, _true_nu_dir_z;

    bool  _is_nu_mu_cc;
    bool  _has_lambda_to_ppi;
    bool  _is_signal_event;

    int _n_lambda_candidates;
    std::vector<int>   _lambda_trackid;
    std::vector<int>   _lambda_pdg;
    std::vector<int>   _lambda_mother_trackid;
    std::vector<int>   _lambda_is_primary;
    std::vector<float> _lambda_mass;
    std::vector<float> _lambda_energy;
    std::vector<float> _lambda_px, _lambda_py, _lambda_pz;
    std::vector<float> _lambda_vx, _lambda_vy, _lambda_vz;
    std::vector<float> _lambda_endx, _lambda_endy, _lambda_endz;
    std::vector<float> _lambda_len;
    std::vector<float> _lambda_ct;
    std::vector<float> _lambda_dist_to_nu_vtx;
    std::vector<int>   _lambda_decay_in_fid;

    std::vector<unsigned int> _lambda_origin_bits;
    std::vector<int>          _lambda_from_sigma0;
    std::vector<int>          _lambda_heavy_feed;
    std::vector<int>          _lambda_parent_pdg;
    std::vector<int>          _lambda_grandparent_pdg;

    std::vector<int>   _lambda_proton_trackid;
    std::vector<int>   _lambda_pion_trackid;
    std::vector<float> _proton_energy, _pion_energy;
    std::vector<float> _proton_px, _proton_py, _proton_pz;
    std::vector<float> _pion_px,   _pion_py,   _pion_pz;
    std::vector<float> _proton_vx, _proton_vy, _proton_vz;
    std::vector<float> _pion_vx,   _pion_vy,   _pion_vz;
    std::vector<float> _ppi_invariant_mass;
    std::vector<float> _ppi_opening_angle;

    int   _sel_index;
    int   _sel_lambda_trackid;
    int   _sel_lambda_pdg;
    float _sel_lambda_mass;
    float _sel_lambda_energy;
    float _sel_lambda_px, _sel_lambda_py, _sel_lambda_pz;
    float _sel_lambda_vx, _sel_lambda_vy, _sel_lambda_vz;
    float _sel_lambda_endx, _sel_lambda_endy, _sel_lambda_endz;
    float _sel_lambda_len, _sel_lambda_ct;
    float _sel_lambda_dist_to_nu_vtx;
    bool  _sel_lambda_decay_in_fid;

    bool  _sel_lambda_from_sigma0;
    bool  _sel_lambda_heavy_feed;
    int   _sel_lambda_parent_pdg;
    int   _sel_lambda_grandparent_pdg;

    int   _sel_proton_trackid, _sel_pion_trackid;
    float _sel_proton_energy, _sel_pion_energy;
    float _sel_proton_px, _sel_proton_py, _sel_proton_pz;
    float _sel_pion_px,   _sel_pion_py,   _sel_pion_pz;
    float _sel_ppi_invariant_mass, _sel_ppi_opening_angle;

    float _sel_proton_purity, _sel_proton_completeness;
    float _sel_pion_purity,   _sel_pion_completeness;

    int   _mu_truth_trackid;
    int   _mu_truth_pdg;
    float _mu_truth_E;
    float _mu_truth_px, _mu_truth_py, _mu_truth_pz;
    float _mu_truth_pt, _mu_truth_theta;
    float _mu_truth_vx, _mu_truth_vy, _mu_truth_vz;

    bool  _pr_eligible_event;
    bool  _pr_valid_assignment;
    bool  _pr_pass_event;
    int   _pr_n_pfps_in_slice;
    int   _pr_best_sum_nshared;

    int _pr_mu_pfp_index, _pr_p_pfp_index, _pr_pi_pfp_index;

    int   _pr_mu_nshared, _pr_mu_nhits_O, _pr_mu_nhits_T;
    int   _pr_p_nshared,  _pr_p_nhits_O,  _pr_p_nhits_T;
    int   _pr_pi_nshared, _pr_pi_nhits_O, _pr_pi_nhits_T;
    float _pr_mu_purity, _pr_mu_completeness;
    float _pr_p_purity,  _pr_p_completeness;
    float _pr_pi_purity, _pr_pi_completeness;

    int  _n_sigma0_truth;
    int  _n_kshort_truth;
    bool _has_sigma0_truth;
    bool _has_kshort_truth;
    int  _n_lambda_from_heavy;

    template <class T> static T nan() { return std::numeric_limits<T>::quiet_NaN(); }
    static float ThreeDistance(float x1,float y1,float z1,float x2,float y2,float z2) {
        const float dx = x1-x2, dy = y1-y2, dz = z1-z2;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    struct DecayMatch {
        bool ok = false;
        int  p_trkid  = -1;
        int  pi_trkid = -1;
    };

    DecayMatch MatchLambdaToPPi(const art::Ptr<simb::MCParticle>& lam,
                                const std::map<int, art::Ptr<simb::MCParticle>>& mp,
                                bool allowAnti) const;

    void SelectBestCandidate();

    void   FindTruthMuon(const art::ValidHandle<std::vector<simb::MCParticle>>& mcp_h);
    size_t CountTruthHitsInSlice(int trackid,
                                 const std::vector<art::Ptr<recob::Hit>>& hits,
                                 std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& assocMCPart) const;
};

SignalAnalysis::SignalAnalysis(fhicl::ParameterSet const& p) { this->configure(p); }

void SignalAnalysis::configure(const fhicl::ParameterSet& p) {
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fMCRproducer = p.get<art::InputTag>("MCRproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fBKTproducer = p.get<art::InputTag>("BKTproducer");

    fFidvolXstart = p.get<double>("fidvolXstart", 10);
    fFidvolXend   = p.get<double>("fidvolXend",   10);
    fFidvolYstart = p.get<double>("fidvolYstart", 15);
    fFidvolYend   = p.get<double>("fidvolYend",   15);
    fFidvolZstart = p.get<double>("fidvolZstart", 10);
    fFidvolZend   = p.get<double>("fidvolZend",   50);

    fRequireNuMu             = p.get<bool>("RequireNuMu", true);
    fRequireCC               = p.get<bool>("RequireCC",   true);
    fAllowAntiLambda         = p.get<bool>("AllowAntiLambda", false);
    fRequireSpecificDecay    = p.get<bool>("RequireSpecificDecay", true);
    fRequirePrimaryLambda    = p.get<bool>("RequirePrimaryLambda", false);
    fChooseClosestToNuVertex = p.get<bool>("ChooseClosestToNuVertex", true);
    fMaxAncestorSteps        = p.get<int> ("MaxAncestorSteps", 16);
    fProtonKEThreshold       = p.get<float>("ProtonKEThreshold",  0.0f);
    fPionKEThreshold         = p.get<float>("PionKEThreshold",    0.0f);

    fAcceptAntiNu            = p.get<bool>("AcceptAntiNu", true);
    fPRMinCompleteness       = p.get<float>("PatternRecoMinCompleteness", 0.1f);
    fPRMinPurity             = p.get<float>("PatternRecoMinPurity",       0.5f);
}

void SignalAnalysis::setBranches(TTree* t) {
    t->Branch("nu_pdg",           &_nu_pdg,           "nu_pdg/I");
    t->Branch("ccnc",             &_ccnc,             "ccnc/I");
    t->Branch("interaction_mode", &_interaction_mode, "interaction_mode/I");
    t->Branch("interaction_type", &_interaction_type, "interaction_type/I");
    t->Branch("nu_E",             &_nu_E,             "nu_E/F");
    t->Branch("nu_theta",         &_nu_theta,         "nu_theta/F");
    t->Branch("nu_pt",            &_nu_pt,            "nu_pt/F");
    t->Branch("nu_vx",            &_nu_vx,            "nu_vx/F");
    t->Branch("nu_vy",            &_nu_vy,            "nu_vy/F");
    t->Branch("nu_vz",            &_nu_vz,            "nu_vz/F");
    t->Branch("nu_t",             &_nu_t,             "nu_t/F");
    t->Branch("nu_sce_vx",        &_nu_sce_vx,        "nu_sce_vx/F");
    t->Branch("nu_sce_vy",        &_nu_sce_vy,        "nu_sce_vy/F");
    t->Branch("nu_sce_vz",        &_nu_sce_vz,        "nu_sce_vz/F");
    t->Branch("nu_vtx_in_fid",    &_nu_vtx_in_fid,    "nu_vtx_in_fid/O");

    t->Branch("true_nu_dir_x", &_true_nu_dir_x, "true_nu_dir_x/F");
    t->Branch("true_nu_dir_y", &_true_nu_dir_y, "true_nu_dir_y/F");
    t->Branch("true_nu_dir_z", &_true_nu_dir_z, "true_nu_dir_z/F");

    t->Branch("is_nu_mu_cc",      &_is_nu_mu_cc,      "is_nu_mu_cc/O");
    t->Branch("has_lambda_to_ppi",&_has_lambda_to_ppi,"has_lambda_to_ppi/O");
    t->Branch("is_signal_event",  &_is_signal_event,  "is_signal_event/O");

    t->Branch("n_lambda_candidates", &_n_lambda_candidates, "n_lambda_candidates/I");
    t->Branch("lambda_trackid",       "std::vector<int>",   &_lambda_trackid);
    t->Branch("lambda_pdg",           "std::vector<int>",   &_lambda_pdg);
    t->Branch("lambda_mother_trackid","std::vector<int>",   &_lambda_mother_trackid);
    t->Branch("lambda_is_primary",    "std::vector<int>",   &_lambda_is_primary);
    t->Branch("lambda_mass",          "std::vector<float>", &_lambda_mass);
    t->Branch("lambda_energy",        "std::vector<float>", &_lambda_energy);
    t->Branch("lambda_px",            "std::vector<float>", &_lambda_px);
    t->Branch("lambda_py",            "std::vector<float>", &_lambda_py);
    t->Branch("lambda_pz",            "std::vector<float>", &_lambda_pz);
    t->Branch("lambda_vx",            "std::vector<float>", &_lambda_vx);
    t->Branch("lambda_vy",            "std::vector<float>", &_lambda_vy);
    t->Branch("lambda_vz",            "std::vector<float>", &_lambda_vz);
    t->Branch("lambda_endx",          "std::vector<float>", &_lambda_endx);
    t->Branch("lambda_endy",          "std::vector<float>", &_lambda_endy);
    t->Branch("lambda_endz",          "std::vector<float>", &_lambda_endz);
    t->Branch("lambda_len",           "std::vector<float>", &_lambda_len);
    t->Branch("lambda_ct",            "std::vector<float>", &_lambda_ct);
    t->Branch("lambda_dist_to_nu_vtx","std::vector<float>", &_lambda_dist_to_nu_vtx);
    t->Branch("lambda_decay_in_fid",  "std::vector<int>",   &_lambda_decay_in_fid);

    t->Branch("lambda_origin_bits",     "std::vector<unsigned int>", &_lambda_origin_bits);
    t->Branch("lambda_from_sigma0",     "std::vector<int>",           &_lambda_from_sigma0);
    t->Branch("lambda_heavy_feed",      "std::vector<int>",           &_lambda_heavy_feed);
    t->Branch("lambda_parent_pdg",      "std::vector<int>",           &_lambda_parent_pdg);
    t->Branch("lambda_grandparent_pdg", "std::vector<int>",           &_lambda_grandparent_pdg);

    t->Branch("lambda_proton_trackid","std::vector<int>",   &_lambda_proton_trackid);
    t->Branch("lambda_pion_trackid",  "std::vector<int>",   &_lambda_pion_trackid);
    t->Branch("proton_energy",        "std::vector<float>", &_proton_energy);
    t->Branch("proton_px",            "std::vector<float>", &_proton_px);
    t->Branch("proton_py",            "std::vector<float>", &_proton_py);
    t->Branch("proton_pz",            "std::vector<float>", &_proton_pz);
    t->Branch("proton_vx",            "std::vector<float>", &_proton_vx);
    t->Branch("proton_vy",            "std::vector<float>", &_proton_vy);
    t->Branch("proton_vz",            "std::vector<float>", &_proton_vz);
    t->Branch("pion_energy",          "std::vector<float>", &_pion_energy);
    t->Branch("pion_px",              "std::vector<float>", &_pion_px);
    t->Branch("pion_py",              "std::vector<float>", &_pion_py);
    t->Branch("pion_pz",              "std::vector<float>", &_pion_pz);
    t->Branch("pion_vx",              "std::vector<float>", &_pion_vx);
    t->Branch("pion_vy",              "std::vector<float>", &_pion_vy);
    t->Branch("pion_vz",              "std::vector<float>", &_pion_vz);
    t->Branch("ppi_invariant_mass",   "std::vector<float>", &_ppi_invariant_mass);
    t->Branch("ppi_opening_angle",    "std::vector<float>", &_ppi_opening_angle);

    t->Branch("sel_index",            &_sel_index,            "sel_index/I");
    t->Branch("sel_lambda_trackid",   &_sel_lambda_trackid,   "sel_lambda_trackid/I");
    t->Branch("sel_lambda_pdg",       &_sel_lambda_pdg,       "sel_lambda_pdg/I");
    t->Branch("sel_lambda_mass",      &_sel_lambda_mass,      "sel_lambda_mass/F");
    t->Branch("sel_lambda_energy",    &_sel_lambda_energy,    "sel_lambda_energy/F");
    t->Branch("sel_lambda_px",        &_sel_lambda_px,        "sel_lambda_px/F");
    t->Branch("sel_lambda_py",        &_sel_lambda_py,        "sel_lambda_py/F");
    t->Branch("sel_lambda_pz",        &_sel_lambda_pz,        "sel_lambda_pz/F");
    t->Branch("sel_lambda_vx",        &_sel_lambda_vx,        "sel_lambda_vx/F");
    t->Branch("sel_lambda_vy",        &_sel_lambda_vy,        "sel_lambda_vy/F");
    t->Branch("sel_lambda_vz",        &_sel_lambda_vz,        "sel_lambda_vz/F");
    t->Branch("sel_lambda_endx",      &_sel_lambda_endx,      "sel_lambda_endx/F");
    t->Branch("sel_lambda_endy",      &_sel_lambda_endy,      "sel_lambda_endy/F");
    t->Branch("sel_lambda_endz",      &_sel_lambda_endz,      "sel_lambda_endz/F");
    t->Branch("sel_lambda_len",       &_sel_lambda_len,       "sel_lambda_len/F");
    t->Branch("sel_lambda_ct",        &_sel_lambda_ct,        "sel_lambda_ct/F");
    t->Branch("sel_lambda_dist_to_nu_vtx",&_sel_lambda_dist_to_nu_vtx,"sel_lambda_dist_to_nu_vtx/F");
    t->Branch("sel_lambda_decay_in_fid",&_sel_lambda_decay_in_fid,"sel_lambda_decay_in_fid/O");

    t->Branch("sel_lambda_from_sigma0",&_sel_lambda_from_sigma0,"sel_lambda_from_sigma0/O");
    t->Branch("sel_lambda_heavy_feed", &_sel_lambda_heavy_feed, "sel_lambda_heavy_feed/O");
    t->Branch("sel_lambda_parent_pdg", &_sel_lambda_parent_pdg, "sel_lambda_parent_pdg/I");
    t->Branch("sel_lambda_grandparent_pdg",&_sel_lambda_grandparent_pdg,"sel_lambda_grandparent_pdg/I");

    t->Branch("sel_proton_trackid",   &_sel_proton_trackid,   "sel_proton_trackid/I");
    t->Branch("sel_pion_trackid",     &_sel_pion_trackid,     "sel_pion_trackid/I");
    t->Branch("sel_proton_energy",    &_sel_proton_energy,    "sel_proton_energy/F");
    t->Branch("sel_pion_energy",      &_sel_pion_energy,      "sel_pion_energy/F");
    t->Branch("sel_proton_px",        &_sel_proton_px,        "sel_proton_px/F");
    t->Branch("sel_proton_py",        &_sel_proton_py,        "sel_proton_py/F");
    t->Branch("sel_proton_pz",        &_sel_proton_pz,        "sel_proton_pz/F");
    t->Branch("sel_pion_px",          &_sel_pion_px,          "sel_pion_px/F");
    t->Branch("sel_pion_py",          &_sel_pion_py,          "sel_pion_py/F");
    t->Branch("sel_pion_pz",          &_sel_pion_pz,          "sel_pion_pz/F");
    t->Branch("sel_ppi_invariant_mass",&_sel_ppi_invariant_mass,"sel_ppi_invariant_mass/F");
    t->Branch("sel_ppi_opening_angle",&_sel_ppi_opening_angle,"sel_ppi_opening_angle/F");

    t->Branch("sel_proton_purity",     &_sel_proton_purity,     "sel_proton_purity/F");
    t->Branch("sel_proton_completeness",&_sel_proton_completeness,"sel_proton_completeness/F");
    t->Branch("sel_pion_purity",       &_sel_pion_purity,       "sel_pion_purity/F");
    t->Branch("sel_pion_completeness", &_sel_pion_completeness, "sel_pion_completeness/F");

    t->Branch("mu_truth_trackid", &_mu_truth_trackid, "mu_truth_trackid/I");
    t->Branch("mu_truth_pdg",     &_mu_truth_pdg,     "mu_truth_pdg/I");
    t->Branch("mu_truth_E",       &_mu_truth_E,       "mu_truth_E/F");
    t->Branch("mu_truth_px",      &_mu_truth_px,      "mu_truth_px/F");
    t->Branch("mu_truth_py",      &_mu_truth_py,      "mu_truth_py/F");
    t->Branch("mu_truth_pz",      &_mu_truth_pz,      "mu_truth_pz/F");
    t->Branch("mu_truth_pt",      &_mu_truth_pt,      "mu_truth_pt/F");
    t->Branch("mu_truth_theta",   &_mu_truth_theta,   "mu_truth_theta/F");
    t->Branch("mu_truth_vx",      &_mu_truth_vx,      "mu_truth_vx/F");
    t->Branch("mu_truth_vy",      &_mu_truth_vy,      "mu_truth_vy/F");
    t->Branch("mu_truth_vz",      &_mu_truth_vz,      "mu_truth_vz/F");

    t->Branch("pr_eligible_event",    &_pr_eligible_event,    "pr_eligible_event/O");
    t->Branch("pr_valid_assignment",  &_pr_valid_assignment,  "pr_valid_assignment/O");
    t->Branch("pr_pass_event",        &_pr_pass_event,        "pr_pass_event/O");
    t->Branch("pr_n_pfps_in_slice",   &_pr_n_pfps_in_slice,   "pr_n_pfps_in_slice/I");
    t->Branch("pr_best_sum_nshared",  &_pr_best_sum_nshared,  "pr_best_sum_nshared/I");

    t->Branch("pr_mu_pfp_index", &_pr_mu_pfp_index, "pr_mu_pfp_index/I");
    t->Branch("pr_p_pfp_index",  &_pr_p_pfp_index,  "pr_p_pfp_index/I");
    t->Branch("pr_pi_pfp_index", &_pr_pi_pfp_index, "pr_pi_pfp_index/I");

    t->Branch("pr_mu_nshared",   &_pr_mu_nshared,   "pr_mu_nshared/I");
    t->Branch("pr_mu_nhits_O",   &_pr_mu_nhits_O,   "pr_mu_nhits_O/I");
    t->Branch("pr_mu_nhits_T",   &_pr_mu_nhits_T,   "pr_mu_nhits_T/I");
    t->Branch("pr_mu_purity",    &_pr_mu_purity,    "pr_mu_purity/F");
    t->Branch("pr_mu_completeness",&_pr_mu_completeness,"pr_mu_completeness/F");

    t->Branch("pr_p_nshared",    &_pr_p_nshared,    "pr_p_nshared/I");
    t->Branch("pr_p_nhits_O",    &_pr_p_nhits_O,    "pr_p_nhits_O/I");
    t->Branch("pr_p_nhits_T",    &_pr_p_nhits_T,    "pr_p_nhits_T/I");
    t->Branch("pr_p_purity",     &_pr_p_purity,     "pr_p_purity/F");
    t->Branch("pr_p_completeness",&_pr_p_completeness,"pr_p_completeness/F");

    t->Branch("pr_pi_nshared",   &_pr_pi_nshared,   "pr_pi_nshared/I");
    t->Branch("pr_pi_nhits_O",   &_pr_pi_nhits_O,   "pr_pi_nhits_O/I");
    t->Branch("pr_pi_nhits_T",   &_pr_pi_nhits_T,   "pr_pi_nhits_T/I");
    t->Branch("pr_pi_purity",    &_pr_pi_purity,    "pr_pi_purity/F");
    t->Branch("pr_pi_completeness",&_pr_pi_completeness,"pr_pi_completeness/F");

    t->Branch("n_sigma0_truth",  &_n_sigma0_truth,  "n_sigma0_truth/I");
    t->Branch("has_sigma0_truth",&_has_sigma0_truth,"has_sigma0_truth/O");
    t->Branch("n_kshort_truth",  &_n_kshort_truth,  "n_kshort_truth/I");
    t->Branch("has_kshort_truth",&_has_kshort_truth,"has_kshort_truth/O");
    t->Branch("n_lambda_from_heavy",&_n_lambda_from_heavy,"n_lambda_from_heavy/I");
}

void SignalAnalysis::resetTTree(TTree*) {
    _nu_pdg = 0; _ccnc = -1; _interaction_mode = -1; _interaction_type = -1;
    _nu_E = nan<float>(); _nu_theta = nan<float>(); _nu_pt = nan<float>();
    _nu_vx = nan<float>(); _nu_vy = nan<float>(); _nu_vz = nan<float>(); _nu_t = nan<float>();
    _nu_sce_vx = nan<float>(); _nu_sce_vy = nan<float>(); _nu_sce_vz = nan<float>();
    _nu_vtx_in_fid = false;
    _true_nu_dir_x = nan<float>(); _true_nu_dir_y = nan<float>(); _true_nu_dir_z = nan<float>();

    _is_nu_mu_cc = false; _has_lambda_to_ppi = false; _is_signal_event = false;

    _n_lambda_candidates = 0;
    _lambda_trackid.clear(); _lambda_pdg.clear(); _lambda_mother_trackid.clear();
    _lambda_is_primary.clear(); _lambda_mass.clear(); _lambda_energy.clear();
    _lambda_px.clear(); _lambda_py.clear(); _lambda_pz.clear();
    _lambda_vx.clear(); _lambda_vy.clear(); _lambda_vz.clear();
    _lambda_endx.clear(); _lambda_endy.clear(); _lambda_endz.clear();
    _lambda_len.clear(); _lambda_ct.clear(); _lambda_dist_to_nu_vtx.clear();
    _lambda_decay_in_fid.clear();

    _lambda_origin_bits.clear();
    _lambda_from_sigma0.clear();
    _lambda_heavy_feed.clear();
    _lambda_parent_pdg.clear();
    _lambda_grandparent_pdg.clear();

    _lambda_proton_trackid.clear(); _lambda_pion_trackid.clear();
    _proton_energy.clear(); _pion_energy.clear();
    _proton_px.clear(); _proton_py.clear(); _proton_pz.clear();
    _pion_px.clear();   _pion_py.clear();   _pion_pz.clear();
    _proton_vx.clear(); _proton_vy.clear(); _proton_vz.clear();
    _pion_vx.clear();   _pion_vy.clear();   _pion_vz.clear();
    _ppi_invariant_mass.clear(); _ppi_opening_angle.clear();

    _sel_index = -1;
    _sel_lambda_trackid = -1; _sel_lambda_pdg = 0;
    _sel_lambda_mass = nan<float>(); _sel_lambda_energy = nan<float>();
    _sel_lambda_px = nan<float>(); _sel_lambda_py = nan<float>(); _sel_lambda_pz = nan<float>();
    _sel_lambda_vx = nan<float>(); _sel_lambda_vy = nan<float>(); _sel_lambda_vz = nan<float>();
    _sel_lambda_endx = nan<float>(); _sel_lambda_endy = nan<float>(); _sel_lambda_endz = nan<float>();
    _sel_lambda_len = nan<float>(); _sel_lambda_ct = nan<float>(); _sel_lambda_dist_to_nu_vtx = nan<float>();
    _sel_lambda_decay_in_fid = false;

    _sel_lambda_from_sigma0 = false;
    _sel_lambda_heavy_feed  = false;
    _sel_lambda_parent_pdg = 0;
    _sel_lambda_grandparent_pdg = 0;

    _sel_proton_trackid = -1; _sel_pion_trackid = -1;
    _sel_proton_energy = nan<float>(); _sel_pion_energy = nan<float>();
    _sel_proton_px = nan<float>(); _sel_proton_py = nan<float>(); _sel_proton_pz = nan<float>();
    _sel_pion_px   = nan<float>(); _sel_pion_py   = nan<float>(); _sel_pion_pz   = nan<float>();
    _sel_ppi_invariant_mass = nan<float>(); _sel_ppi_opening_angle = nan<float>();

    _sel_proton_purity = 0.f; _sel_proton_completeness = 0.f;
    _sel_pion_purity   = 0.f; _sel_pion_completeness   = 0.f;

    _mu_truth_trackid = -1; _mu_truth_pdg = 0;
    _mu_truth_E = nan<float>(); _mu_truth_px = nan<float>(); _mu_truth_py = nan<float>(); _mu_truth_pz = nan<float>();
    _mu_truth_pt = nan<float>(); _mu_truth_theta = nan<float>();
    _mu_truth_vx = nan<float>(); _mu_truth_vy = nan<float>(); _mu_truth_vz = nan<float>();

    _pr_eligible_event = false; _pr_valid_assignment = false; _pr_pass_event = false;
    _pr_n_pfps_in_slice = 0; _pr_best_sum_nshared = 0;
    _pr_mu_pfp_index = -1; _pr_p_pfp_index = -1; _pr_pi_pfp_index = -1;

    _pr_mu_nshared = 0; _pr_mu_nhits_O = 0; _pr_mu_nhits_T = 0;
    _pr_p_nshared  = 0; _pr_p_nhits_O  = 0; _pr_p_nhits_T  = 0;
    _pr_pi_nshared = 0; _pr_pi_nhits_O = 0; _pr_pi_nhits_T = 0;

    _pr_mu_purity = 0.f; _pr_mu_completeness = 0.f;
    _pr_p_purity  = 0.f; _pr_p_completeness  = 0.f;
    _pr_pi_purity = 0.f; _pr_pi_completeness = 0.f;

    _n_sigma0_truth = 0;
    _n_kshort_truth = 0;
    _has_sigma0_truth = false;
    _has_kshort_truth = false;
    _n_lambda_from_heavy = 0;
}

SignalAnalysis::DecayMatch
SignalAnalysis::MatchLambdaToPPi(const art::Ptr<simb::MCParticle>& lam,
                                       const std::map<int, art::Ptr<simb::MCParticle>>& mp,
                                       bool allowAnti) const {
    DecayMatch ret;
    const int lam_pdg = lam->PdgCode();
    if (std::abs(lam_pdg) != 3122) return ret;
    if (!allowAnti && lam_pdg != 3122) return ret;

    const int want_p  = (lam_pdg > 0) ? 2212 : -2212;
    const int want_pi = (lam_pdg > 0) ? -211 :  211;

    int p_tid  = -1;
    int pi_tid = -1;
    for (int i = 0; i < lam->NumberDaughters(); ++i) {
        const int d_tid = lam->Daughter(i);
        auto it = mp.find(d_tid);
        if (it == mp.end()) continue;
        const auto& d = it->second;
        if (d->PdgCode() == want_p)  p_tid  = d->TrackId();
        if (d->PdgCode() == want_pi) pi_tid = d->TrackId();
    }
    if (p_tid >= 0 && pi_tid >= 0) {
        ret.ok = true;
        ret.p_trkid  = p_tid;
        ret.pi_trkid = pi_tid;
    }
    return ret;
}

void SignalAnalysis::SelectBestCandidate() {
    _sel_index = -1;
    if (_n_lambda_candidates <= 0) return;

    if (fChooseClosestToNuVertex) {
        auto it = std::min_element(_lambda_dist_to_nu_vtx.begin(), _lambda_dist_to_nu_vtx.end());
        _sel_index = static_cast<int>(std::distance(_lambda_dist_to_nu_vtx.begin(), it));
    } else {
        _sel_index = 0;
    }

    const int i = _sel_index;
    _sel_lambda_trackid = _lambda_trackid[i];
    _sel_lambda_pdg     = _lambda_pdg[i];
    _sel_lambda_mass    = _lambda_mass[i];
    _sel_lambda_energy  = _lambda_energy[i];
    _sel_lambda_px      = _lambda_px[i];
    _sel_lambda_py      = _lambda_py[i];
    _sel_lambda_pz      = _lambda_pz[i];
    _sel_lambda_vx      = _lambda_vx[i];
    _sel_lambda_vy      = _lambda_vy[i];
    _sel_lambda_vz      = _lambda_vz[i];
    _sel_lambda_endx    = _lambda_endx[i];
    _sel_lambda_endy    = _lambda_endy[i];
    _sel_lambda_endz    = _lambda_endz[i];
    _sel_lambda_len     = _lambda_len[i];
    _sel_lambda_ct      = _lambda_ct[i];
    _sel_lambda_dist_to_nu_vtx = _lambda_dist_to_nu_vtx[i];
    _sel_lambda_decay_in_fid   = (_lambda_decay_in_fid[i] != 0);

    _sel_lambda_from_sigma0      = (_lambda_from_sigma0[i] != 0);
    _sel_lambda_heavy_feed       = (_lambda_heavy_feed[i]  != 0);
    _sel_lambda_parent_pdg       = _lambda_parent_pdg[i];
    _sel_lambda_grandparent_pdg  = _lambda_grandparent_pdg[i];

    _sel_proton_trackid = _lambda_proton_trackid[i];
    _sel_pion_trackid   = _lambda_pion_trackid[i];

    _sel_proton_energy = _proton_energy[i];
    _sel_proton_px     = _proton_px[i];
    _sel_proton_py     = _proton_py[i];
    _sel_proton_pz     = _proton_pz[i];

    _sel_pion_energy   = _pion_energy[i];
    _sel_pion_px       = _pion_px[i];
    _sel_pion_py       = _pion_py[i];
    _sel_pion_pz       = _pion_pz[i];

    _sel_ppi_invariant_mass = _ppi_invariant_mass[i];
    _sel_ppi_opening_angle  = _ppi_opening_angle[i];
}

void SignalAnalysis::FindTruthMuon(const art::ValidHandle<std::vector<simb::MCParticle>>& mcp_h) {
    _mu_truth_trackid = -1; _mu_truth_pdg = 0;
    double bestE = -1.0;
    for (size_t i = 0; i < mcp_h->size(); ++i) {
        const auto& p = mcp_h->at(i);
        if (std::abs(p.PdgCode()) != 13) continue;
        if (!(p.Process() == "primary" && p.StatusCode() == 1)) continue;
        if (p.E() > bestE) {
            bestE = p.E();
            _mu_truth_trackid = p.TrackId();
            _mu_truth_pdg     = p.PdgCode();
            _mu_truth_E  = p.E();
            _mu_truth_px = p.Px(); _mu_truth_py = p.Py(); _mu_truth_pz = p.Pz();
            _mu_truth_vx = p.Vx(); _mu_truth_vy = p.Vy(); _mu_truth_vz = p.Vz();
        }
    }
    if (_mu_truth_trackid >= 0) {
        _mu_truth_pt = std::sqrt(_mu_truth_px*_mu_truth_px + _mu_truth_py*_mu_truth_py);
        TVector3 nu_dir(_true_nu_dir_x, _true_nu_dir_y, _true_nu_dir_z);
        TVector3 mu_p(_mu_truth_px, _mu_truth_py, _mu_truth_pz);
        _mu_truth_theta = (nu_dir.Mag() > 0 && mu_p.Mag() > 0) ? nu_dir.Angle(mu_p) : nan<float>();
    } else {
        _mu_truth_E = nan<float>(); _mu_truth_px = _mu_truth_py = _mu_truth_pz = nan<float>();
        _mu_truth_pt = nan<float>(); _mu_truth_theta = nan<float>();
        _mu_truth_vx = _mu_truth_vy = _mu_truth_vz = nan<float>();
    }
}

size_t SignalAnalysis::CountTruthHitsInSlice(
    int trackid,
    const std::vector<art::Ptr<recob::Hit>>& hits,
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& assocMCPart) const
{
    if (trackid < 0) return 0;
    size_t n = 0;
    for (auto const& h : hits) {
        auto parts = assocMCPart->at(h.key());
        for (auto const& mcp : parts) {
            if (mcp->TrackId() == trackid) { ++n; break; }
        }
    }
    return n;
}

void SignalAnalysis::analyseEvent(const art::Event& event, bool is_data) {
    this->resetTTree(nullptr);
    if (is_data) return;

    auto const& mct_h = event.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    if (!mct_h.isValid() || mct_h->empty()) return;
    const auto& mct = mct_h->front();

    if (mct.NeutrinoSet()) {
        const auto nuinfo = mct.GetNeutrino();
        const auto& nu    = nuinfo.Nu();
        _nu_pdg = nu.PdgCode();
        _ccnc   = nuinfo.CCNC();
        _interaction_mode = nuinfo.Mode();
        _interaction_type = nuinfo.InteractionType();
        _nu_E    = nu.Trajectory().E(0);
        _nu_theta= nu.Momentum().Vect().Theta();
        _nu_pt   = nuinfo.Pt();
        _nu_t    = nu.T();
        _nu_vx   = nu.Vx(); _nu_vy = nu.Vy(); _nu_vz = nu.Vz();

        const double pmag = std::sqrt(nu.Px()*nu.Px() + nu.Py()*nu.Py() + nu.Pz()*nu.Pz());
        if (pmag > 0) {
            _true_nu_dir_x = nu.Px()/pmag;
            _true_nu_dir_y = nu.Py()/pmag;
            _true_nu_dir_z = nu.Pz()/pmag;
        }

        float sce[3];
        common::True2RecoMappingXYZ(_nu_t, _nu_vx, _nu_vy, _nu_vz, sce);
        _nu_sce_vx = sce[0]; _nu_sce_vy = sce[1]; _nu_sce_vz = sce[2];
        double vtx[3] = {_nu_vx,_nu_vy,_nu_vz};
        _nu_vtx_in_fid = common::isFiducial(vtx, fFidvolXstart, fFidvolYstart, fFidvolZstart,
                                                  fFidvolXend,   fFidvolYend,   fFidvolZend);
    }

    const bool is_numu_like = (std::abs(_nu_pdg) == 14);
    _is_nu_mu_cc = (!fRequireNuMu || is_numu_like) && (!fRequireCC || _ccnc == 0);

    auto const& mcp_h = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mp;
    for (size_t i = 0; i < mcp_h->size(); ++i)
        mp[mcp_h->at(i).TrackId()] = art::Ptr<simb::MCParticle>(mcp_h, i);

    _n_sigma0_truth = 0;
    _n_kshort_truth = 0;
    for (size_t i = 0; i < mcp_h->size(); ++i) {
        const auto& p = mcp_h->at(i);
        const int apdg = std::abs(p.PdgCode());
        if (apdg == 3212) ++_n_sigma0_truth;
        if (p.PdgCode() == 310) ++_n_kshort_truth;
    }
    _has_sigma0_truth = (_n_sigma0_truth > 0);
    _has_kshort_truth = (_n_kshort_truth > 0);

    auto isSigma0 = [](int pdg){ return std::abs(pdg) == 3212; };
    auto isHeavyStrangeBaryon = [](int pdg){ int a=std::abs(pdg); return (a==3312 || a==3322 || a==3334); };
    auto isCharmHadron = [](int pdg){ int a=std::abs(pdg); return ((a>400 && a<500) || (a>4000 && a<5000)); };
    auto isTau = [](int pdg){ return std::abs(pdg)==15; };

    for (size_t i = 0; i < mcp_h->size(); ++i) {
        const auto lam_ptr = art::Ptr<simb::MCParticle>(mcp_h, i);
        const auto& lam    = *lam_ptr;
        if (std::abs(lam.PdgCode()) != 3122) continue;
        if (!fAllowAntiLambda && lam.PdgCode() != 3122) continue;
        if (fRequirePrimaryLambda && lam.Process() != "primary") continue;

        DecayMatch dm = MatchLambdaToPPi(lam_ptr, mp, fAllowAntiLambda);
        if (fRequireSpecificDecay && !dm.ok) continue;

        const simb::MCParticle* p  = (dm.p_trkid  >= 0 && mp.count(dm.p_trkid))  ? mp.at(dm.p_trkid).get()  : nullptr;
        const simb::MCParticle* pi = (dm.pi_trkid >= 0 && mp.count(dm.pi_trkid)) ? mp.at(dm.pi_trkid).get() : nullptr;
        if (p && fProtonKEThreshold > 0.f) {
            const float pKE = p->E() - p->Mass();
            if (pKE < fProtonKEThreshold) continue;
        }
        if (pi && fPionKEThreshold > 0.f) {
            const float piKE = pi->E() - pi->Mass();
            if (piKE < fPionKEThreshold) continue;
        }

        const float L = ThreeDistance(lam.EndX(), lam.EndY(), lam.EndZ(), lam.Vx(), lam.Vy(), lam.Vz());
        const float p_mag = std::sqrt(lam.Px()*lam.Px() + lam.Py()*lam.Py() + lam.Pz()*lam.Pz());
        const float ct = (p_mag > 0.f) ? (L * lam.Mass() / p_mag) : nan<float>();
        const float dist2nu = ThreeDistance(lam.Vx(), lam.Vy(), lam.Vz(), _nu_vx, _nu_vy, _nu_vz);
        double dpos[3] = {lam.EndX(), lam.EndY(), lam.EndZ()};
        const bool decay_in_fid = common::isFiducial(dpos, fFidvolXstart, fFidvolYstart, fFidvolZstart,
                                                            fFidvolXend,   fFidvolYend,   fFidvolZend);

        _lambda_trackid.push_back(lam.TrackId());
        _lambda_pdg.push_back(lam.PdgCode());
        _lambda_mother_trackid.push_back(lam.Mother());
        _lambda_is_primary.push_back(lam.Process() == "primary" ? 1 : 0);
        _lambda_mass.push_back(lam.Mass());
        _lambda_energy.push_back(lam.E());
        _lambda_px.push_back(lam.Px());
        _lambda_py.push_back(lam.Py());
        _lambda_pz.push_back(lam.Pz());
        _lambda_vx.push_back(lam.Vx());
        _lambda_vy.push_back(lam.Vy());
        _lambda_vz.push_back(lam.Vz());
        _lambda_endx.push_back(lam.EndX());
        _lambda_endy.push_back(lam.EndY());
        _lambda_endz.push_back(lam.EndZ());
        _lambda_len.push_back(L);
        _lambda_ct.push_back(ct);
        _lambda_dist_to_nu_vtx.push_back(dist2nu);
        _lambda_decay_in_fid.push_back(decay_in_fid ? 1 : 0);

        unsigned int bits = 0u;
        if (lam.Process() == "primary") bits |= (1u<<0); else bits |= (1u<<1);

        bool from_sigma0 = false;
        bool heavy_feed  = false;
        int  parent_pdg = 0;
        int  grandparent_pdg = 0;

        if (mp.count(lam.Mother())) {
            auto cur = mp.at(lam.Mother());
            parent_pdg = cur->PdgCode();

            if (isSigma0(parent_pdg)) {
                from_sigma0 = true;
                bits |= (1u<<2);
            }

            int steps = 0;
            while (cur && steps < fMaxAncestorSteps) {
                if (steps == 0) {
                } else if (steps == 1) {
                    grandparent_pdg = cur->PdgCode();
                }

                const int a = cur->PdgCode();
                if (isHeavyStrangeBaryon(a) || isCharmHadron(a) || isTau(a)) {
                    heavy_feed = true;
                    bits |= (1u<<3);
                }

                if (!mp.count(cur->Mother())) break;
                cur = mp.at(cur->Mother());
                ++steps;
            }
        }

        _lambda_origin_bits.push_back(bits);
        _lambda_from_sigma0.push_back(from_sigma0 ? 1 : 0);
        _lambda_heavy_feed.push_back(heavy_feed ? 1 : 0);
        _lambda_parent_pdg.push_back(parent_pdg);
        _lambda_grandparent_pdg.push_back(grandparent_pdg);

        if (heavy_feed) ++_n_lambda_from_heavy;

        _lambda_proton_trackid.push_back(dm.p_trkid);
        _lambda_pion_trackid.push_back(dm.pi_trkid);

        auto fill_d = [&](const simb::MCParticle* d,
                          std::vector<float>& E, std::vector<float>& px,
                          std::vector<float>& py, std::vector<float>& pz,
                          std::vector<float>& vx, std::vector<float>& vy,
                          std::vector<float>& vz) {
            if (!d) { E.push_back(nan<float>()); px.push_back(nan<float>());
                      py.push_back(nan<float>()); pz.push_back(nan<float>());
                      vx.push_back(nan<float>()); vy.push_back(nan<float>()); vz.push_back(nan<float>());
                      return;
            }
            E.push_back(d->E());
            px.push_back(d->Px()); py.push_back(d->Py()); pz.push_back(d->Pz());
            vx.push_back(d->Vx()); vy.push_back(d->Vy()); vz.push_back(d->Vz());
        };
        fill_d(p,  _proton_energy, _proton_px, _proton_py, _proton_pz, _proton_vx, _proton_vy, _proton_vz);
        fill_d(pi, _pion_energy,   _pion_px,   _pion_py,   _pion_pz,   _pion_vx,   _pion_vy,   _pion_vz);

        float invM = nan<float>(), open = nan<float>();
        if (p && pi) {
            TLorentzVector Pp, Ppi;
            Pp.SetXYZM(p->Px(),  p->Py(),  p->Pz(),  p->Mass());
            Ppi.SetXYZM(pi->Px(), pi->Py(), pi->Pz(), pi->Mass());
            invM = (Pp + Ppi).M();
            TVector3 vp(p->Px(), p->Py(), p->Pz());
            TVector3 vpi(pi->Px(), pi->Py(), pi->Pz());
            if (vp.Mag() > 0 && vpi.Mag() > 0) open = vp.Angle(vpi);
        }
        _ppi_invariant_mass.push_back(invM);
        _ppi_opening_angle.push_back(open);
    }

    _n_lambda_candidates = static_cast<int>(_lambda_trackid.size());
    _has_lambda_to_ppi = (_n_lambda_candidates > 0);

    if (_has_lambda_to_ppi) SelectBestCandidate();

    const bool accept_nu_flavor = fAcceptAntiNu ? (std::abs(_nu_pdg)==14) : (_nu_pdg==14);
    _is_nu_mu_cc = (!fRequireNuMu || accept_nu_flavor) && (!fRequireCC || _ccnc == 0);
    _is_signal_event = _is_nu_mu_cc && _has_lambda_to_ppi;

    FindTruthMuon(mcp_h);

    bool any_lambda_decay_infid = false;
    for (auto v : _lambda_decay_in_fid) if (v) { any_lambda_decay_infid = true; break; }
    _pr_eligible_event = _is_nu_mu_cc && _nu_vtx_in_fid && _has_lambda_to_ppi && any_lambda_decay_infid;
}

void SignalAnalysis::analyseSlice(const art::Event& event,
                                        std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                                        bool is_data, bool) {
    if (is_data) return;

    _pr_n_pfps_in_slice = static_cast<int>(slice_pfp_vec.size());
    _pr_valid_assignment = false;
    _pr_pass_event = false;

    auto const& cluster_h = event.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
    art::FindManyP<recob::Hit> assocHits(cluster_h, event, fCLSproducer);

    std::vector<std::vector<art::Ptr<recob::Hit>>> pfp_hits(slice_pfp_vec.size());
    std::vector<art::Ptr<recob::Hit>> inputHits;
    inputHits.reserve(1024);

    for (size_t ip = 0; ip < slice_pfp_vec.size(); ++ip) {
        auto clusters = slice_pfp_vec[ip].get<recob::Cluster>();
        for (auto const& c : clusters) {
            auto hits = assocHits.at(c.key());
            pfp_hits[ip].insert(pfp_hits[ip].end(), hits.begin(), hits.end());
            inputHits.insert(inputHits.end(), hits.begin(), hits.end());
        }
    }

    auto const& hit_h = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    auto assocMCPart = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_h, event, fBKTproducer);

    {
        auto const& mcshower_h = event.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer);
        auto const& mctrack_h  = event.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer);

        std::vector<recob::Hit> inputHits_deref;
        inputHits_deref.reserve(inputHits.size());
        for (auto const& hptr : inputHits) inputHits_deref.push_back(*hptr);

        std::vector<common::BtPart> btparts_v =
            common::initBacktrackingParticleVec(*mcshower_h, *mctrack_h, inputHits_deref, assocMCPart);

        auto computeMetrics = [&](int trackid, float& purity, float& completeness) {
            purity = 0.f; completeness = 0.f;
            if (trackid < 0) return false;

            auto it = std::find_if(btparts_v.begin(), btparts_v.end(),
                                   [trackid](const common::BtPart& bt) {
                                       return std::find(bt.tids.begin(), bt.tids.end(), trackid) != bt.tids.end();
                                   });
            if (it == btparts_v.end()) return false;

            const auto& bt = *it;
            size_t n_shared = 0;
            for (auto const& hptr : inputHits) {
                auto parts = assocMCPart->at(hptr.key());
                for (auto const& mcp : parts) {
                    if (std::find(bt.tids.begin(), bt.tids.end(), mcp->TrackId()) != bt.tids.end()) {
                        ++n_shared; break;
                    }
                }
            }
            completeness = bt.nhits > 0 ? static_cast<float>(n_shared) / bt.nhits : 0.f;
            purity       = !inputHits.empty() ? static_cast<float>(n_shared) / inputHits.size() : 0.f;
            return true;
        };

        (void)computeMetrics(_sel_proton_trackid, _sel_proton_purity, _sel_proton_completeness);
        (void)computeMetrics(_sel_pion_trackid,   _sel_pion_purity,   _sel_pion_completeness);
    }

    if (_mu_truth_trackid < 0 || _sel_proton_trackid < 0 || _sel_pion_trackid < 0) return;

    const int Tids[3] = {_mu_truth_trackid, _sel_proton_trackid, _sel_pion_trackid};
    int N_hits_T[3] = {0,0,0};
    for (int ti=0; ti<3; ++ti)
        N_hits_T[ti] = static_cast<int>(CountTruthHitsInSlice(Tids[ti], inputHits, assocMCPart));

    const int N = static_cast<int>(pfp_hits.size());
    if (N < 3) return;

    std::vector<std::array<int,3>> S(N, {0,0,0});
    std::vector<int> N_hits_O(N, 0);

    for (int i = 0; i < N; ++i) {
        N_hits_O[i] = static_cast<int>(pfp_hits[i].size());
        for (auto const& h : pfp_hits[i]) {
            auto parts = assocMCPart->at(h.key());
            bool has[3] = {false,false,false};
            for (auto const& mcp : parts) {
                const int tid = mcp->TrackId();
                for (int t=0; t<3; ++t) if (!has[t] && tid == Tids[t]) has[t] = true;
            }
            for (int t=0; t<3; ++t) if (has[t]) ++S[i][t];
        }
    }

    int best_i=-1, best_j=-1, best_k=-1, best_sum=-1;
    for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (j!=i) for (int k=0;k<N;++k) if (k!=i && k!=j) {
        const int sum = S[i][0] + S[j][1] + S[k][2];
        if (sum > best_sum) { best_sum = sum; best_i=i; best_j=j; best_k=k; }
    }
    if (best_sum < 0) return;

    _pr_valid_assignment = true;
    _pr_best_sum_nshared = best_sum;
    _pr_mu_pfp_index = best_i;
    _pr_p_pfp_index  = best_j;
    _pr_pi_pfp_index = best_k;

    auto fill_pair = [&](int i, int t, int& nshared, int& nhitsO, int& nhitsT, float& purity, float& comp) {
        nshared = S[i][t];
        nhitsO  = N_hits_O[i];
        nhitsT  = N_hits_T[t];
        purity  = (nhitsO>0) ? static_cast<float>(nshared)/nhitsO : 0.f;
        comp    = (nhitsT>0) ? static_cast<float>(nshared)/nhitsT : 0.f;
    };

    fill_pair(best_i, 0, _pr_mu_nshared, _pr_mu_nhits_O, _pr_mu_nhits_T, _pr_mu_purity, _pr_mu_completeness);
    fill_pair(best_j, 1, _pr_p_nshared,  _pr_p_nhits_O,  _pr_p_nhits_T,  _pr_p_purity,  _pr_p_completeness);
    fill_pair(best_k, 2, _pr_pi_nshared, _pr_pi_nhits_O, _pr_pi_nhits_T, _pr_pi_purity, _pr_pi_completeness);

    const bool pass_mu = (_pr_mu_completeness > fPRMinCompleteness) && (_pr_mu_purity > fPRMinPurity);
    const bool pass_p  = (_pr_p_completeness  > fPRMinCompleteness) && (_pr_p_purity  > fPRMinPurity);
    const bool pass_pi = (_pr_pi_completeness > fPRMinCompleteness) && (_pr_pi_purity > fPRMinPurity);

    _pr_pass_event = (pass_mu && pass_p && pass_pi);
}

DEFINE_ART_CLASS_TOOL(SignalAnalysis)

}

#endif

