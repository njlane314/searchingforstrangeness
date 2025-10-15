#ifndef ANALYSIS_TRUTH_CXX
#define ANALYSIS_TRUTH_CXX

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
#include "larpandora/LArPandoraEventBuilding/LArPandoraSliceIdHelper.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"
#include "Common/GeometryUtils.h"
#include "Common/ParticleScattering.h"
#include "Common/SpaceChargeCorrections.h"
#include "Common/TruthContainment.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"
#include <cmath>
#include <limits>
#include <string>
#include <vector>
// Explicit STL deps used below
#include <map>
#include <algorithm>

namespace analysis {

class TruthAnalysis : public AnalysisToolBase {
public:
    explicit TruthAnalysis(fhicl::ParameterSet const& p);
    virtual ~TruthAnalysis() = default;
    void configure(const fhicl::ParameterSet& pset) override;
    void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec, bool is_data, bool is_selected) override;
    void analyseEvent(const art::Event& event, bool is_data) override;
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    art::InputTag fMCTproducer;
    art::InputTag fMCFproducer;
    art::InputTag fMCPproducer;
    art::InputTag fHITproducer;
    art::InputTag fBKTproducer;
    art::InputTag fMCRproducer;
    art::InputTag fCLSproducer;

    float fFidvolXstart;
    float fFidvolXend;
    float fFidvolYstart;
    float fFidvolYend;
    float fFidvolZstart;
    float fFidvolZend;

    int _neutrino_pdg;
    int _interaction_ccnc;
    int _interaction_mode;
    int _interaction_type;
    float _neutrino_energy;
    float _neutrino_theta;
    float _neutrino_pt;
    int _target_nucleus_pdg;
    int _hit_nucleon_pdg;
    float _kinematic_W;
    float _kinematic_X;
    float _kinematic_Y;
    float _kinematic_Q_squared;
    float _neutrino_momentum_x;
    float _neutrino_momentum_y;
    float _neutrino_momentum_z;
    float _neutrino_vertex_x;
    float _neutrino_vertex_y;
    float _neutrino_vertex_z;
    float _neutrino_vertex_wire_u;
    float _neutrino_vertex_wire_v;
    float _neutrino_vertex_wire_w;
    float _neutrino_vertex_time;
    float _neutrino_sce_vertex_x;
    float _neutrino_sce_vertex_y;
    float _neutrino_sce_vertex_z;
    float _lepton_energy;
    float _true_neutrino_momentum_x;
    float _true_neutrino_momentum_y;
    float _true_neutrino_momentum_z;
    float _flux_path_length;
    int _flux_parent_pdg;
    int _flux_hadron_pdg;
    int _flux_decay_mode;
    double _flux_decay_vtx_x;
    double _flux_decay_vtx_y;
    double _flux_decay_vtx_z;
    double _flux_decay_mom_x;
    double _flux_decay_mom_y;
    double _flux_decay_mom_z;
    double _numi_baseline;
    float _numi_off_axis_angle;
    double _bnb_baseline;
    float _bnb_off_axis_angle;
    bool _is_vertex_in_fiducial;
    int _count_mu_minus;
    int _count_mu_plus;
    int _count_e_minus;
    int _count_e_plus;
    int _count_pi_zero;
    int _count_pi_plus;
    int _count_pi_minus;
    int _count_kaon_plus;
    int _count_kaon_minus;
    int _count_kaon_zero;
    int _count_proton;
    int _count_neutron;
    int _count_gamma;
    int _count_lambda;
    int _count_sigma_plus;
    int _count_sigma_zero;
    int _count_sigma_minus;
    std::vector<int> _mc_particle_pdg;
    std::vector<int> _mc_particle_trackid;
    std::vector<float> _mc_particle_energy;
    std::vector<uint> _mc_elastic_scatters;
    std::vector<uint> _mc_inelastic_scatters;
    std::vector<float> _mc_momentum_x;
    std::vector<float> _mc_momentum_y;
    std::vector<float> _mc_momentum_z;
    std::vector<float> _mc_end_momentum;
    std::vector<float> _mc_start_vertex_x;
    std::vector<float> _mc_start_vertex_y;
    std::vector<float> _mc_start_vertex_z;
    std::vector<float> _mc_end_vertex_x;
    std::vector<float> _mc_end_vertex_y;
    std::vector<float> _mc_end_vertex_z;
    std::vector<std::string> _mc_particle_final_state;
    std::vector<float> _mc_completeness;
    std::vector<float> _mc_purity;
    std::vector<std::vector<int>> _mc_daughter_pdg;
    std::vector<std::vector<float>> _mc_daughter_energy;
    std::vector<std::string> _mc_daughter_process_flat;
    std::vector<size_t>    _mc_daughter_process_idx;
    std::vector<std::vector<float>> _mc_daughter_mom_x;
    std::vector<std::vector<float>> _mc_daughter_mom_y;
    std::vector<std::vector<float>> _mc_daughter_mom_z;
    std::vector<std::vector<float>> _mc_daughter_vtx_x;
    std::vector<std::vector<float>> _mc_daughter_vtx_y;
    std::vector<std::vector<float>> _mc_daughter_vtx_z;
    float _true_transverse_momentum;
    float _true_visible_transverse_momentum;
    float _true_total_momentum;
    float _true_visible_total_momentum;
    float _true_visible_energy;
};

TruthAnalysis::TruthAnalysis(fhicl::ParameterSet const& p) {
    this->configure(p);
}

void TruthAnalysis::configure(fhicl::ParameterSet const& p) {
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCFproducer = p.get<art::InputTag>("MCFproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fBKTproducer = p.get<art::InputTag>("BKTproducer");
    fMCRproducer = p.get<art::InputTag>("MCRproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");

    fFidvolXstart = p.get<double>("fidvolXstart", 10);
    fFidvolXend = p.get<double>("fidvolXend", 10);
    fFidvolYstart = p.get<double>("fidvolYstart", 15);
    fFidvolYend = p.get<double>("fidvolYend", 15);
    fFidvolZstart = p.get<double>("fidvolZstart", 10);
    fFidvolZend = p.get<double>("fidvolZend", 50);

}

void TruthAnalysis::setBranches(TTree* _tree) {
    // ---- Interaction & neutrino summary ----
    _tree->Branch("nu_pdg",              &_neutrino_pdg,          "nu_pdg/I");
    _tree->Branch("int_ccnc",            &_interaction_ccnc,      "int_ccnc/I");
    _tree->Branch("int_mode",            &_interaction_mode,      "int_mode/I");
    _tree->Branch("int_type",            &_interaction_type,      "int_type/I");
    _tree->Branch("nu_E",                &_neutrino_energy,       "nu_E/F");
    _tree->Branch("nu_theta",            &_neutrino_theta,        "nu_theta/F");
    _tree->Branch("nu_pt",               &_neutrino_pt,           "nu_pt/F");
    _tree->Branch("tgt_pdg",             &_target_nucleus_pdg,    "tgt_pdg/I");
    _tree->Branch("struck_nucleon_pdg",  &_hit_nucleon_pdg,       "struck_nucleon_pdg/I");

    _tree->Branch("kin_W",               &_kinematic_W,           "kin_W/F");
    _tree->Branch("bjorken_x",           &_kinematic_X,           "bjorken_x/F");
    _tree->Branch("inelasticity_y",      &_kinematic_Y,           "inelasticity_y/F");
    _tree->Branch("Q2",                  &_kinematic_Q_squared,   "Q2/F");

    // ---- Neutrino 4-mom + vertex ----
    _tree->Branch("nu_px",               &_neutrino_momentum_x,   "nu_px/F");
    _tree->Branch("nu_py",               &_neutrino_momentum_y,   "nu_py/F");
    _tree->Branch("nu_pz",               &_neutrino_momentum_z,   "nu_pz/F");

    _tree->Branch("nu_vtx_x",            &_neutrino_vertex_x,     "nu_vtx_x/F");
    _tree->Branch("nu_vtx_y",            &_neutrino_vertex_y,     "nu_vtx_y/F");
    _tree->Branch("nu_vtx_z",            &_neutrino_vertex_z,     "nu_vtx_z/F");
    _tree->Branch("nu_vtx_wire_u",       &_neutrino_vertex_wire_u,"nu_vtx_wire_u/F");
    _tree->Branch("nu_vtx_wire_v",       &_neutrino_vertex_wire_v,"nu_vtx_wire_v/F");
    _tree->Branch("nu_vtx_wire_w",       &_neutrino_vertex_wire_w,"nu_vtx_wire_w/F");
    _tree->Branch("nu_vtx_t",            &_neutrino_vertex_time,  "nu_vtx_t/F");

    _tree->Branch("nu_vtx_sce_x",        &_neutrino_sce_vertex_x, "nu_vtx_sce_x/F");
    _tree->Branch("nu_vtx_sce_y",        &_neutrino_sce_vertex_y, "nu_vtx_sce_y/F");
    _tree->Branch("nu_vtx_sce_z",        &_neutrino_sce_vertex_z, "nu_vtx_sce_z/F");

    _tree->Branch("lep_E",               &_lepton_energy,         "lep_E/F");

    // (kept) duplicate of nu_p{xyz} for backward compat, renamed clearly
    _tree->Branch("nu_true_px",          &_true_neutrino_momentum_x, "nu_true_px/F");
    _tree->Branch("nu_true_py",          &_true_neutrino_momentum_y, "nu_true_py/F");
    _tree->Branch("nu_true_pz",          &_true_neutrino_momentum_z, "nu_true_pz/F");

    // ---- Flux / beam geometry ----
    _tree->Branch("flux_path_len",       &_flux_path_length,      "flux_path_len/F");
    _tree->Branch("flux_parent_pdg",     &_flux_parent_pdg,       "flux_parent_pdg/I");
    _tree->Branch("flux_hadron_pdg",     &_flux_hadron_pdg,       "flux_hadron_pdg/I");
    _tree->Branch("flux_decay_mode",     &_flux_decay_mode,       "flux_decay_mode/I");
    _tree->Branch("flux_decay_vtx_x",    &_flux_decay_vtx_x,      "flux_decay_vtx_x/D");
    _tree->Branch("flux_decay_vtx_y",    &_flux_decay_vtx_y,      "flux_decay_vtx_y/D");
    _tree->Branch("flux_decay_vtx_z",    &_flux_decay_vtx_z,      "flux_decay_vtx_z/D");
    _tree->Branch("flux_decay_px",       &_flux_decay_mom_x,      "flux_decay_px/D");
    _tree->Branch("flux_decay_py",       &_flux_decay_mom_y,      "flux_decay_py/D");
    _tree->Branch("flux_decay_pz",       &_flux_decay_mom_z,      "flux_decay_pz/D");

    _tree->Branch("numi_L",              &_numi_baseline,         "numi_L/D");
    _tree->Branch("numi_offaxis_deg",    &_numi_off_axis_angle,   "numi_offaxis_deg/F");
    _tree->Branch("bnb_L",               &_bnb_baseline,          "bnb_L/D");
    _tree->Branch("bnb_offaxis_deg",     &_bnb_off_axis_angle,    "bnb_offaxis_deg/F");

    _tree->Branch("nu_vtx_in_fv",        &_is_vertex_in_fiducial, "nu_vtx_in_fv/O");

    // ---- Final-state counts (threshold-free) ----
    _tree->Branch("n_mu_minus",          &_count_mu_minus,   "n_mu_minus/I");
    _tree->Branch("n_mu_plus",           &_count_mu_plus,    "n_mu_plus/I");
    _tree->Branch("n_e_minus",           &_count_e_minus,    "n_e_minus/I");
    _tree->Branch("n_e_plus",            &_count_e_plus,     "n_e_plus/I");
    _tree->Branch("n_pi0",               &_count_pi_zero,    "n_pi0/I");
    _tree->Branch("n_pi_plus",           &_count_pi_plus,    "n_pi_plus/I");
    _tree->Branch("n_pi_minus",          &_count_pi_minus,   "n_pi_minus/I");
    _tree->Branch("n_K_plus",            &_count_kaon_plus,  "n_K_plus/I");
    _tree->Branch("n_K_minus",           &_count_kaon_minus, "n_K_minus/I");
    _tree->Branch("n_K0",                &_count_kaon_zero,  "n_K0/I");
    _tree->Branch("n_p",                 &_count_proton,     "n_p/I");
    _tree->Branch("n_n",                 &_count_neutron,    "n_n/I");
    _tree->Branch("n_gamma",             &_count_gamma,      "n_gamma/I");
    _tree->Branch("n_lambda",            &_count_lambda,     "n_lambda/I");
    _tree->Branch("n_sigma_plus",        &_count_sigma_plus, "n_sigma_plus/I");
    _tree->Branch("n_sigma0",            &_count_sigma_zero, "n_sigma0/I");
    _tree->Branch("n_sigma_minus",       &_count_sigma_minus,"n_sigma_minus/I");

    // ---- Primary particles ----
    _tree->Branch("prim_pdg",            "std::vector<int>",               &_mc_particle_pdg);
    _tree->Branch("prim_track_id",       "std::vector<int>",               &_mc_particle_trackid);
    _tree->Branch("prim_E",              "std::vector<float>",             &_mc_particle_energy);
    _tree->Branch("prim_n_scat_elastic", "std::vector<uint>",              &_mc_elastic_scatters);
    _tree->Branch("prim_n_scat_inelastic","std::vector<uint>",             &_mc_inelastic_scatters);
    _tree->Branch("prim_px",             "std::vector<float>",             &_mc_momentum_x);
    _tree->Branch("prim_py",             "std::vector<float>",             &_mc_momentum_y);
    _tree->Branch("prim_pz",             "std::vector<float>",             &_mc_momentum_z);
    _tree->Branch("prim_p_end",          "std::vector<float>",             &_mc_end_momentum);
    _tree->Branch("prim_vx",             "std::vector<float>",             &_mc_start_vertex_x);
    _tree->Branch("prim_vy",             "std::vector<float>",             &_mc_start_vertex_y);
    _tree->Branch("prim_vz",             "std::vector<float>",             &_mc_start_vertex_z);
    _tree->Branch("prim_end_x",          "std::vector<float>",             &_mc_end_vertex_x);
    _tree->Branch("prim_end_y",          "std::vector<float>",             &_mc_end_vertex_y);
    _tree->Branch("prim_end_z",          "std::vector<float>",             &_mc_end_vertex_z);
    _tree->Branch("prim_final_state",    "std::vector<std::string>",       &_mc_particle_final_state);

    // Backtracking metrics
    _tree->Branch("prim_bt_completeness","std::vector<float>",             &_mc_completeness);
    _tree->Branch("prim_bt_purity",      "std::vector<float>",             &_mc_purity);

    // ---- Direct daughters ----
    _tree->Branch("prim_dau_pdg",        "std::vector<std::vector<int>>",  &_mc_daughter_pdg);
    _tree->Branch("prim_dau_E",          "std::vector<std::vector<float>>",&_mc_daughter_energy);
    _tree->Branch("prim_dau_process_flat","std::vector<std::string>",      &_mc_daughter_process_flat);
    _tree->Branch("prim_dau_process_idx","std::vector<size_t>",            &_mc_daughter_process_idx);
    _tree->Branch("prim_dau_px",         "std::vector<std::vector<float>>",&_mc_daughter_mom_x);
    _tree->Branch("prim_dau_py",         "std::vector<std::vector<float>>",&_mc_daughter_mom_y);
    _tree->Branch("prim_dau_pz",         "std::vector<std::vector<float>>",&_mc_daughter_mom_z);
    _tree->Branch("prim_dau_vx",         "std::vector<std::vector<float>>",&_mc_daughter_vtx_x);
    _tree->Branch("prim_dau_vy",         "std::vector<std::vector<float>>",&_mc_daughter_vtx_y);
    _tree->Branch("prim_dau_vz",         "std::vector<std::vector<float>>",&_mc_daughter_vtx_z);

    // ---- Event-level vector sums ----
    _tree->Branch("sum_pt_true",         &_true_transverse_momentum,        "sum_pt_true/F");
    _tree->Branch("sum_pt_vis",          &_true_visible_transverse_momentum,"sum_pt_vis/F");
    _tree->Branch("sum_p_true",          &_true_total_momentum,             "sum_p_true/F");
    _tree->Branch("sum_p_vis",           &_true_visible_total_momentum,     "sum_p_vis/F");
    _tree->Branch("sum_E_vis",           &_true_visible_energy,             "sum_E_vis/F");
}

void TruthAnalysis::resetTTree(TTree* tree) {
    _neutrino_pdg = 0;
    _interaction_ccnc = -1;
    _interaction_mode = -1;
    _interaction_type = -1;
    _neutrino_energy = std::numeric_limits<float>::quiet_NaN();
    _neutrino_theta = std::numeric_limits<float>::quiet_NaN();
    _neutrino_pt = std::numeric_limits<float>::quiet_NaN();
    _target_nucleus_pdg = 0;
    _hit_nucleon_pdg = 0;
    _kinematic_W = std::numeric_limits<float>::quiet_NaN();
    _kinematic_X = std::numeric_limits<float>::quiet_NaN();
    _kinematic_Y = std::numeric_limits<float>::quiet_NaN();
    _kinematic_Q_squared = std::numeric_limits<float>::quiet_NaN();
    _neutrino_momentum_x = std::numeric_limits<float>::quiet_NaN();
    _neutrino_momentum_y = std::numeric_limits<float>::quiet_NaN();
    _neutrino_momentum_z = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_x = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_y = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_z = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_wire_u = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_wire_v = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_wire_w = std::numeric_limits<float>::quiet_NaN();
    _neutrino_vertex_time = std::numeric_limits<float>::quiet_NaN();
    _neutrino_sce_vertex_x = std::numeric_limits<float>::quiet_NaN();
    _neutrino_sce_vertex_y = std::numeric_limits<float>::quiet_NaN();
    _neutrino_sce_vertex_z = std::numeric_limits<float>::quiet_NaN();
    _lepton_energy = std::numeric_limits<float>::quiet_NaN();
    _true_neutrino_momentum_x = std::numeric_limits<float>::quiet_NaN();
    _true_neutrino_momentum_y = std::numeric_limits<float>::quiet_NaN();
    _true_neutrino_momentum_z = std::numeric_limits<float>::quiet_NaN();
    _flux_path_length = std::numeric_limits<float>::quiet_NaN();
    _flux_parent_pdg = 0;
    _flux_hadron_pdg = 0;
    _flux_decay_mode = -1;
    _flux_decay_vtx_x = std::numeric_limits<double>::quiet_NaN();
    _flux_decay_vtx_y = std::numeric_limits<double>::quiet_NaN();
    _flux_decay_vtx_z = std::numeric_limits<double>::quiet_NaN();
    _flux_decay_mom_x = std::numeric_limits<double>::quiet_NaN();
    _flux_decay_mom_y = std::numeric_limits<double>::quiet_NaN();
    _flux_decay_mom_z = std::numeric_limits<double>::quiet_NaN();
    _numi_baseline = std::numeric_limits<double>::quiet_NaN();
    _numi_off_axis_angle = std::numeric_limits<float>::quiet_NaN();
    _bnb_baseline = std::numeric_limits<double>::quiet_NaN();
    _bnb_off_axis_angle = std::numeric_limits<float>::quiet_NaN();
    _is_vertex_in_fiducial = false;
    _count_mu_minus = 0;
    _count_mu_plus = 0;
    _count_e_minus = 0;
    _count_e_plus = 0;
    _count_pi_zero = 0;
    _count_pi_plus = 0;
    _count_pi_minus = 0;
    _count_kaon_plus = 0;
    _count_kaon_minus = 0;
    _count_kaon_zero = 0;
    _count_proton = 0;
    _count_neutron = 0;
    _count_gamma = 0;
    _count_lambda = 0;
    _count_sigma_plus = 0;
    _count_sigma_zero = 0;
    _count_sigma_minus = 0;
    _mc_particle_pdg.clear();
    _mc_particle_trackid.clear();
    _mc_particle_energy.clear();
    _mc_elastic_scatters.clear();
    _mc_inelastic_scatters.clear();
    _mc_momentum_x.clear();
    _mc_momentum_y.clear();
    _mc_momentum_z.clear();
    _mc_end_momentum.clear();
    _mc_start_vertex_x.clear();
    _mc_start_vertex_y.clear();
    _mc_start_vertex_z.clear();
    _mc_end_vertex_x.clear();
    _mc_end_vertex_y.clear();
    _mc_end_vertex_z.clear();
    _mc_particle_final_state.clear();
    _mc_completeness.clear();
    _mc_purity.clear();
    _mc_daughter_pdg.clear();
    _mc_daughter_energy.clear();
    _mc_daughter_process_flat.clear();
    _mc_daughter_process_idx.clear();
    _mc_daughter_mom_x.clear();
    _mc_daughter_mom_y.clear();
    _mc_daughter_mom_z.clear();
    _mc_daughter_vtx_x.clear();
    _mc_daughter_vtx_y.clear();
    _mc_daughter_vtx_z.clear();
    _true_transverse_momentum = 0;
    _true_visible_transverse_momentum = 0;
    _true_total_momentum = 0;
    _true_visible_total_momentum = 0;
    _true_visible_energy = 0;
}

void TruthAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec, bool is_data, bool is_selected) {
    if (is_data) return;

    auto const& cluster_h = event.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
    art::FindManyP<recob::Hit> assocHits(cluster_h, event, fCLSproducer);

    std::vector<art::Ptr<recob::Hit>> inputHits;
    for (const auto& pfp : slice_pfp_vec) {
        auto clusters = pfp.get<recob::Cluster>();
        for (const auto& cluster : clusters) {
            auto hits = assocHits.at(cluster.key());
            inputHits.insert(inputHits.end(), hits.begin(), hits.end());
        }
    }

    auto const& hit_h = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    auto assocMCPart = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_h, event, fBKTproducer);

    auto const& mcshower_h = event.getValidHandle<std::vector<sim::MCShower>>(fMCRproducer);
    auto const& mctrack_h = event.getValidHandle<std::vector<sim::MCTrack>>(fMCRproducer);

    std::vector<recob::Hit> inputHits_deref;
    inputHits_deref.reserve(inputHits.size());
    for (const auto& hit_ptr : inputHits) {
        inputHits_deref.push_back(*hit_ptr);
    }

    std::vector<common::BtPart> btparts_v = common::initBacktrackingParticleVec(*mcshower_h, *mctrack_h, inputHits_deref, assocMCPart);

    auto computeMetrics = [&](int trackid, float& purity, float& completeness) -> bool {
        auto it = std::find_if(btparts_v.begin(), btparts_v.end(), [trackid](const common::BtPart& bt) {
             return std::find(bt.tids.begin(), bt.tids.end(), trackid) != bt.tids.end();
        });

        if (it == btparts_v.end()) return false;

        const auto& bt = *it;
        size_t n_shared = 0;

        for (const auto& hit_ptr : inputHits) {
            std::vector<art::Ptr<simb::MCParticle>> particle_vec = assocMCPart->at(hit_ptr.key());
            for (const auto& mcp : particle_vec) {
                if (std::find(bt.tids.begin(), bt.tids.end(), mcp->TrackId()) != bt.tids.end()) {
                    n_shared++;
                    break; 
                }
            }
        }

        completeness = bt.nhits > 0 ? static_cast<float>(n_shared) / bt.nhits : 0.f;
        purity = inputHits.size() > 0 ? static_cast<float>(n_shared) / inputHits.size() : 0.f;
        return true;
    };

    for (size_t i = 0; i < _mc_particle_pdg.size(); ++i) {
        float purity = 0.0;
        float completeness = 0.0;
        if (computeMetrics(_mc_particle_trackid[i], purity, completeness)) {
            _mc_purity[i] = purity;
            _mc_completeness[i] = completeness;
        } else {
            _mc_purity[i] = 0.0;
            _mc_completeness[i] = 0.0;
        }
    }

}

void TruthAnalysis::analyseEvent(const art::Event& event, bool is_data) {
    if (is_data) {
        this->resetTTree(nullptr);
        return;
    }

    auto const& mct_h = event.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    auto mct = mct_h->at(0);

    auto const& mcflux_h = event.getValidHandle<std::vector<simb::MCFlux>>(fMCFproducer);
    if (mcflux_h.isValid() && !mcflux_h->empty()) {
        auto flux = mcflux_h->at(0);
        
        _flux_parent_pdg = flux.fptype;
        _flux_hadron_pdg = flux.ftptype;
        _flux_decay_mode = flux.fndecay;
        _flux_decay_vtx_x = flux.fvx;
        _flux_decay_vtx_y = flux.fvy;
        _flux_decay_vtx_z = flux.fvz;
        _flux_decay_mom_x = flux.fpdpx;
        _flux_decay_mom_y = flux.fpdpy;
        _flux_decay_mom_z = flux.fpdpz;
        _flux_path_length = flux.fdk2gen + flux.fgen2vtx;
        
        TVector3 decay_position(flux.fvx, flux.fvy, flux.fvz);

        TVector3 numi_target_to_detector_beam(5502, 7259, 67270);
        TVector3 numi_baseline_vector = numi_target_to_detector_beam - decay_position;
        _numi_baseline = numi_baseline_vector.Mag() / 100.0;
        if (numi_baseline_vector.Mag() > 0) {
            _numi_off_axis_angle = std::acos(numi_baseline_vector.Unit().Dot(numi_target_to_detector_beam.Unit())) * 180.0 / M_PI;
        }

        TVector3 bnb_target_to_detector_beam(0.0, 0.0, 47000.0);
        TVector3 bnb_baseline_vector = bnb_target_to_detector_beam - decay_position;
        _bnb_baseline = bnb_baseline_vector.Mag() / 100.0;
        if (bnb_baseline_vector.Mag() > 0) {
            _bnb_off_axis_angle = std::acos(bnb_baseline_vector.Unit().Dot(bnb_target_to_detector_beam.Unit())) * 180.0 / M_PI;
        }
    }

    if (mct.NeutrinoSet()) {
        auto neutrino = mct.GetNeutrino();
        auto nu = neutrino.Nu();
        _interaction_ccnc = neutrino.CCNC();
        _interaction_type = neutrino.InteractionType();
        _interaction_mode = neutrino.Mode();
        _neutrino_pdg = nu.PdgCode();
        _neutrino_energy = nu.Trajectory().E(0);
        _lepton_energy = neutrino.Lepton().E();
        _neutrino_vertex_time = nu.T();
        _neutrino_vertex_x = nu.Vx();
        _neutrino_vertex_y = nu.Vy();
        _neutrino_vertex_z = nu.Vz();
        _true_neutrino_momentum_x = nu.Px();
        _true_neutrino_momentum_y = nu.Py();
        _true_neutrino_momentum_z = nu.Pz();
        float neutrino_vertex_sce[3];
        common::True2RecoMappingXYZ(_neutrino_vertex_time, _neutrino_vertex_x, _neutrino_vertex_y, _neutrino_vertex_z, neutrino_vertex_sce);
        _neutrino_sce_vertex_x = neutrino_vertex_sce[0];
        _neutrino_sce_vertex_y = neutrino_vertex_sce[1];
        _neutrino_sce_vertex_z = neutrino_vertex_sce[2];
        _neutrino_theta = neutrino.Theta();
        _neutrino_pt = neutrino.Pt();
        _target_nucleus_pdg = neutrino.Target();
        _hit_nucleon_pdg = neutrino.HitNuc();
        _kinematic_W = neutrino.W();
        _kinematic_X = neutrino.X();
        _kinematic_Y = neutrino.Y();
        _kinematic_Q_squared = neutrino.QSqr();
        _neutrino_momentum_x = nu.Px();
        _neutrino_momentum_y = nu.Py();
        _neutrino_momentum_z = nu.Pz();
        double vertex[3] = {_neutrino_vertex_x, _neutrino_vertex_y, _neutrino_vertex_z};
        _is_vertex_in_fiducial = common::isFiducial(vertex, fFidvolXstart, fFidvolYstart, fFidvolZstart,
                                                    fFidvolXend, fFidvolYend, fFidvolZend);
    }

    TLorentzVector total_momentum;
    TLorentzVector total_visible_momentum;
    size_t n_particles = mct.NParticles();
    for (size_t i = 0; i < n_particles; i++) {
        auto const& particle = mct.GetParticle(i);
        if (!(particle.StatusCode() == 1 && particle.Process() == "primary")) continue;
        total_momentum += particle.Momentum(0);
        total_visible_momentum += particle.Momentum(0);
    }
    _true_transverse_momentum = total_momentum.Perp();
    _true_visible_transverse_momentum = total_visible_momentum.Perp();
    _true_total_momentum = total_momentum.Mag();
    _true_visible_total_momentum = total_visible_momentum.Mag();
    _true_visible_energy = total_visible_momentum.E();

    auto const& mcp_h = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcParticleMap;
    for (size_t i = 0; i < mcp_h->size(); ++i) {
        mcParticleMap[mcp_h->at(i).TrackId()] = art::Ptr<simb::MCParticle>(mcp_h, i);
    }

    // Simple, threshold-free counting of final-state primary particles by PDG.
    for (int i = 0; i < mct.NParticles(); ++i) {
        const simb::MCParticle& particle = mct.GetParticle(i);
        if (!(particle.StatusCode() == 1 && particle.Process() == "primary")) continue;
        const int pdg = particle.PdgCode();
        switch (pdg) {
            case 13:   ++_count_mu_minus;   break;
            case -13:  ++_count_mu_plus;    break;
            case 11:   ++_count_e_minus;    break;
            case -11:  ++_count_e_plus;     break;
            case 22:   ++_count_gamma;      break;
            case 111:  ++_count_pi_zero;    break;
            case 211:  ++_count_pi_plus;    break;
            case -211: ++_count_pi_minus;   break;
            case 321:  ++_count_kaon_plus;  break;
            case -321: ++_count_kaon_minus; break;
            case 2212: ++_count_proton;     break;
            case 2112: ++_count_neutron;    break;
            case 3222: ++_count_sigma_plus; break;
            case 3212: ++_count_sigma_zero; break;
            case 3112: ++_count_sigma_minus;break;
            default:
                if (pdg == 130 || pdg == 310 || pdg == 311) {
                    ++_count_kaon_zero;
                } else if (std::abs(pdg) == 3122) {
                    ++_count_lambda; // count Λ and anti-Λ
                }
                break;
        }
    }

    for (size_t p = 0; p < mcp_h->size(); p++) {
        auto const& mcp = mcp_h->at(p);
        if (!(mcp.Process() == "primary" && mcp.StatusCode() == 1)) continue;
        
        const art::Ptr<simb::MCParticle> mcp_ptr(mcp_h, p);
        int trackid = mcp.TrackId();
        
        _mc_particle_trackid.push_back(trackid);
        _mc_particle_energy.push_back(mcp.E());
        _mc_particle_pdg.push_back(mcp.PdgCode());

        auto n_elastic = 0u;
        auto n_inelastic = 0u;
        art::Ptr<simb::MCParticle> final_scattered_particle;
        common::GetNScatters(mcp_h, mcp_ptr, final_scattered_particle, n_elastic, n_inelastic);
        
        _mc_elastic_scatters.push_back(n_elastic);
        _mc_inelastic_scatters.push_back(n_inelastic);
        _mc_momentum_x.push_back(mcp.Px());
        _mc_momentum_y.push_back(mcp.Py());
        _mc_momentum_z.push_back(mcp.Pz());
        _mc_end_momentum.push_back(final_scattered_particle->Momentum(std::max(0u, final_scattered_particle->NumberTrajectoryPoints() - 2)).Vect().Mag());
        _mc_start_vertex_x.push_back(mcp.Vx());
        _mc_start_vertex_y.push_back(mcp.Vy());
        _mc_start_vertex_z.push_back(mcp.Vz());
        _mc_end_vertex_x.push_back(mcp.EndX());
        _mc_end_vertex_y.push_back(mcp.EndY());
        _mc_end_vertex_z.push_back(mcp.EndZ());
        _mc_completeness.push_back(std::numeric_limits<float>::quiet_NaN());
        _mc_purity.push_back(std::numeric_limits<float>::quiet_NaN());

        std::string final_state;
        if (mcp.NumberDaughters() > 0) {
            bool has_same_daughter = false;
            for (int i = 0; i < mcp.NumberDaughters(); ++i) {
                if (mcParticleMap.count(mcp.Daughter(i))) {
                    const auto& daughter = *(mcParticleMap.at(mcp.Daughter(i)));
                    if (daughter.PdgCode() == mcp.PdgCode()) {
                        has_same_daughter = true;
                        break;
                    }
                }
            }
            if (mcp.EndProcess() == "Decay" && !has_same_daughter) {
                final_state = "decayed";
            } else {
                final_state = "interacted";
            }
        } else if (_mc_end_momentum.back() < 1e-3) {
            final_state = "stopped";
        } else {
            final_state = "escaped";
        }
        _mc_particle_final_state.push_back(final_state);

        std::vector<int> daughter_pdgs;
        std::vector<float> daughter_energies;
        std::vector<std::string> daughter_processes;
        std::vector<float> daughter_moms_x, daughter_moms_y, daughter_moms_z;
        std::vector<float> daughter_vtxs_x, daughter_vtxs_y, daughter_vtxs_z;
        
        _mc_daughter_process_idx.push_back(_mc_daughter_process_flat.size());
        for (int i = 0; i < mcp.NumberDaughters(); ++i) {
            if (mcParticleMap.count(mcp.Daughter(i))) {
                const auto& daughter = *(mcParticleMap.at(mcp.Daughter(i)));
                daughter_pdgs.push_back(daughter.PdgCode());
                daughter_energies.push_back(daughter.E());
                daughter_processes.push_back(daughter.Process());
                _mc_daughter_process_flat.push_back(daughter.Process());
                daughter_moms_x.push_back(daughter.Px());
                daughter_moms_y.push_back(daughter.Py());
                daughter_moms_z.push_back(daughter.Pz());
                daughter_vtxs_x.push_back(daughter.Vx());
                daughter_vtxs_y.push_back(daughter.Vy());
                daughter_vtxs_z.push_back(daughter.Vz());
            }
        }
        _mc_daughter_pdg.push_back(daughter_pdgs);
        _mc_daughter_energy.push_back(daughter_energies);
        _mc_daughter_mom_x.push_back(daughter_moms_x);
        _mc_daughter_mom_y.push_back(daughter_moms_y);
        _mc_daughter_mom_z.push_back(daughter_moms_z);
        _mc_daughter_vtx_x.push_back(daughter_vtxs_x);
        _mc_daughter_vtx_y.push_back(daughter_vtxs_y);
        _mc_daughter_vtx_z.push_back(daughter_vtxs_z);

        // NOTE: recursive all-chain collection removed.
    }
}

DEFINE_ART_CLASS_TOOL(TruthAnalysis)
}

#endif
