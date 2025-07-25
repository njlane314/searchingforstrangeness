#ifndef ANALYSIS_TRUTH_CXX
#define ANALYSIS_TRUTH_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraSliceIdHelper.h"
#include "TTree.h"
#include <limits>
#include <vector>
#include <string>
#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Containment.h"
#include "../CommonDefs/Scatters.h"

namespace analysis {

class TruthAnalysis : public AnalysisToolBase {
public:
    explicit TruthAnalysis(fhicl::ParameterSet const& p);
    virtual ~TruthAnalysis() = default;
    void configure(const fhicl::ParameterSet& pset) override;
    void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;
    void analyseEvent(const art::Event& event, bool is_data) override;
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    art::InputTag fMCTproducer;
    art::InputTag fMCFLUXproducer;
    art::InputTag fMCPproducer;
    art::InputTag fHITproducer;
    art::InputTag fHITTASproducer;
    art::InputTag fMCSHRproducer;
    art::InputTag fMCTRKproducer;
    art::InputTag fCLSproducer;
    float fFidvolXstart;
    float fFidvolXend;
    float fFidvolYstart;
    float fFidvolYend;
    float fFidvolZstart;
    float fFidvolZend;
    float fProtonKEThreshold;
    float fPionKEThreshold;
    float fMuonKEThreshold;
    float fElectronKEThreshold;
    float fGammaKEThreshold;
    float fNeutralPionKEThreshold;
    float fHadronKEThreshold;

    struct ParticleCounter {
        std::vector<int> pdgs;
        bool use_abs;
        float threshold;
        int* counter;
    };
    std::vector<ParticleCounter> particle_counters;

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
    std::vector<std::vector<std::string>> _mc_daughter_process;
    std::vector<std::vector<float>> _mc_daughter_mom_x;
    std::vector<std::vector<float>> _mc_daughter_mom_y;
    std::vector<std::vector<float>> _mc_daughter_mom_z;
    std::vector<std::vector<float>> _mc_daughter_vtx_x;
    std::vector<std::vector<float>> _mc_daughter_vtx_y;
    std::vector<std::vector<float>> _mc_daughter_vtx_z;
    std::vector<int> _mc_allchain_primary_index;
    std::vector<int> _mc_allchain_trackid;
    std::vector<int> _mc_allchain_pdg;
    std::vector<float> _mc_allchain_energy;
    std::vector<uint> _mc_allchain_elastic_scatters;
    std::vector<uint> _mc_allchain_inelastic_scatters;
    std::vector<float> _mc_allchain_momentum_x;
    std::vector<float> _mc_allchain_momentum_y;
    std::vector<float> _mc_allchain_momentum_z;
    std::vector<float> _mc_allchain_end_momentum;
    std::vector<float> _mc_allchain_start_vertex_x;
    std::vector<float> _mc_allchain_start_vertex_y;
    std::vector<float> _mc_allchain_start_vertex_z;
    std::vector<float> _mc_allchain_end_vertex_x;
    std::vector<float> _mc_allchain_end_vertex_y;
    std::vector<float> _mc_allchain_end_vertex_z;
    std::vector<int> _mc_allchain_parent_trackid;
    std::vector<std::string> _mc_allchain_process;
    std::vector<std::string> _mc_allchain_final_state;
    std::vector<float> _mc_allchain_completeness;
    std::vector<float> _mc_allchain_purity;
    float _true_transverse_momentum;
    float _true_visible_transverse_momentum;
    float _true_total_momentum;
    float _true_visible_total_momentum;
    float _true_visible_energy;
    void CollectDescendants(const art::ValidHandle<std::vector<simb::MCParticle>>& mcp_h, const std::map<int, art::Ptr<simb::MCParticle>>& mcParticleMap, const art::Ptr<simb::MCParticle>& part, int primary_index);
};

TruthAnalysis::TruthAnalysis(fhicl::ParameterSet const& p) {
    this->configure(p);
}

void TruthAnalysis::configure(fhicl::ParameterSet const& p) {
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCFLUXproducer = p.get<art::InputTag>("MCFLUXproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fHITTASproducer = p.get<art::InputTag>("HITTASproducer");
    fMCSHRproducer = p.get<art::InputTag>("MCSHRproducer");
    fMCTRKproducer = p.get<art::InputTag>("MCTRKproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fFidvolXstart = p.get<double>("fidvolXstart", 10);
    fFidvolXend = p.get<double>("fidvolXend", 10);
    fFidvolYstart = p.get<double>("fidvolYstart", 15);
    fFidvolYend = p.get<double>("fidvolYend", 15);
    fFidvolZstart = p.get<double>("fidvolZstart", 10);
    fFidvolZend = p.get<double>("fidvolZend", 50);
    fProtonKEThreshold = p.get<float>("ProtonKEThreshold", 0.035);
    fPionKEThreshold = p.get<float>("PionKEThreshold", 0.020);
    fMuonKEThreshold = p.get<float>("MuonKEThreshold", 0.015);
    fElectronKEThreshold = p.get<float>("ElectronKEThreshold", 0.030);
    fGammaKEThreshold = p.get<float>("GammaKEThreshold", 0.030);
    fNeutralPionKEThreshold = p.get<float>("NeutralPionKEThreshold", 0.030);
    fHadronKEThreshold = p.get<float>("HadronKEThreshold", 0.035);
    particle_counters = {
        {{13}, false, fMuonKEThreshold, &_count_mu_minus},
        {{-13}, false, fMuonKEThreshold, &_count_mu_plus},
        {{11}, false, fElectronKEThreshold, &_count_e_minus},
        {{-11}, false, fElectronKEThreshold, &_count_e_plus},
        {{22}, false, fGammaKEThreshold, &_count_gamma},
        {{111}, false, fNeutralPionKEThreshold, &_count_pi_zero},
        {{211}, false, fPionKEThreshold, &_count_pi_plus},
        {{-211}, false, fPionKEThreshold, &_count_pi_minus},
        {{321}, false, fHadronKEThreshold, &_count_kaon_plus},
        {{-321}, false, fHadronKEThreshold, &_count_kaon_minus},
        {{130, 310, 311}, false, -1.f, &_count_kaon_zero},
        {{2212}, false, fProtonKEThreshold, &_count_proton},
        {{2112}, false, -1.f, &_count_neutron},
        {{3122}, true, -1.f, &_count_lambda},
        {{3222}, false, fHadronKEThreshold, &_count_sigma_plus},
        {{3212}, false, -1.f, &_count_sigma_zero},
        {{3112}, false, fHadronKEThreshold, &_count_sigma_minus}
    };
}

void TruthAnalysis::setBranches(TTree* _tree) {
    _tree->Branch("neutrino_pdg", &_neutrino_pdg, "neutrino_pdg/I");
    _tree->Branch("interaction_ccnc", &_interaction_ccnc, "interaction_ccnc/I");
    _tree->Branch("interaction_mode", &_interaction_mode, "interaction_mode/I");
    _tree->Branch("interaction_type", &_interaction_type, "interaction_type/I");
    _tree->Branch("neutrino_energy", &_neutrino_energy, "neutrino_energy/F");
    _tree->Branch("neutrino_theta", &_neutrino_theta, "neutrino_theta/F");
    _tree->Branch("neutrino_pt", &_neutrino_pt, "neutrino_pt/F");
    _tree->Branch("target_nucleus_pdg", &_target_nucleus_pdg, "target_nucleus_pdg/I");
    _tree->Branch("hit_nucleon_pdg", &_hit_nucleon_pdg, "hit_nucleon_pdg/I");
    _tree->Branch("kinematic_W", &_kinematic_W, "kinematic_W/F");
    _tree->Branch("kinematic_X", &_kinematic_X, "kinematic_X/F");
    _tree->Branch("kinematic_Y", &_kinematic_Y, "kinematic_Y/F");
    _tree->Branch("kinematic_Q_squared", &_kinematic_Q_squared, "kinematic_Q_squared/F");
    _tree->Branch("neutrino_momentum_x", &_neutrino_momentum_x, "neutrino_momentum_x/F");
    _tree->Branch("neutrino_momentum_y", &_neutrino_momentum_y, "neutrino_momentum_y/F");
    _tree->Branch("neutrino_momentum_z", &_neutrino_momentum_z, "neutrino_momentum_z/F");
    _tree->Branch("neutrino_vertex_x", &_neutrino_vertex_x, "neutrino_vertex_x/F");
    _tree->Branch("neutrino_vertex_y", &_neutrino_vertex_y, "neutrino_vertex_y/F");
    _tree->Branch("neutrino_vertex_z", &_neutrino_vertex_z, "neutrino_vertex_z/F");
    _tree->Branch("neutrino_vertex_wire_u", &_neutrino_vertex_wire_u, "neutrino_vertex_wire_u/F");
    _tree->Branch("neutrino_vertex_wire_v", &_neutrino_vertex_wire_v, "neutrino_vertex_wire_v/F");
    _tree->Branch("neutrino_vertex_wire_w", &_neutrino_vertex_wire_w, "neutrino_vertex_wire_w/F");
    _tree->Branch("neutrino_vertex_time", &_neutrino_vertex_time, "neutrino_vertex_time/F");
    _tree->Branch("neutrino_sce_vertex_x", &_neutrino_sce_vertex_x, "neutrino_sce_vertex_x/F");
    _tree->Branch("neutrino_sce_vertex_y", &_neutrino_sce_vertex_y, "neutrino_sce_vertex_y/F");
    _tree->Branch("neutrino_sce_vertex_z", &_neutrino_sce_vertex_z, "neutrino_sce_vertex_z/F");
    _tree->Branch("lepton_energy", &_lepton_energy, "lepton_energy/F");
    _tree->Branch("true_neutrino_momentum_x", &_true_neutrino_momentum_x, "true_neutrino_momentum_x/F");
    _tree->Branch("true_neutrino_momentum_y", &_true_neutrino_momentum_y, "true_neutrino_momentum_y/F");
    _tree->Branch("true_neutrino_momentum_z", &_true_neutrino_momentum_z, "true_neutrino_momentum_z/F");
    _tree->Branch("flux_path_length", &_flux_path_length, "flux_path_length/F");
    _tree->Branch("flux_parent_pdg", &_flux_parent_pdg, "flux_parent_pdg/I");
    _tree->Branch("flux_hadron_pdg", &_flux_hadron_pdg, "flux_hadron_pdg/I");
    _tree->Branch("flux_decay_mode", &_flux_decay_mode, "flux_decay_mode/I");
    _tree->Branch("flux_decay_vtx_x", &_flux_decay_vtx_x, "flux_decay_vtx_x/D");
    _tree->Branch("flux_decay_vtx_y", &_flux_decay_vtx_y, "flux_decay_vtx_y/D");
    _tree->Branch("flux_decay_vtx_z", &_flux_decay_vtx_z, "flux_decay_vtx_z/D");
    _tree->Branch("flux_decay_mom_x", &_flux_decay_mom_x, "flux_decay_mom_x/D");
    _tree->Branch("flux_decay_mom_y", &_flux_decay_mom_y, "flux_decay_mom_y/D");
    _tree->Branch("flux_decay_mom_z", &_flux_decay_mom_z, "flux_decay_mom_z/D");
    _tree->Branch("numi_baseline", &_numi_baseline, "numi_baseline/D");
    _tree->Branch("numi_off_axis_angle", &_numi_off_axis_angle, "numi_off_axis_angle/F");
    _tree->Branch("bnb_baseline", &_bnb_baseline, "bnb_baseline/D");
    _tree->Branch("bnb_off_axis_angle", &_bnb_off_axis_angle, "bnb_off_axis_angle/F");
    _tree->Branch("is_vertex_in_fiducial", &_is_vertex_in_fiducial, "is_vertex_in_fiducial/O");
    _tree->Branch("count_mu_minus", &_count_mu_minus, "count_mu_minus/I");
    _tree->Branch("count_mu_plus", &_count_mu_plus, "count_mu_plus/I");
    _tree->Branch("count_e_minus", &_count_e_minus, "count_e_minus/I");
    _tree->Branch("count_e_plus", &_count_e_plus, "count_e_plus/I");
    _tree->Branch("count_pi_zero", &_count_pi_zero, "count_pi_zero/I");
    _tree->Branch("count_pi_plus", &_count_pi_plus, "count_pi_plus/I");
    _tree->Branch("count_pi_minus", &_count_pi_minus, "count_pi_minus/I");
    _tree->Branch("count_kaon_plus", &_count_kaon_plus, "count_kaon_plus/I");
    _tree->Branch("count_kaon_minus", &_count_kaon_minus, "count_kaon_minus/I");
    _tree->Branch("count_kaon_zero", &_count_kaon_zero, "count_kaon_zero/I");
    _tree->Branch("count_proton", &_count_proton, "count_proton/I");
    _tree->Branch("count_neutron", &_count_neutron, "count_neutron/I");
    _tree->Branch("count_gamma", &_count_gamma, "count_gamma/I");
    _tree->Branch("count_lambda", &_count_lambda, "count_lambda/I");
    _tree->Branch("count_sigma_plus", &_count_sigma_plus, "count_sigma_plus/I");
    _tree->Branch("count_sigma_zero", &_count_sigma_zero, "count_sigma_zero/I");
    _tree->Branch("count_sigma_minus", &_count_sigma_minus, "count_sigma_minus/I");
    _tree->Branch("mc_particle_pdg", "std::vector<int>", &_mc_particle_pdg);
    _tree->Branch("mc_particle_trackid", "std::vector<int>", &_mc_particle_trackid);
    _tree->Branch("mc_particle_energy", "std::vector<float>", &_mc_particle_energy);
    _tree->Branch("mc_elastic_scatters", "std::vector<uint>", &_mc_elastic_scatters);
    _tree->Branch("mc_inelastic_scatters", "std::vector<uint>", &_mc_inelastic_scatters);
    _tree->Branch("mc_momentum_x", "std::vector<float>", &_mc_momentum_x);
    _tree->Branch("mc_momentum_y", "std::vector<float>", &_mc_momentum_y);
    _tree->Branch("mc_momentum_z", "std::vector<float>", &_mc_momentum_z);
    _tree->Branch("mc_end_momentum", "std::vector<float>", &_mc_end_momentum);
    _tree->Branch("mc_start_vertex_x", "std::vector<float>", &_mc_start_vertex_x);
    _tree->Branch("mc_start_vertex_y", "std::vector<float>", &_mc_start_vertex_y);
    _tree->Branch("mc_start_vertex_z", "std::vector<float>", &_mc_start_vertex_z);
    _tree->Branch("mc_end_vertex_x", "std::vector<float>", &_mc_end_vertex_x);
    _tree->Branch("mc_end_vertex_y", "std::vector<float>", &_mc_end_vertex_y);
    _tree->Branch("mc_end_vertex_z", "std::vector<float>", &_mc_end_vertex_z);
    _tree->Branch("mc_particle_final_state", "std::vector<std::string>", &_mc_particle_final_state);
    _tree->Branch("mc_completeness", "std::vector<float>", &_mc_completeness);
    _tree->Branch("mc_purity", "std::vector<float>", &_mc_purity);
    _tree->Branch("mc_daughter_pdg", "std::vector<std::vector<int>>", &_mc_daughter_pdg);
    _tree->Branch("mc_daughter_energy", "std::vector<std::vector<float>>", &_mc_daughter_energy);
    _tree->Branch("mc_daughter_process", "std::vector<std::vector<std::string>>", &_mc_daughter_process);
    _tree->Branch("mc_daughter_mom_x", "std::vector<std::vector<float>>", &_mc_daughter_mom_x);
    _tree->Branch("mc_daughter_mom_y", "std::vector<std::vector<float>>", &_mc_daughter_mom_y);
    _tree->Branch("mc_daughter_mom_z", "std::vector<std::vector<float>>", &_mc_daughter_mom_z);
    _tree->Branch("mc_daughter_vtx_x", "std::vector<std::vector<float>>", &_mc_daughter_vtx_x);
    _tree->Branch("mc_daughter_vtx_y", "std::vector<std::vector<float>>", &_mc_daughter_vtx_y);
    _tree->Branch("mc_daughter_vtx_z", "std::vector<std::vector<float>>", &_mc_daughter_vtx_z);
    _tree->Branch("mc_allchain_primary_index", "std::vector<int>", &_mc_allchain_primary_index);
    _tree->Branch("mc_allchain_trackid", "std::vector<int>", &_mc_allchain_trackid);
    _tree->Branch("mc_allchain_pdg", "std::vector<int>", &_mc_allchain_pdg);
    _tree->Branch("mc_allchain_energy", "std::vector<float>", &_mc_allchain_energy);
    _tree->Branch("mc_allchain_elastic_scatters", "std::vector<uint>", &_mc_allchain_elastic_scatters);
    _tree->Branch("mc_allchain_inelastic_scatters", "std::vector<uint>", &_mc_allchain_inelastic_scatters);
    _tree->Branch("mc_allchain_momentum_x", "std::vector<float>", &_mc_allchain_momentum_x);
    _tree->Branch("mc_allchain_momentum_y", "std::vector<float>", &_mc_allchain_momentum_y);
    _tree->Branch("mc_allchain_momentum_z", "std::vector<float>", &_mc_allchain_momentum_z);
    _tree->Branch("mc_allchain_end_momentum", "std::vector<float>", &_mc_allchain_end_momentum);
    _tree->Branch("mc_allchain_start_vertex_x", "std::vector<float>", &_mc_allchain_start_vertex_x);
    _tree->Branch("mc_allchain_start_vertex_y", "std::vector<float>", &_mc_allchain_start_vertex_y);
    _tree->Branch("mc_allchain_start_vertex_z", "std::vector<float>", &_mc_allchain_start_vertex_z);
    _tree->Branch("mc_allchain_end_vertex_x", "std::vector<float>", &_mc_allchain_end_vertex_x);
    _tree->Branch("mc_allchain_end_vertex_y", "std::vector<float>", &_mc_allchain_end_vertex_y);
    _tree->Branch("mc_allchain_end_vertex_z", "std::vector<float>", &_mc_allchain_end_vertex_z);
    _tree->Branch("mc_allchain_parent_trackid", "std::vector<int>", &_mc_allchain_parent_trackid);
    _tree->Branch("mc_allchain_process", "std::vector<std::string>", &_mc_allchain_process);
    _tree->Branch("mc_allchain_final_state", "std::vector<std::string>", &_mc_allchain_final_state);
    _tree->Branch("mc_allchain_completeness", "std::vector<float>", &_mc_allchain_completeness);
    _tree->Branch("mc_allchain_purity", "std::vector<float>", &_mc_allchain_purity);
    _tree->Branch("true_transverse_momentum", &_true_transverse_momentum, "true_transverse_momentum/F");
    _tree->Branch("true_visible_transverse_momentum", &_true_visible_transverse_momentum, "true_visible_transverse_momentum/F");
    _tree->Branch("true_total_momentum", &_true_total_momentum, "true_total_momentum/F");
    _tree->Branch("true_visible_total_momentum", &_true_visible_total_momentum, "true_visible_total_momentum/F");
    _tree->Branch("true_visible_energy", &_true_visible_energy, "true_visible_energy/F");
}

void TruthAnalysis::resetTTree(TTree* tree) {
    _neutrino_pdg = std::numeric_limits<int>::lowest();
    _interaction_ccnc = std::numeric_limits<int>::lowest();
    _interaction_mode = std::numeric_limits<int>::lowest();
    _interaction_type = std::numeric_limits<int>::lowest();
    _neutrino_energy = std::numeric_limits<float>::lowest();
    _neutrino_theta = std::numeric_limits<float>::lowest();
    _neutrino_pt = std::numeric_limits<float>::lowest();
    _target_nucleus_pdg = std::numeric_limits<int>::lowest();
    _hit_nucleon_pdg = std::numeric_limits<int>::lowest();
    _kinematic_W = std::numeric_limits<float>::lowest();
    _kinematic_X = std::numeric_limits<float>::lowest();
    _kinematic_Y = std::numeric_limits<float>::lowest();
    _kinematic_Q_squared = std::numeric_limits<float>::lowest();
    _neutrino_momentum_x = std::numeric_limits<float>::lowest();
    _neutrino_momentum_y = std::numeric_limits<float>::lowest();
    _neutrino_momentum_z = std::numeric_limits<float>::lowest();
    _neutrino_vertex_x = std::numeric_limits<float>::lowest();
    _neutrino_vertex_y = std::numeric_limits<float>::lowest();
    _neutrino_vertex_z = std::numeric_limits<float>::lowest();
    _neutrino_vertex_wire_u = std::numeric_limits<float>::lowest();
    _neutrino_vertex_wire_v = std::numeric_limits<float>::lowest();
    _neutrino_vertex_wire_w = std::numeric_limits<float>::lowest();
    _neutrino_vertex_time = std::numeric_limits<float>::lowest();
    _neutrino_sce_vertex_x = std::numeric_limits<float>::lowest();
    _neutrino_sce_vertex_y = std::numeric_limits<float>::lowest();
    _neutrino_sce_vertex_z = std::numeric_limits<float>::lowest();
    _lepton_energy = std::numeric_limits<float>::lowest();
    _true_neutrino_momentum_x = std::numeric_limits<float>::lowest();
    _true_neutrino_momentum_y = std::numeric_limits<float>::lowest();
    _true_neutrino_momentum_z = std::numeric_limits<float>::lowest();
    _flux_path_length = std::numeric_limits<float>::lowest();
    _flux_parent_pdg = std::numeric_limits<int>::lowest();
    _flux_hadron_pdg = std::numeric_limits<int>::lowest();
    _flux_decay_mode = std::numeric_limits<int>::lowest();
    _flux_decay_vtx_x = std::numeric_limits<double>::lowest();
    _flux_decay_vtx_y = std::numeric_limits<double>::lowest();
    _flux_decay_vtx_z = std::numeric_limits<double>::lowest();
    _flux_decay_mom_x = std::numeric_limits<double>::lowest();
    _flux_decay_mom_y = std::numeric_limits<double>::lowest();
    _flux_decay_mom_z = std::numeric_limits<double>::lowest();
    _numi_baseline = std::numeric_limits<double>::lowest();
    _numi_off_axis_angle = -999.0;
    _bnb_baseline = std::numeric_limits<double>::lowest();
    _bnb_off_axis_angle = -999.0;
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
    _mc_daughter_process.clear();
    _mc_daughter_mom_x.clear();
    _mc_daughter_mom_y.clear();
    _mc_daughter_mom_z.clear();
    _mc_daughter_vtx_x.clear();
    _mc_daughter_vtx_y.clear();
    _mc_daughter_vtx_z.clear();
    _mc_allchain_primary_index.clear();
    _mc_allchain_trackid.clear();
    _mc_allchain_pdg.clear();
    _mc_allchain_energy.clear();
    _mc_allchain_elastic_scatters.clear();
    _mc_allchain_inelastic_scatters.clear();
    _mc_allchain_momentum_x.clear();
    _mc_allchain_momentum_y.clear();
    _mc_allchain_momentum_z.clear();
    _mc_allchain_end_momentum.clear();
    _mc_allchain_start_vertex_x.clear();
    _mc_allchain_start_vertex_y.clear();
    _mc_allchain_start_vertex_z.clear();
    _mc_allchain_end_vertex_x.clear();
    _mc_allchain_end_vertex_y.clear();
    _mc_allchain_end_vertex_z.clear();
    _mc_allchain_parent_trackid.clear();
    _mc_allchain_process.clear();
    _mc_allchain_final_state.clear();
    _mc_allchain_completeness.clear();
    _mc_allchain_purity.clear();
    _true_transverse_momentum = 0;
    _true_visible_transverse_momentum = 0;
    _true_total_momentum = 0;
    _true_visible_total_momentum = 0;
    _true_visible_energy = 0;
}

void TruthAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
    if (is_data) return;

    auto const& cluster_h = event.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
    art::FindManyP<recob::Hit> assocHits(cluster_h, event, fCLSproducer);

    std::vector<art::Ptr<recob::Hit>> inputHits;
    for (const auto& pfp : slice_pfp_v) {
        auto clusters = pfp.get<recob::Cluster>();
        for (const auto& cluster : clusters) {
            auto hits = assocHits.at(cluster.key());
            inputHits.insert(inputHits.end(), hits.begin(), hits.end());
        }
    }

    auto const& hit_h = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    auto assocMCPart = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_h, event, fHITTASproducer);

    auto const& mcshower_h = event.getValidHandle<std::vector<sim::MCShower>>(fMCSHRproducer);
    auto const& mctrack_h = event.getValidHandle<std::vector<sim::MCTrack>>(fMCTRKproducer);

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

    for (size_t j = 0; j < _mc_allchain_pdg.size(); ++j) {
        float purity = 0.0;
        float completeness = 0.0;
        if (computeMetrics(_mc_allchain_trackid[j], purity, completeness)) {
            _mc_allchain_purity[j] = purity;
            _mc_allchain_completeness[j] = completeness;
        } else {
            _mc_allchain_purity[j] = 0.0;
            _mc_allchain_completeness[j] = 0.0;
        }
    }
}

void TruthAnalysis::analyseEvent(art::Event const& e, bool is_data) {
    if (is_data) {
        this->resetTTree(nullptr);
        return;
    }

    auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    auto mct = mct_h->at(0);

    auto const& mcflux_h = e.getValidHandle<std::vector<simb::MCFlux>>(fMCFLUXproducer);
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
        if (particle.StatusCode() != 1) continue;
        total_momentum += particle.Momentum(0);
        total_visible_momentum += particle.Momentum(0);
    }
    _true_transverse_momentum = total_momentum.Perp();
    _true_visible_transverse_momentum = total_visible_momentum.Perp();
    _true_total_momentum = total_momentum.Mag();
    _true_visible_total_momentum = total_visible_momentum.Mag();
    _true_visible_energy = total_visible_momentum.E();

    auto const& mcp_h = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mcParticleMap;
    for (size_t i = 0; i < mcp_h->size(); ++i) {
        mcParticleMap[mcp_h->at(i).TrackId()] = art::Ptr<simb::MCParticle>(mcp_h, i);
    }

    for (int i = 0; i < mct.NParticles(); ++i) {
        const simb::MCParticle& particle = mct.GetParticle(i);
        if (particle.StatusCode() != 1) continue;
        int pdg = particle.PdgCode();

        for (const auto& pc : particle_counters) {
            int check_pdg = pc.use_abs ? std::abs(pdg) : pdg;
            auto it = std::find(pc.pdgs.begin(), pc.pdgs.end(), check_pdg);
            if (it == pc.pdgs.end()) continue;

            float ke = particle.E() - particle.Mass();
            if (pc.threshold < 0.f || ke > pc.threshold) {
                (*pc.counter)++;
                break;
            }
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
        _mc_completeness.push_back(std::numeric_limits<float>::lowest());
        _mc_purity.push_back(std::numeric_limits<float>::lowest());

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
        
        for (int i = 0; i < mcp.NumberDaughters(); ++i) {
            if (mcParticleMap.count(mcp.Daughter(i))) {
                const auto& daughter = *(mcParticleMap.at(mcp.Daughter(i)));
                daughter_pdgs.push_back(daughter.PdgCode());
                daughter_energies.push_back(daughter.E());
                daughter_processes.push_back(daughter.Process());
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
        _mc_daughter_process.push_back(daughter_processes);
        _mc_daughter_mom_x.push_back(daughter_moms_x);
        _mc_daughter_mom_y.push_back(daughter_moms_y);
        _mc_daughter_mom_z.push_back(daughter_moms_z);
        _mc_daughter_vtx_x.push_back(daughter_vtxs_x);
        _mc_daughter_vtx_y.push_back(daughter_vtxs_y);
        _mc_daughter_vtx_z.push_back(daughter_vtxs_z);

        CollectDescendants(mcp_h, mcParticleMap, mcp_ptr, p);
    }
}

void TruthAnalysis::CollectDescendants(const art::ValidHandle<std::vector<simb::MCParticle>>& mcp_h, const std::map<int, art::Ptr<simb::MCParticle>>& mcParticleMap, const art::Ptr<simb::MCParticle>& part, int primary_index) {
    for (int i = 0; i < part->NumberDaughters(); ++i) {
        int dau_tid = part->Daughter(i);
        if (!mcParticleMap.count(dau_tid)) continue;

        auto const& dau = mcParticleMap.at(dau_tid);

        _mc_allchain_primary_index.push_back(primary_index);
        _mc_allchain_trackid.push_back(dau_tid);
        _mc_allchain_pdg.push_back(dau->PdgCode());
        _mc_allchain_energy.push_back(dau->E());
        _mc_allchain_momentum_x.push_back(dau->Px());
        _mc_allchain_momentum_y.push_back(dau->Py());
        _mc_allchain_momentum_z.push_back(dau->Pz());
        _mc_allchain_start_vertex_x.push_back(dau->Vx());
        _mc_allchain_start_vertex_y.push_back(dau->Vy());
        _mc_allchain_start_vertex_z.push_back(dau->Vz());
        _mc_allchain_end_vertex_x.push_back(dau->EndX());
        _mc_allchain_end_vertex_y.push_back(dau->EndY());
        _mc_allchain_end_vertex_z.push_back(dau->EndZ());
        _mc_allchain_parent_trackid.push_back(part->TrackId());
        _mc_allchain_process.push_back(dau->Process());

        auto n_elastic = 0u;
        auto n_inelastic = 0u;
        art::Ptr<simb::MCParticle> final_scattered;
        common::GetNScatters(mcp_h, dau, final_scattered, n_elastic, n_inelastic);
        _mc_allchain_elastic_scatters.push_back(n_elastic);
        _mc_allchain_inelastic_scatters.push_back(n_inelastic);
        _mc_allchain_end_momentum.push_back(final_scattered->Momentum(std::max(0u, final_scattered->NumberTrajectoryPoints() - 2)).Vect().Mag());

        std::string final_state;
        if (dau->NumberDaughters() > 0) {
            bool has_same_daughter = false;
            for (int j = 0; j < dau->NumberDaughters(); ++j) {
                if (mcParticleMap.count(dau->Daughter(j))) {
                    const auto& grand_daughter = *(mcParticleMap.at(dau->Daughter(j)));
                    if (grand_daughter.PdgCode() == dau->PdgCode()) {
                        has_same_daughter = true;
                        break;
                    }
                }
            }
            if (dau->EndProcess() == "Decay" && !has_same_daughter) {
                final_state = "decayed";
            } else {
                final_state = "interacted";
            }
        } else if (_mc_allchain_end_momentum.back() < 1e-3) {
            final_state = "stopped";
        } else {
            final_state = "escaped";
        }
        _mc_allchain_final_state.push_back(final_state);

        _mc_allchain_completeness.push_back(std::numeric_limits<float>::lowest());
        _mc_allchain_purity.push_back(std::numeric_limits<float>::lowest());

        CollectDescendants(mcp_h, mcParticleMap, dau, primary_index);
    }
}

DEFINE_ART_CLASS_TOOL(TruthAnalysis)
}

#endif 