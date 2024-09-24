#ifndef ANALYSIS_TRACKCALORIMETRY_CXX
#define ANALYSIS_TRACKCALORIMETRY_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "CommonFunctions/Types.h"
#include "CommonFunctions/Backtracking.h"
#include "CommonFunctions/Scores.h"
#include "CommonFunctions/Identification.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Geometry.h"
#include "CommonFunctions/Calibration.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace analysis
{

class TrackAnalysis : public AnalysisToolBase
{

public:
    TrackAnalysis(const fhicl::ParameterSet &pset);
    ~TrackAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);
    void analyzeEvent(art::Event const &e, bool is_data) override;
    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected) override;
    void SaveTruth(art::Event const &e);
    void fillDefault();
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    float CalculateTrackTrunkdEdxByHits(const std::vector<float> &dEdx_per_hit);
    float CalculateTrackTrunkdEdxByRange(const std::vector<float> &dEdx_per_hit, const std::vector<float> &rr_per_hit);
    void CalculateTrackDeflections(const art::Ptr<recob::Track> &trk, std::vector<float> &mean_v, std::vector<float> &stdev_v, std::vector<float> &separation_mean_v);

    const trkf::TrackMomentumCalculator _trkmom;
    const trkf::TrajectoryMCSFitter _mcsfitter;

    TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
    TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);

    art::InputTag _CALOproducer;
    art::InputTag _PIDproducer;
    art::InputTag _TRKproducer;
    art::InputTag _BacktrackTag;
    art::InputTag _Hproducer;
    art::InputTag _CLSproducer;

    bool _RecalibrateHits;
    float _EnergyThresholdForHits;
    std::vector<float> _ADCtoE; 
    float _EndSpacepointDistance;

    std::vector<size_t> _trk_pfp_id_v;
    std::vector<float> _trk_score_v;

    std::vector<float> _trk_start_x_v;
    std::vector<float> _trk_start_y_v;
    std::vector<float> _trk_start_z_v;

    std::vector<float> _trk_sce_start_x_v;
    std::vector<float> _trk_sce_start_y_v;
    std::vector<float> _trk_sce_start_z_v;

    std::vector<float> _trk_distance_v;

    std::vector<float> _trk_theta_v;
    std::vector<float> _trk_phi_v;

    std::vector<float> _trk_dir_x_v;
    std::vector<float> _trk_dir_y_v;
    std::vector<float> _trk_dir_z_v;

    std::vector<float> _trk_end_x_v;
    std::vector<float> _trk_end_y_v;
    std::vector<float> _trk_end_z_v;

    std::vector<float> _trk_sce_end_x_v;
    std::vector<float> _trk_sce_end_y_v;
    std::vector<float> _trk_sce_end_z_v;

    std::vector<float> _trk_len_v;

    std::vector<float> _trk_bragg_p_v; // Largest bragg PID value under the proton hypothesis between forward & backward fit in the w plane
    std::vector<float> _trk_bragg_mu_v; // ... under the muon hypothesis ...
    std::vector<float> _trk_bragg_pion_v; // ... under the pion hypothesis ...
    std::vector<float> _trk_bragg_mip_v; // Bragg PID value under the MIP hypothesis
    std::vector<float> _trk_bragg_p_alt_dir_v; // Bragg PID value for the alternative direction
    std::vector<float> _trk_bragg_mu_alt_dir_v; // Bragg PID value for the alternative direction
    std::vector<float> _trk_bragg_pion_alt_dir_v; // Bragg PID value for the alternative direction
    std::vector<bool> _trk_bragg_p_fwd_preferred_v; // Whether _trk_bragg_p_v uses the forward fit
    std::vector<bool> _trk_bragg_mu_fwd_preferred_v; // Whether _trk_bragg_mu_v uses the forward fit
    std::vector<bool> _trk_bragg_pion_fwd_preferred_v; // Whether _trk_bragg_pion_v uses the forward fit
    std::vector<float> _trk_pid_chipr_v;
    std::vector<float> _trk_pid_chika_v;
    std::vector<float> _trk_pid_chipi_v;
    std::vector<float> _trk_pid_chimu_v;
    std::vector<float> _trk_pida_v;

    std::vector<float> _trk_bragg_p_u_v; // Same as above but in the u plane
    std::vector<float> _trk_bragg_mu_u_v;
    std::vector<float> _trk_bragg_pion_u_v;
    std::vector<float> _trk_bragg_mip_u_v;
    std::vector<float> _trk_bragg_p_alt_dir_u_v;
    std::vector<float> _trk_bragg_mu_alt_dir_u_v;
    std::vector<float> _trk_bragg_pion_alt_dir_u_v;
    std::vector<bool> _trk_bragg_p_fwd_preferred_u_v;
    std::vector<bool> _trk_bragg_mu_fwd_preferred_u_v;
    std::vector<bool> _trk_bragg_pion_fwd_preferred_u_v;
    std::vector<float> _trk_pid_chipr_u_v;
    std::vector<float> _trk_pid_chika_u_v;
    std::vector<float> _trk_pid_chipi_u_v;
    std::vector<float> _trk_pid_chimu_u_v;
    std::vector<float> _trk_pida_u_v;

    std::vector<float> _trk_bragg_p_v_v; // Same as above but in the v plane
    std::vector<float> _trk_bragg_mu_v_v;
    std::vector<float> _trk_bragg_pion_v_v;
    std::vector<float> _trk_bragg_mip_v_v;
    std::vector<float> _trk_bragg_p_alt_dir_v_v;
    std::vector<float> _trk_bragg_mu_alt_dir_v_v;
    std::vector<float> _trk_bragg_pion_alt_dir_v_v;
    std::vector<bool> _trk_bragg_p_fwd_preferred_v_v;
    std::vector<bool> _trk_bragg_mu_fwd_preferred_v_v;
    std::vector<bool> _trk_bragg_pion_fwd_preferred_v_v;
    std::vector<float> _trk_pid_chipr_v_v;
    std::vector<float> _trk_pid_chika_v_v;
    std::vector<float> _trk_pid_chipi_v_v;
    std::vector<float> _trk_pid_chimu_v_v;
    std::vector<float> _trk_pida_v_v;

    std::vector<float> _trk_mcs_muon_mom_v;
    std::vector<float> _trk_range_muon_mom_v;
    std::vector<float> _trk_energy_proton_v;
    std::vector<float> _trk_energy_muon_v;
    std::vector<float> _trk_calo_energy_u_v;
    std::vector<float> _trk_calo_energy_v_v;
    std::vector<float> _trk_calo_energy_y_v;

    std::vector<float> _trk_trunk_dEdx_u_v;
    std::vector<float> _trk_trunk_dEdx_v_v;
    std::vector<float> _trk_trunk_dEdx_y_v;

    std::vector<float> _trk_trunk_rr_dEdx_u_v;
    std::vector<float> _trk_trunk_rr_dEdx_v_v;
    std::vector<float> _trk_trunk_rr_dEdx_y_v;

    std::vector<int> _trk_nhits_u_v;
    std::vector<int> _trk_nhits_v_v;
    std::vector<int> _trk_nhits_y_v;

    std::vector<float> _trk_avg_deflection_mean_v;
    std::vector<float> _trk_avg_deflection_stdev_v;
    std::vector<float> _trk_avg_deflection_separation_mean_v;

    std::vector<int> _trk_end_spacepoints_v;
};

TrackAnalysis::TrackAnalysis(const fhicl::ParameterSet &p) : _mcsfitter(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(p.get<fhicl::ParameterSet>("mcsfitmu")))
{
    _CALOproducer = p.get<art::InputTag>("CALOproducer", "pandoracali");
    _PIDproducer = p.get<art::InputTag>("PIDproducer", "pandoracalipid");
    _TRKproducer = p.get<art::InputTag>("TRKproducer", "pandora");
    _BacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
    _Hproducer = p.get<art::InputTag>("Hproducer", "gaushit");
    _CLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
    _EnergyThresholdForHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
    _RecalibrateHits = p.get<bool>("RecalibrateHits", false);
    _ADCtoE = p.get<std::vector<float>>("ADCtoE");
    _EndSpacepointDistance = p.get<float>("EndSpacepointDistance", 5.0);
}

void TrackAnalysis::configure(fhicl::ParameterSet const &p)
{
}

void TrackAnalysis::analyzeEvent(art::Event const &e, bool is_data)
{
}

void TrackAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool is_data, bool selected)
{
    common::ProxyCaloColl_t const &calo_proxy = proxy::getCollection<std::vector<recob::Track>>(e, _TRKproducer,
                                                proxy::withAssociated<anab::Calorimetry>(_CALOproducer));

    common::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, _TRKproducer,
                                              proxy::withAssociated<anab::ParticleID>(_PIDproducer));

    TVector3 nuvtx;
    for (auto pfp : slice_pfp_v)
    {
        if (pfp->IsPrimary())
        {
            double xyz[3] = {};
            auto vtx = pfp.get<recob::Vertex>();
            if (vtx.size() == 1)
            {
                vtx.at(0)->XYZ(xyz);
                nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
            }

            break;
        }
    }

    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> mcp_bkth_assoc;

    if (!is_data)
    {
        art::ValidHandle<std::vector<recob::Hit>> in_hits = e.getValidHandle<std::vector<recob::Hit>>(_Hproducer);
        mcp_bkth_assoc = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(in_hits, e, _BacktrackTag));
    }

    auto sp_handle = e.getValidHandle<std::vector<recob::SpacePoint>>(_CLSproducer);
    std::vector< art::Ptr<recob::SpacePoint> > sp_v;
    for (size_t i_sp = 0; i_sp < sp_handle->size(); i_sp++) {
        sp_v.emplace_back(sp_handle, i_sp);
    }

    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
    {
        auto pfp = slice_pfp_v[i_pfp];
        if (pfp->IsPrimary())
            continue;

        auto trk_v = pfp.get<recob::Track>();

        if (trk_v.size() == 1)
        {
            _trk_score_v.push_back(common::GetTrackShowerScore(pfp));
            auto trk = trk_v.at(0);

            auto trk_prxy_temp = pid_proxy[trk.key()];
            auto pid_prxy_v = trk_prxy_temp.get<anab::ParticleID>();

            auto get_max_pid = [&pid_prxy_v](const int &pdg, const unsigned int plane) 
            {
                return std::max(common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pdg, plane),
                                common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pdg, plane));
            };
            auto get_min_pid = [&pid_prxy_v](const int &pdg, const unsigned int plane) 
            {
                return std::min(common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pdg, plane),
                                common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pdg, plane));
            };
            auto is_fwd_pid_preferred = [&pid_prxy_v](const int &pdg, const unsigned int plane) 
            {
                return common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pdg, plane) >
                    common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pdg, plane);
            };

            float bragg_p = get_max_pid(2212, 2);
            float bragg_mu = get_max_pid(13, 2);
            float bragg_pion = get_max_pid(211, 2);
            float bragg_mip = common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
            float bragg_p_alt_dir = get_min_pid(2212, 2);
            float bragg_mu_alt_dir = get_min_pid(13, 2);
            float bragg_pion_alt_dir = get_min_pid(211, 2);
            bool bragg_p_fwd_preferred = is_fwd_pid_preferred(2212, 2);
            bool bragg_mu_fwd_preferred = is_fwd_pid_preferred(13, 2);
            bool bragg_pion_fwd_preferred = is_fwd_pid_preferred(211, 2);

            float pid_chipr = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 2);
            float pid_chimu = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 2);
            float pid_chipi = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 2);
            float pid_chika = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 2);

            float pida_mean = common::PID(pid_prxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);

            _trk_bragg_p_v.push_back(bragg_p);
            _trk_bragg_mu_v.push_back(bragg_mu);
            _trk_bragg_pion_v.push_back(bragg_pion);
            _trk_bragg_mip_v.push_back(bragg_mip);
            _trk_bragg_p_alt_dir_v.push_back(bragg_p_alt_dir);
            _trk_bragg_mu_alt_dir_v.push_back(bragg_mu_alt_dir);
            _trk_bragg_pion_alt_dir_v.push_back(bragg_pion_alt_dir);
            _trk_bragg_p_fwd_preferred_v.push_back(bragg_p_fwd_preferred);
            _trk_bragg_mu_fwd_preferred_v.push_back(bragg_mu_fwd_preferred);
            _trk_bragg_pion_fwd_preferred_v.push_back(bragg_pion_fwd_preferred);

            _trk_pid_chipr_v.push_back(pid_chipr);
            _trk_pid_chimu_v.push_back(pid_chimu);
            _trk_pid_chipi_v.push_back(pid_chipi);
            _trk_pid_chika_v.push_back(pid_chika);
            _trk_pida_v.push_back(pida_mean);

            //u plane
            float bragg_p_u = get_max_pid(2212, 0);
            float bragg_mu_u = get_max_pid(13, 0);
            float bragg_pion_u = get_max_pid(211, 0);
            float bragg_mip_u = common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 0);
            float bragg_p_alt_dir_u = get_min_pid(2212, 0);
            float bragg_mu_alt_dir_u = get_min_pid(13, 0);
            float bragg_pion_alt_dir_u = get_min_pid(211, 0);
            bool bragg_p_fwd_preferred_u = is_fwd_pid_preferred(2212, 0);
            bool bragg_mu_fwd_preferred_u = is_fwd_pid_preferred(13, 0);
            bool bragg_pion_fwd_preferred_u = is_fwd_pid_preferred(211, 0);

            float pid_chipr_u = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
            float pid_chimu_u = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
            float pid_chipi_u = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
            float pid_chika_u = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);

            float pida_mean_u = common::PID(pid_prxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 0);

            _trk_bragg_p_u_v.push_back(bragg_p_u);
            _trk_bragg_mu_u_v.push_back(bragg_mu_u);
            _trk_bragg_pion_u_v.push_back(bragg_pion_u);
            _trk_bragg_mip_u_v.push_back(bragg_mip_u);
            _trk_bragg_p_alt_dir_u_v.push_back(bragg_p_alt_dir_u);
            _trk_bragg_mu_alt_dir_u_v.push_back(bragg_mu_alt_dir_u);
            _trk_bragg_pion_alt_dir_u_v.push_back(bragg_pion_alt_dir_u);
            _trk_bragg_p_fwd_preferred_u_v.push_back(bragg_p_fwd_preferred_u);
            _trk_bragg_mu_fwd_preferred_u_v.push_back(bragg_mu_fwd_preferred_u);
            _trk_bragg_pion_fwd_preferred_u_v.push_back(bragg_pion_fwd_preferred_u);

            _trk_pid_chipr_u_v.push_back(pid_chipr_u);
            _trk_pid_chimu_u_v.push_back(pid_chimu_u);
            _trk_pid_chipi_u_v.push_back(pid_chipi_u);
            _trk_pid_chika_u_v.push_back(pid_chika_u);
            _trk_pida_u_v.push_back(pida_mean_u);

            //v plane
            float bragg_p_v = get_max_pid(2212, 1);
            float bragg_mu_v = get_max_pid(13, 1);
            float bragg_pion_v = get_max_pid(211, 1);
            float bragg_mip_v = common::PID(pid_prxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 1);
            float bragg_p_alt_dir_v = get_min_pid(2212, 1);
            float bragg_mu_alt_dir_v = get_min_pid(13, 1);
            float bragg_pion_alt_dir_v = get_min_pid(211, 1);
            bool bragg_p_fwd_preferred_v = is_fwd_pid_preferred(2212, 1);
            bool bragg_mu_fwd_preferred_v = is_fwd_pid_preferred(13, 1);
            bool bragg_pion_fwd_preferred_v = is_fwd_pid_preferred(211, 1);

            float pid_chipr_v = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 1);
            float pid_chimu_v = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 1);
            float pid_chipi_v = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 1);
            float pid_chika_v = common::PID(pid_prxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 1);

            float pida_mean_v = common::PID(pid_prxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 1);

            _trk_bragg_p_v_v.push_back(bragg_p_v);
            _trk_bragg_mu_v_v.push_back(bragg_mu_v);
            _trk_bragg_pion_v_v.push_back(bragg_pion_v);
            _trk_bragg_mip_v_v.push_back(bragg_mip_v);
            _trk_bragg_p_alt_dir_v_v.push_back(bragg_p_alt_dir_v);
            _trk_bragg_mu_alt_dir_v_v.push_back(bragg_mu_alt_dir_v);
            _trk_bragg_pion_alt_dir_v_v.push_back(bragg_pion_alt_dir_v);
            _trk_bragg_p_fwd_preferred_v_v.push_back(bragg_p_fwd_preferred_v);
            _trk_bragg_mu_fwd_preferred_v_v.push_back(bragg_mu_fwd_preferred_v);
            _trk_bragg_pion_fwd_preferred_v_v.push_back(bragg_pion_fwd_preferred_v);

            _trk_pid_chipr_v_v.push_back(pid_chipr_v);
            _trk_pid_chimu_v_v.push_back(pid_chimu_v);
            _trk_pid_chipi_v_v.push_back(pid_chipi_v);
            _trk_pid_chika_v_v.push_back(pid_chika_v);
            _trk_pida_v_v.push_back(pida_mean_v);

            float mcs_momentum_muon = _mcsfitter.fitMcs(trk->Trajectory(), 13).bestMomentum();
            float range_momentum_muon = _trkmom.GetTrackMomentum(common::GetSCECorrTrackLength(trk), 13);
            float energy_proton = std::sqrt(std::pow(_trkmom.GetTrackMomentum(common::GetSCECorrTrackLength(trk), 2212), 2) + std::pow(proton->Mass(), 2)) - proton->Mass();
            float energy_muon = std::sqrt(std::pow(mcs_momentum_muon, 2) + std::pow(muon->Mass(), 2)) - muon->Mass();

            _trk_mcs_muon_mom_v.push_back(mcs_momentum_muon);
            _trk_range_muon_mom_v.push_back(range_momentum_muon);
            _trk_energy_proton_v.push_back(energy_proton);
            _trk_energy_muon_v.push_back(energy_muon);
            _trk_calo_energy_u_v.push_back(-1);
            _trk_calo_energy_v_v.push_back(-1);
            _trk_calo_energy_y_v.push_back(-1);

            _trk_dir_x_v.push_back(trk->StartDirection().X());
            _trk_dir_y_v.push_back(trk->StartDirection().Y());
            _trk_dir_z_v.push_back(trk->StartDirection().Z());

            _trk_start_x_v.push_back(trk->Start().X());
            _trk_start_y_v.push_back(trk->Start().Y());
            _trk_start_z_v.push_back(trk->Start().Z());

            float _trk_start_sce[3];
            common::ApplySCECorrectionXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z(), _trk_start_sce);
            _trk_sce_start_x_v.push_back(_trk_start_sce[0]);
            _trk_sce_start_y_v.push_back(_trk_start_sce[1]);
            _trk_sce_start_z_v.push_back(_trk_start_sce[2]);

            _trk_end_x_v.push_back(trk->End().X());
            _trk_end_y_v.push_back(trk->End().Y());
            _trk_end_z_v.push_back(trk->End().Z());

            float _trk_end_sce[3];
            common::ApplySCECorrectionXYZ(trk->End().X(), trk->End().Y(), trk->End().Z(), _trk_end_sce);
            _trk_sce_end_x_v.push_back(_trk_end_sce[0]);
            _trk_sce_end_y_v.push_back(_trk_end_sce[1]);
            _trk_sce_end_z_v.push_back(_trk_end_sce[2]);

            _trk_theta_v.push_back(trk->Theta());
            _trk_phi_v.push_back(trk->Phi());

            _trk_len_v.push_back(common::GetSCECorrTrackLength(trk));

            TVector3 trk_vtx_v;
            trk_vtx_v.SetXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
            trk_vtx_v -= nuvtx;
            _trk_distance_v.push_back(trk_vtx_v.Mag());

            _trk_pfp_id_v.push_back(slice_pfp_v.at(i_pfp)->Self());

            _trk_nhits_u_v.push_back(0);
            _trk_nhits_v_v.push_back(0);
            _trk_nhits_y_v.push_back(0);
            _trk_trunk_dEdx_u_v.push_back(std::numeric_limits<float>::lowest());
            _trk_trunk_dEdx_v_v.push_back(std::numeric_limits<float>::lowest());
            _trk_trunk_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());

            _trk_trunk_rr_dEdx_u_v.push_back(std::numeric_limits<float>::lowest());
            _trk_trunk_rr_dEdx_v_v.push_back(std::numeric_limits<float>::lowest());
            _trk_trunk_rr_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());

            auto calo_v = calo_proxy[trk.key()].get<anab::Calorimetry>();
            for (auto const &calo : calo_v)
            {
                auto const &plane = calo->PlaneID().Plane;
                //auto const &dqdx_values = calo->dQdx();
                auto const &dedx_values = calo->dEdx();
                auto const &rr = calo->ResidualRange();
                auto const &pitch = calo->TrkPitchVec();
                auto const& xyz_v = calo->XYZ();
                std::vector<std::vector<float>> par_values;
                par_values.push_back(rr);
                par_values.push_back(pitch);

                float calo_energy = 0;
                std::vector<float> dqdx_values_corrected, dedx_values_corrected;
                for (size_t i = 0; i < dqdx_values_corrected.size(); i++)
                {
                    float aux_dedx;
                    aux_dedx = common::ModBoxCorrection(dqdx_values_corrected[i]*_ADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
                    dedx_values_corrected.push_back(aux_dedx);
                    calo_energy += aux_dedx * pitch[i];
                }

                float trk_nhits = dedx_values.size();
                float trk_trunk_dEdx = CalculateTrackTrunkdEdxByHits(dedx_values);
                float trk_trunk_rr_dEdx = CalculateTrackTrunkdEdxByRange(dedx_values, rr);

                if (plane == 0)
                {
                    _trk_calo_energy_u_v.back() = calo_energy;
                    _trk_nhits_u_v.back() = trk_nhits;
                    _trk_trunk_dEdx_u_v.back() = trk_trunk_dEdx;
                    _trk_trunk_rr_dEdx_u_v.back() = trk_trunk_rr_dEdx;
                }
                else if (plane == 1)
                {
                    _trk_calo_energy_v_v.back() = calo_energy;
                    _trk_nhits_v_v.back() = trk_nhits;
                    _trk_trunk_dEdx_v_v.back() = trk_trunk_dEdx;
                    _trk_trunk_rr_dEdx_v_v.back() = trk_trunk_rr_dEdx;
                }
                else if (plane == 2)
                {
                    _trk_calo_energy_y_v.back() = calo_energy;
                    _trk_nhits_y_v.back() = trk_nhits;
                    _trk_trunk_dEdx_y_v.back() = trk_trunk_dEdx;
                    _trk_trunk_rr_dEdx_y_v.back() = trk_trunk_rr_dEdx;
                }
            }

            CalculateTrackDeflections(trk, _trk_avg_deflection_mean_v, _trk_avg_deflection_stdev_v, _trk_avg_deflection_separation_mean_v);

            int nPoints = 0;
            float distSquared = _EndSpacepointDistance*_EndSpacepointDistance;
            TVector3 trkEnd(_trk_end_sce[0], _trk_end_sce[1], _trk_end_sce[2]);
            for (auto &sp : sp_v) {
                float _sp_sce[3];
                common::ApplySCECorrectionXYZ(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2], _sp_sce);
                TVector3 spacePoint(_sp_sce[0], _sp_sce[1], _sp_sce[2]);
                if ((trkEnd - spacePoint).Mag2() < distSquared) 
                    nPoints++;
            }

            _trk_end_spacepoints_v.push_back(nPoints);
        }
        else
        {
            fillDefault();
        }
    } 
}

void TrackAnalysis::fillDefault()
{
    _trk_pfp_id_v.push_back(std::numeric_limits<int>::lowest());
    _trk_score_v.push_back(std::numeric_limits<float>::lowest());

    _trk_distance_v.push_back(std::numeric_limits<float>::lowest());

    _trk_theta_v.push_back(std::numeric_limits<float>::lowest());
    _trk_phi_v.push_back(std::numeric_limits<float>::lowest());

    _trk_dir_x_v.push_back(std::numeric_limits<float>::lowest());
    _trk_dir_y_v.push_back(std::numeric_limits<float>::lowest());
    _trk_dir_z_v.push_back(std::numeric_limits<float>::lowest());

    _trk_start_x_v.push_back(std::numeric_limits<float>::lowest());
    _trk_start_y_v.push_back(std::numeric_limits<float>::lowest());
    _trk_start_z_v.push_back(std::numeric_limits<float>::lowest());

    _trk_sce_start_x_v.push_back(std::numeric_limits<float>::lowest());
    _trk_sce_start_y_v.push_back(std::numeric_limits<float>::lowest());
    _trk_sce_start_z_v.push_back(std::numeric_limits<float>::lowest());

    _trk_end_x_v.push_back(std::numeric_limits<float>::lowest());
    _trk_end_y_v.push_back(std::numeric_limits<float>::lowest());
    _trk_end_z_v.push_back(std::numeric_limits<float>::lowest());

    _trk_sce_end_x_v.push_back(std::numeric_limits<float>::lowest());
    _trk_sce_end_y_v.push_back(std::numeric_limits<float>::lowest());
    _trk_sce_end_z_v.push_back(std::numeric_limits<float>::lowest());

    _trk_len_v.push_back(std::numeric_limits<float>::lowest());

    _trk_bragg_p_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mip_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_p_alt_dir_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_alt_dir_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_alt_dir_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_p_fwd_preferred_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_fwd_preferred_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_fwd_preferred_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chipr_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chika_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chipi_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chimu_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pida_v.push_back(std::numeric_limits<float>::lowest());

    _trk_bragg_p_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mip_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_p_alt_dir_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_alt_dir_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_alt_dir_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_p_fwd_preferred_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_fwd_preferred_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_fwd_preferred_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chipr_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chika_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chipi_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chimu_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pida_u_v.push_back(std::numeric_limits<float>::lowest());

    _trk_bragg_p_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mip_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_p_alt_dir_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_alt_dir_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_alt_dir_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_p_fwd_preferred_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_mu_fwd_preferred_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_bragg_pion_fwd_preferred_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chipr_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chika_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chipi_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pid_chimu_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_pida_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_mcs_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
    _trk_range_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
    _trk_energy_proton_v.push_back(std::numeric_limits<float>::lowest());
    _trk_energy_muon_v.push_back(std::numeric_limits<float>::lowest());
    _trk_calo_energy_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_calo_energy_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_calo_energy_y_v.push_back(std::numeric_limits<float>::lowest());

    _trk_trunk_dEdx_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_trunk_dEdx_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_trunk_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());

    _trk_trunk_rr_dEdx_u_v.push_back(std::numeric_limits<float>::lowest());
    _trk_trunk_rr_dEdx_v_v.push_back(std::numeric_limits<float>::lowest());
    _trk_trunk_rr_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());

    _trk_nhits_u_v.push_back(std::numeric_limits<int>::lowest());
    _trk_nhits_v_v.push_back(std::numeric_limits<int>::lowest());
    _trk_nhits_y_v.push_back(std::numeric_limits<int>::lowest());

    _trk_avg_deflection_mean_v.push_back(std::numeric_limits<float>::lowest());
    _trk_avg_deflection_stdev_v.push_back(std::numeric_limits<float>::lowest());
    _trk_avg_deflection_separation_mean_v.push_back(std::numeric_limits<float>::lowest());

    _trk_end_spacepoints_v.push_back(std::numeric_limits<int>::lowest());
}

void TrackAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("trk_bragg_p_v", "std::vector<float>", &_trk_bragg_p_v);
    _tree->Branch("trk_bragg_mu_v", "std::vector<float>", &_trk_bragg_mu_v);
    _tree->Branch("trk_bragg_pion_v", "std::vector<float>", &_trk_bragg_pion_v);
    _tree->Branch("trk_bragg_mip_v", "std::vector<float>", &_trk_bragg_mip_v);
    _tree->Branch("trk_bragg_p_alt_dir_v", "std::vector<float>", &_trk_bragg_p_alt_dir_v);
    _tree->Branch("trk_bragg_mu_alt_dir_v", "std::vector<float>", &_trk_bragg_mu_alt_dir_v);
    _tree->Branch("trk_bragg_pion_alt_dir_v", "std::vector<float>", &_trk_bragg_pion_alt_dir_v);
    _tree->Branch("trk_bragg_p_fwd_preferred_v", "std::vector< bool >", &_trk_bragg_p_fwd_preferred_v);
    _tree->Branch("trk_bragg_mu_fwd_preferred_v", "std::vector< bool >", &_trk_bragg_mu_fwd_preferred_v);
    _tree->Branch("trk_bragg_pion_fwd_preferred_v", "std::vector< bool >", &_trk_bragg_pion_fwd_preferred_v);
    _tree->Branch("trk_pida_v", "std::vector<float>", &_trk_pida_v);
    _tree->Branch("trk_pid_chipr_v", "std::vector<float>", &_trk_pid_chipr_v);
    _tree->Branch("trk_pid_chipi_v", "std::vector<float>", &_trk_pid_chipi_v);
    _tree->Branch("trk_pid_chika_v", "std::vector<float>", &_trk_pid_chika_v);
    _tree->Branch("trk_pid_chimu_v", "std::vector<float>", &_trk_pid_chimu_v);

    _tree->Branch("trk_bragg_p_u_v", "std::vector<float>", &_trk_bragg_p_u_v);
    _tree->Branch("trk_bragg_mu_u_v", "std::vector<float>", &_trk_bragg_mu_u_v);
    _tree->Branch("trk_bragg_pion_u_v", "std::vector<float>", &_trk_bragg_pion_u_v);
    _tree->Branch("trk_bragg_mip_u_v", "std::vector<float>", &_trk_bragg_mip_u_v);
    _tree->Branch("trk_bragg_p_alt_dir_u_v", "std::vector<float>", &_trk_bragg_p_alt_dir_u_v);
    _tree->Branch("trk_bragg_mu_alt_dir_u_v", "std::vector<float>", &_trk_bragg_mu_alt_dir_u_v);
    _tree->Branch("trk_bragg_pion_alt_dir_u_v", "std::vector<float>", &_trk_bragg_pion_alt_dir_u_v);
    _tree->Branch("trk_bragg_p_fwd_preferred_u_v", "std::vector< bool >", &_trk_bragg_p_fwd_preferred_u_v);
    _tree->Branch("trk_bragg_mu_fwd_preferred_u_v", "std::vector< bool >", &_trk_bragg_mu_fwd_preferred_u_v);
    _tree->Branch("trk_bragg_pion_fwd_preferred_u_v", "std::vector< bool >", &_trk_bragg_pion_fwd_preferred_u_v);
    _tree->Branch("trk_pida_u_v", "std::vector<float>", &_trk_pida_u_v);
    _tree->Branch("trk_pid_chipr_u_v", "std::vector<float>", &_trk_pid_chipr_u_v);
    _tree->Branch("trk_pid_chipi_u_v", "std::vector<float>", &_trk_pid_chipi_u_v);
    _tree->Branch("trk_pid_chika_u_v", "std::vector<float>", &_trk_pid_chika_u_v);
    _tree->Branch("trk_pid_chimu_u_v", "std::vector<float>", &_trk_pid_chimu_u_v);

    _tree->Branch("trk_bragg_p_v_v", "std::vector<float>", &_trk_bragg_p_v_v);
    _tree->Branch("trk_bragg_mu_v_v", "std::vector<float>", &_trk_bragg_mu_v_v);
    _tree->Branch("trk_bragg_pion_v_v", "std::vector<float>", &_trk_bragg_pion_v_v);
    _tree->Branch("trk_bragg_mip_v_v", "std::vector<float>", &_trk_bragg_mip_v_v);
    _tree->Branch("trk_bragg_p_alt_dir_v_v", "std::vector<float>", &_trk_bragg_p_alt_dir_v_v);
    _tree->Branch("trk_bragg_mu_alt_dir_v_v", "std::vector<float>", &_trk_bragg_mu_alt_dir_v_v);
    _tree->Branch("trk_bragg_pion_alt_dir_v_v", "std::vector<float>", &_trk_bragg_pion_alt_dir_v_v);
    _tree->Branch("trk_bragg_p_fwd_preferred_v_v", "std::vector< bool >", &_trk_bragg_p_fwd_preferred_v_v);
    _tree->Branch("trk_bragg_mu_fwd_preferred_v_v", "std::vector< bool >", &_trk_bragg_mu_fwd_preferred_v_v);
    _tree->Branch("trk_bragg_pion_fwd_preferred_v_v", "std::vector< bool >", &_trk_bragg_pion_fwd_preferred_v_v);
    _tree->Branch("trk_pida_v_v", "std::vector<float>", &_trk_pida_v_v);
    _tree->Branch("trk_pid_chipr_v_v", "std::vector<float>", &_trk_pid_chipr_v_v);
    _tree->Branch("trk_pid_chipi_v_v", "std::vector<float>", &_trk_pid_chipi_v_v);
    _tree->Branch("trk_pid_chika_v_v", "std::vector<float>", &_trk_pid_chika_v_v);
    _tree->Branch("trk_pid_chimu_v_v", "std::vector<float>", &_trk_pid_chimu_v_v);

    _tree->Branch("trk_pfp_id_v", "std::vector<size_t>", &_trk_pfp_id_v);
    _tree->Branch("trk_score_v", "std::vector<float>", &_trk_score_v);
    _tree->Branch("trk_dir_x_v", "std::vector<float>", &_trk_dir_x_v);
    _tree->Branch("trk_dir_y_v", "std::vector<float>", &_trk_dir_y_v);
    _tree->Branch("trk_dir_z_v", "std::vector<float>", &_trk_dir_z_v);

    _tree->Branch("trk_start_x_v", "std::vector<float>", &_trk_start_x_v);
    _tree->Branch("trk_start_y_v", "std::vector<float>", &_trk_start_y_v);
    _tree->Branch("trk_start_z_v", "std::vector<float>", &_trk_start_z_v);

    _tree->Branch("trk_sce_start_x_v", "std::vector<float>", &_trk_sce_start_x_v);
    _tree->Branch("trk_sce_start_y_v", "std::vector<float>", &_trk_sce_start_y_v);
    _tree->Branch("trk_sce_start_z_v", "std::vector<float>", &_trk_sce_start_z_v);

    _tree->Branch("trk_end_x_v", "std::vector<float>", &_trk_end_x_v);
    _tree->Branch("trk_end_y_v", "std::vector<float>", &_trk_end_y_v);
    _tree->Branch("trk_end_z_v", "std::vector<float>", &_trk_end_z_v);

    _tree->Branch("trk_sce_end_x_v", "std::vector<float>", &_trk_sce_end_x_v);
    _tree->Branch("trk_sce_end_y_v", "std::vector<float>", &_trk_sce_end_y_v);
    _tree->Branch("trk_sce_end_z_v", "std::vector<float>", &_trk_sce_end_z_v);

    _tree->Branch("trk_distance_v", "std::vector<float>", &_trk_distance_v);
    _tree->Branch("trk_theta_v", "std::vector<float>", &_trk_theta_v);
    _tree->Branch("trk_phi_v", "std::vector<float>", &_trk_phi_v);

    _tree->Branch("trk_len_v", "std::vector<float>", &_trk_len_v);
    _tree->Branch("trk_mcs_muon_mom_v", "std::vector<float>", &_trk_mcs_muon_mom_v);
    _tree->Branch("trk_range_muon_mom_v", "std::vector<float>", &_trk_range_muon_mom_v);
    _tree->Branch("trk_energy_proton_v", "std::vector<float>", &_trk_energy_proton_v);
    _tree->Branch("trk_energy_muon_v", "std::vector<float>", &_trk_energy_muon_v);
    _tree->Branch("trk_calo_energy_u_v", "std::vector<float>", &_trk_calo_energy_u_v);
    _tree->Branch("trk_calo_energy_v_v", "std::vector<float>", &_trk_calo_energy_v_v);
    _tree->Branch("trk_calo_energy_y_v", "std::vector<float>", &_trk_calo_energy_y_v);

    _tree->Branch("trk_trunk_dEdx_u_v", "std::vector<float>", &_trk_trunk_dEdx_u_v);
    _tree->Branch("trk_trunk_dEdx_v_v", "std::vector<float>", &_trk_trunk_dEdx_v_v);
    _tree->Branch("trk_trunk_dEdx_y_v", "std::vector<float>", &_trk_trunk_dEdx_y_v);

    _tree->Branch("trk_trunk_rr_dEdx_u_v", "std::vector<float>", &_trk_trunk_rr_dEdx_u_v);
    _tree->Branch("trk_trunk_rr_dEdx_v_v", "std::vector<float>", &_trk_trunk_rr_dEdx_v_v);
    _tree->Branch("trk_trunk_rr_dEdx_y_v", "std::vector<float>", &_trk_trunk_rr_dEdx_y_v);

    _tree->Branch("trk_nhits_u_v", "std::vector<int>", &_trk_nhits_u_v);
    _tree->Branch("trk_nhits_v_v", "std::vector<int>", &_trk_nhits_v_v);
    _tree->Branch("trk_nhits_y_v", "std::vector<int>", &_trk_nhits_y_v);

    _tree->Branch("trk_avg_deflection_mean_v", "std::vector<float>", &_trk_avg_deflection_mean_v);
    _tree->Branch("trk_avg_deflection_stdev_v", "std::vector<float>", &_trk_avg_deflection_stdev_v);
    _tree->Branch("trk_avg_deflection_separation_mean_v", "std::vector<float>", &_trk_avg_deflection_separation_mean_v);

    _tree->Branch("trk_end_spacepoints_v", "std::vector<int>", &_trk_end_spacepoints_v);
}

void TrackAnalysis::resetTTree(TTree *_tree)
{
    _trk_bragg_p_v.clear();
    _trk_bragg_mu_v.clear();
    _trk_bragg_pion_v.clear();
    _trk_bragg_mip_v.clear();
    _trk_bragg_p_alt_dir_v.clear();
    _trk_bragg_mu_alt_dir_v.clear();
    _trk_bragg_pion_alt_dir_v.clear();
    _trk_bragg_p_fwd_preferred_v.clear();
    _trk_bragg_mu_fwd_preferred_v.clear();
    _trk_bragg_pion_fwd_preferred_v.clear();
    _trk_pida_v.clear();
    _trk_pid_chipr_v.clear();
    _trk_pid_chika_v.clear();
    _trk_pid_chipi_v.clear();
    _trk_pid_chimu_v.clear();

    _trk_bragg_p_u_v.clear();
    _trk_bragg_mu_u_v.clear();
    _trk_bragg_pion_u_v.clear();
    _trk_bragg_mip_u_v.clear();
    _trk_bragg_p_alt_dir_u_v.clear();
    _trk_bragg_mu_alt_dir_u_v.clear();
    _trk_bragg_pion_alt_dir_u_v.clear();
    _trk_bragg_p_fwd_preferred_u_v.clear();
    _trk_bragg_mu_fwd_preferred_u_v.clear();
    _trk_bragg_pion_fwd_preferred_u_v.clear();
    _trk_pida_u_v.clear();
    _trk_pid_chipr_u_v.clear();
    _trk_pid_chika_u_v.clear();
    _trk_pid_chipi_u_v.clear();
    _trk_pid_chimu_u_v.clear();

    _trk_bragg_p_v_v.clear();
    _trk_bragg_mu_v_v.clear();
    _trk_bragg_pion_v_v.clear();
    _trk_bragg_mip_v_v.clear();
    _trk_bragg_p_alt_dir_v_v.clear();
    _trk_bragg_mu_alt_dir_v_v.clear();
    _trk_bragg_pion_alt_dir_v_v.clear();
    _trk_bragg_p_fwd_preferred_v_v.clear();
    _trk_bragg_mu_fwd_preferred_v_v.clear();
    _trk_bragg_pion_fwd_preferred_v_v.clear();
    _trk_pida_v_v.clear();
    _trk_pid_chipr_v_v.clear();
    _trk_pid_chika_v_v.clear();
    _trk_pid_chipi_v_v.clear();
    _trk_pid_chimu_v_v.clear();

    _trk_pfp_id_v.clear();

    _trk_start_x_v.clear();
    _trk_start_y_v.clear();
    _trk_start_z_v.clear();

    _trk_sce_start_x_v.clear();
    _trk_sce_start_y_v.clear();
    _trk_sce_start_z_v.clear();

    _trk_end_x_v.clear();
    _trk_end_y_v.clear();
    _trk_end_z_v.clear();

    _trk_sce_end_x_v.clear();
    _trk_sce_end_y_v.clear();
    _trk_sce_end_z_v.clear();

    _trk_dir_x_v.clear();
    _trk_dir_y_v.clear();
    _trk_dir_z_v.clear();
    _trk_distance_v.clear();

    _trk_theta_v.clear();
    _trk_phi_v.clear();

    _trk_len_v.clear();

    _trk_mcs_muon_mom_v.clear();
    _trk_range_muon_mom_v.clear();
    _trk_energy_muon_v.clear();
    _trk_energy_proton_v.clear();
    _trk_calo_energy_u_v.clear();
    _trk_calo_energy_v_v.clear();
    _trk_calo_energy_y_v.clear();

    _trk_trunk_dEdx_u_v.clear();
    _trk_trunk_dEdx_v_v.clear();
    _trk_trunk_dEdx_y_v.clear();

    _trk_trunk_rr_dEdx_u_v.clear();
    _trk_trunk_rr_dEdx_v_v.clear();
    _trk_trunk_rr_dEdx_y_v.clear();

    _trk_nhits_u_v.clear();
    _trk_nhits_v_v.clear();
    _trk_nhits_y_v.clear();

    _trk_avg_deflection_mean_v.clear();
    _trk_avg_deflection_stdev_v.clear();
    _trk_avg_deflection_separation_mean_v.clear();

    _trk_end_spacepoints_v.clear();
}

float TrackAnalysis::CalculateTrackTrunkdEdxByHits(const std::vector<float> &dEdx_values) 
{
  unsigned int trk_nhits = dEdx_values.size();
  int first_hit_id = trk_nhits - 3 - 1;
  int last_hit_id = trk_nhits - (int)(trk_nhits/3) - 1; 

    if (first_hit_id - last_hit_id < 5) 
    {
        return std::numeric_limits<float>::lowest();
    }
    else 
    {
        std::vector<float> trk_trunk_dEdx_values;
        trk_trunk_dEdx_values.reserve(first_hit_id - last_hit_id); 

        for (int i = trk_nhits - 1; i >= 0; i--) {
            if (i > first_hit_id) 
                continue;

            trk_trunk_dEdx_values.push_back(dEdx_values[i]);

            if (i < last_hit_id) 
                break;  
        }

        float median;
        std::sort(trk_trunk_dEdx_values.begin(), trk_trunk_dEdx_values.end());
        if (trk_trunk_dEdx_values.size() % 2 == 0) 
            median = 0.5 * (trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2 - 1] + trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2]);
        else 
            median = trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2];

        double sum = std::accumulate(std::begin(trk_trunk_dEdx_values), std::end(trk_trunk_dEdx_values), 0.0);
        double m =  sum / trk_trunk_dEdx_values.size();
    
        double accum = 0.0;
        std::for_each(std::begin(trk_trunk_dEdx_values), std::end(trk_trunk_dEdx_values), [&](const double d) {accum += (d - m) * (d - m);});
        double stdev = sqrt(accum / (trk_trunk_dEdx_values.size()-1));

        std::vector<float> trk_trunk_dEdx_values_trimmed;
        trk_trunk_dEdx_values_trimmed.reserve(first_hit_id - last_hit_id);
        for (unsigned int i = 0; i < trk_trunk_dEdx_values.size(); i++) {
            if (trk_trunk_dEdx_values[i] <= median + stdev) 
                trk_trunk_dEdx_values_trimmed.push_back(trk_trunk_dEdx_values[i]);
        }

        double sum_trimmed = std::accumulate(std::begin(trk_trunk_dEdx_values_trimmed), std::end(trk_trunk_dEdx_values_trimmed), 0.0);
        double trk_dEdx_trunk =  sum_trimmed / trk_trunk_dEdx_values_trimmed.size();

        return trk_dEdx_trunk;
    }
}

float TrackAnalysis::CalculateTrackTrunkdEdxByRange(const std::vector<float> &dEdx_per_hit, const std::vector<float> &rr_per_hit)
{
    const auto nhits_skip = 3u;
    const auto l_frac = 1.f/3;

    if (rr_per_hit.size() <= nhits_skip)
        return -std::numeric_limits<float>::max();

    std::vector<std::pair<float, unsigned int> > rr_i;
    float max_rr = -std::numeric_limits<float>::max();
    for (unsigned int i = 0; i < rr_per_hit.size(); ++i)
    {
        const auto rr = rr_per_hit.at(i);
        max_rr = std::max(max_rr, rr);
        rr_i.emplace_back(rr, i);
    }

    const auto rr_cutoff = max_rr * l_frac;

    std::sort(rr_i.begin(), rr_i.end(), [](auto &a, auto &b) {
        return a.first > b.first;
    });

    std::vector<float> dEdx_per_hit_start;
    for (unsigned int i = nhits_skip; i < rr_i.size(); ++i)
    {
        const auto entry = rr_i.at(i);
        const auto rr = entry.first;
        const auto hit_id = entry.second;

        if (rr < rr_cutoff)
            continue;

        dEdx_per_hit_start.push_back(dEdx_per_hit.at(hit_id));
    }

    const auto n_hits = dEdx_per_hit_start.size();
    if (n_hits == 0)
        return -std::numeric_limits<float>::max();

    std::sort(dEdx_per_hit_start.begin(), dEdx_per_hit_start.end());
    const auto median = dEdx_per_hit_start.at(n_hits / 2); 

    float total = 0.f;
    for (const auto &dEdx : dEdx_per_hit_start)
        total += dEdx;
    const auto mean = total / static_cast<float>(n_hits);

    float sqr_sum = 0.f;
    for (const auto &dEdx : dEdx_per_hit_start)
        sqr_sum += std::pow(dEdx - mean, 2);
    const auto variance = sqr_sum / static_cast<float>(n_hits);

    float trun_tot = 0.f;
    unsigned int trun_nhits = 0;
    for (const auto &dEdx : dEdx_per_hit_start)
    {
        if (std::pow(dEdx - median, 2) > variance)
            continue;

        trun_tot += dEdx;
        trun_nhits++;
    }

    if (trun_nhits == 0)
        return -std::numeric_limits<float>::max();

    return trun_tot / static_cast<float>(trun_nhits);
}

void TrackAnalysis::CalculateTrackDeflections(const art::Ptr<recob::Track> &trk, std::vector<float> &mean_v, std::vector<float> &stdev_v, std::vector<float> &separation_mean_v)
{
    std::vector<size_t> valid_points;
    auto first_valid_point = trk->FirstValidPoint();
    valid_points.push_back(first_valid_point);
    auto next_valid_point = trk->NextValidPoint(first_valid_point + 1);
    while (next_valid_point != recob::TrackTrajectory::InvalidIndex)
    {
        valid_points.push_back(next_valid_point);
        next_valid_point = trk->NextValidPoint(next_valid_point + 1);
    }

    if (valid_points.size() < 3) 
    {
        mean_v.push_back(0);
        stdev_v.push_back(0);
        separation_mean_v.push_back(0);
    }
    else 
    {
        std::vector<float> theta_v;
        float theta_sum = 0.f;
        float sep_sum = 0.f;
        for (unsigned int i = 1; i < valid_points.size(); ++i) {
            const auto dir = trk->DirectionAtPoint(valid_points.at(i));
            const auto dir_prev = trk->DirectionAtPoint(valid_points.at(i - 1));

            const auto cos_theta = std::min(1.f, std::max(-1.f, static_cast<float>(dir.Dot(dir_prev))));
            const auto theta = std::acos(cos_theta);

            theta_sum += theta;
            theta_v.push_back(theta);

            const TVector3 point(trk->LocationAtPoint(valid_points.at(i)).X(), trk->LocationAtPoint(valid_points.at(i)).Y(), trk->LocationAtPoint(valid_points.at(i)).Z());
            const TVector3 point_prev(trk->LocationAtPoint(valid_points.at(i - 1)).X(), trk->LocationAtPoint(valid_points.at(i - 1)).Y(), trk->LocationAtPoint(valid_points.at(i - 1)).Z());
            const TVector3 separation = point - point_prev; 
            sep_sum += separation.Mag();
        }

        float theta_mean = theta_sum / static_cast<float>(theta_v.size());
        float sep_mean = sep_sum / static_cast<float>(theta_v.size());

        float theta_dif_sum = 0.f;
        for (const auto &theta : theta_v) 
        {
            theta_dif_sum += std::pow(theta - theta_mean, 2);
        }

        const auto variance = theta_dif_sum / static_cast<float>(theta_v.size() - 1);

        mean_v.push_back(theta_mean);
        stdev_v.push_back(std::sqrt(variance));
        separation_mean_v.push_back(sep_mean);
    }
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)
} 
#endif