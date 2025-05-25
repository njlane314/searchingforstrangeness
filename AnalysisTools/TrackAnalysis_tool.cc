#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonDefs/Types.h"

#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/PIDFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Geometry.h"

#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_proton_muon_lookup.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/CalibrationFuncs.h"

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

        void analyseEvent(art::Event const &e, bool fData) override;

        void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

        void SaveTruth(art::Event const &e);

        void fillDefault();

        void setBranches(TTree *_tree) override;

        void resetTTree(TTree *_tree) override;

    private:
        float CalculateTrackTrunkdEdxByHits(const std::vector<float> &dedxPerHit);
        float CalculateTrackTrunkdEdxByRange(const std::vector<float> &dedxPerHit, const std::vector<float> &residualRangePerHit);
        void CalculateTrackDeflections(const art::Ptr<recob::Track> &trk, std::vector<float> &mean_v, std::vector<float> &stdev_v, std::vector<float> &separation_mean_v);

        const trkf::TrackMomentumCalculator _trkmom;
        const trkf::TrajectoryMCSFitter _mcsfitter;

        TParticlePDG *proton = TDatabasePDG::Instance()->GetParticle(2212);
        TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);

        common::LLRPID llr_pid_calculator;
        common::ProtonMuonLookUpParameters protonmuon_parameters;
        common::CorrectionLookUpParameters correction_parameters;

        art::InputTag fCALOproducer;
        art::InputTag fPIDproducer;
        art::InputTag fTRKproducer;
        art::InputTag fBacktrackTag;
        art::InputTag fHproducer;
        art::InputTag fCLSproducer;

        bool fRecalibrateHits;
        float fEnergyThresholdForMCHits;
        std::vector<float> fADCtoE;
        float fEndSpacepointDistance;

        int _run, _sub, _evt;

        std::vector<size_t> _trk_pfp_id_v;

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

        std::vector<float> _trk_bragg_p_v;
        std::vector<float> _trk_bragg_mu_v;
        std::vector<float> _trk_bragg_pion_v;
        std::vector<float> _trk_bragg_mip_v;
        std::vector<float> _trk_bragg_p_alt_dir_v;
        std::vector<float> _trk_bragg_mu_alt_dir_v;
        std::vector<float> _trk_bragg_pion_alt_dir_v;
        std::vector<bool> _trk_bragg_p_fwd_preferred_v;
        std::vector<bool> _trk_bragg_mu_fwd_preferred_v;
        std::vector<bool> _trk_bragg_pion_fwd_preferred_v;
        std::vector<float> _trk_pid_chipr_v;
        std::vector<float> _trk_pid_chika_v;
        std::vector<float> _trk_pid_chipi_v;
        std::vector<float> _trk_pid_chimu_v;
        std::vector<float> _trk_pida_v;

        std::vector<float> _trk_bragg_p_u_v;
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

        std::vector<float> _trk_bragg_p_v_v;
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

        std::vector<float> _trk_llr_pid_u_v;
        std::vector<float> _trk_llr_pid_v_v;
        std::vector<float> _trk_llr_pid_y_v;
        std::vector<float> _trk_llr_pid_v;

        std::vector<float> _trk_llr_pid_score_v;

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
        fCALOproducer = p.get<art::InputTag>("CALOproducer");
        fPIDproducer = p.get<art::InputTag>("PIDproducer");
        fTRKproducer = p.get<art::InputTag>("TRKproducer");
        fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
        fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
        fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
        fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
        fRecalibrateHits = p.get<bool>("RecalibrateHits", false);
        fADCtoE = p.get<std::vector<float>>("ADCtoE");
        fEndSpacepointDistance = p.get<float>("EndSpacepointDistance", 5.0);

        llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
        llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
        llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);

        llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
        llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
        llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);

        llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
        llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
        llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);

        if (fRecalibrateHits)
        {
            llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
            llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);

            llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
            llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);

            llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
            llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
        }
    }

    void TrackAnalysis::configure(fhicl::ParameterSet const &p)
    {
    }

    void TrackAnalysis::analyseEvent(art::Event const &e, bool fData)
    {
        _evt = e.event();
        _sub = e.subRun();
        _run = e.run();
        std::cout << "[TrackAnalysis::analyseEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;
    }

    void TrackAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
    {
        std::cout << "track slice" << std::endl;
        common::ProxyCaloColl_t const &calo_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                            proxy::withAssociated<anab::Calorimetry>(fCALOproducer));

        common::ProxyPIDColl_t const &pid_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer,
                                                                                                            proxy::withAssociated<anab::ParticleID>(fPIDproducer));

        TVector3 nuvtx;
        for (auto pfp : slice_pfp_v)
        {
            if (pfp->IsPrimary())
            {
                double xyz[3] = {};
                auto vtx = pfp.get<recob::Vertex>();
                if (vtx.size() != 1)
                {
                    std::cout << "ERROR. Found neutrino PFP w/ != 1 associated vertices..." << std::endl;
                }
                else
                {
                    vtx.at(0)->XYZ(xyz);
                    nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
                }
                break;
            }
        }

        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

        if (!fData)
        {
            art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
            assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
        }

        auto spacePointHandle = e.getValidHandle<std::vector<recob::SpacePoint>>(fCLSproducer);
        std::vector< art::Ptr<recob::SpacePoint> > spacePointCollection;
        for (size_t i_sp = 0; i_sp < spacePointHandle->size(); i_sp++) {
            spacePointCollection.emplace_back(spacePointHandle, i_sp);
        }

        for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++)
        {
            auto pfp = slice_pfp_v[i_pfp];
            if (pfp->IsPrimary())
                continue;

            auto trk_v = pfp.get<recob::Track>();

            if (trk_v.size() == 1)
            {
                auto trk = trk_v.at(0);

                auto trkpxy2 = pid_proxy[trk.key()];
                auto pidpxy_v = trkpxy2.get<anab::ParticleID>();

                auto GetMaxPID = [&pidpxy_v](const int &pdg, const unsigned int plane) 
                {
                    return std::max(common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pdg, plane),
                                    common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pdg, plane));
                };
                auto GetMinPID = [&pidpxy_v](const int &pdg, const unsigned int plane) 
                {
                    return std::min(common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pdg, plane),
                                    common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pdg, plane));
                };
                auto IsFwdPIDPreferred = [&pidpxy_v](const int &pdg, const unsigned int plane) 
                {
                    return common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, pdg, plane)>
                        common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kBackward, pdg, plane);
                };

                float bragg_p = GetMaxPID(2212, 2);
                float bragg_mu = GetMaxPID(13, 2);
                float bragg_pion = GetMaxPID(211, 2);
                float bragg_mip = common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2);
                float bragg_p_alt_dir = GetMinPID(2212, 2);
                float bragg_mu_alt_dir = GetMinPID(13, 2);
                float bragg_pion_alt_dir = GetMinPID(211, 2);
                bool bragg_p_fwd_preferred = IsFwdPIDPreferred(2212, 2);
                bool bragg_mu_fwd_preferred = IsFwdPIDPreferred(13, 2);
                bool bragg_pion_fwd_preferred = IsFwdPIDPreferred(211, 2);

                float pidchipr = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 2);
                float pidchimu = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 2);
                float pidchipi = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 2);
                float pidchika = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 2);

                float pida_mean = common::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2);

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

                _trk_pid_chipr_v.push_back(pidchipr);
                _trk_pid_chimu_v.push_back(pidchimu);
                _trk_pid_chipi_v.push_back(pidchipi);
                _trk_pid_chika_v.push_back(pidchika);
                _trk_pida_v.push_back(pida_mean);

                float bragg_p_u = GetMaxPID(2212, 0);
                float bragg_mu_u = GetMaxPID(13, 0);
                float bragg_pion_u = GetMaxPID(211, 0);
                float bragg_mip_u = common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 0);
                float bragg_p_alt_dir_u = GetMinPID(2212, 0);
                float bragg_mu_alt_dir_u = GetMinPID(13, 0);
                float bragg_pion_alt_dir_u = GetMinPID(211, 0);
                bool bragg_p_fwd_preferred_u = IsFwdPIDPreferred(2212, 0);
                bool bragg_mu_fwd_preferred_u = IsFwdPIDPreferred(13, 0);
                bool bragg_pion_fwd_preferred_u = IsFwdPIDPreferred(211, 0);

                float pidchipr_u = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 0);
                float pidchimu_u = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 0);
                float pidchipi_u = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 0);
                float pidchika_u = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 0);

                float pida_mean_u = common::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 0);

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

                _trk_pid_chipr_u_v.push_back(pidchipr_u);
                _trk_pid_chimu_u_v.push_back(pidchimu_u);
                _trk_pid_chipi_u_v.push_back(pidchipi_u);
                _trk_pid_chika_u_v.push_back(pidchika_u);
                _trk_pida_u_v.push_back(pida_mean_u);

                float bragg_p_v = GetMaxPID(2212, 1);
                float bragg_mu_v = GetMaxPID(13, 1);
                float bragg_pion_v = GetMaxPID(211, 1);
                float bragg_mip_v = common::PID(pidpxy_v[0], "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 1);
                float bragg_p_alt_dir_v = GetMinPID(2212, 1);
                float bragg_mu_alt_dir_v = GetMinPID(13, 1);
                float bragg_pion_alt_dir_v = GetMinPID(211, 1);
                bool bragg_p_fwd_preferred_v = IsFwdPIDPreferred(2212, 1);
                bool bragg_mu_fwd_preferred_v = IsFwdPIDPreferred(13, 1);
                bool bragg_pion_fwd_preferred_v = IsFwdPIDPreferred(211, 1);

                float pidchipr_v = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 2212, 1);
                float pidchimu_v = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 13, 1);
                float pidchipi_v = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 211, 1);
                float pidchika_v = common::PID(pidpxy_v[0], "Chi2", anab::kGOF, anab::kForward, 321, 1);

                float pida_mean_v = common::PID(pidpxy_v[0], "PIDA_mean", anab::kPIDA, anab::kForward, 0, 1);

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

                _trk_pid_chipr_v_v.push_back(pidchipr_v);
                _trk_pid_chimu_v_v.push_back(pidchimu_v);
                _trk_pid_chipi_v_v.push_back(pidchipi_v);
                _trk_pid_chika_v_v.push_back(pidchika_v);
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

                _trk_pfp_id_v.push_back(i_pfp);

                _trk_llr_pid_u_v.push_back(0);
                _trk_llr_pid_v_v.push_back(0);
                _trk_llr_pid_y_v.push_back(0);
                _trk_llr_pid_v.push_back(0);
                _trk_llr_pid_score_v.push_back(0);

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
                    if (plane < 0 || plane > 2) {
                        std::cerr << "Skipping calorimetry data with invalid plane ID: " << plane << std::endl;
                        continue; 
                    }
                    auto const &dqdx_values = calo->dQdx();
                    auto const &dedx_values = calo->dEdx();
                    auto const &rr = calo->ResidualRange();
                    auto const &pitch = calo->TrkPitchVec();
                    auto const& xyz_v = calo->XYZ();
                    std::vector<std::vector<float>> par_values;
                    par_values.push_back(rr);
                    par_values.push_back(pitch);

                    float calo_energy = 0;
                    std::vector<float> dqdx_values_corrected, dedx_values_corrected;

                    if (fData || !fRecalibrateHits)
                    {
                        dqdx_values_corrected = dqdx_values;
                    }
                    else
                    {
                        dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, trk.value(), assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, false);
                    }

                    for (size_t i = 0; i < dqdx_values_corrected.size(); i++)
                    {
                        float aux_dedx;
                        aux_dedx = common::ModBoxCorrection(dqdx_values_corrected[i]*fADCtoE[plane], xyz_v[i].X(), xyz_v[i].Y(), xyz_v[i].Z());
                        dedx_values_corrected.push_back(aux_dedx);
                        calo_energy += aux_dedx * pitch[i];
                    }

                    float trk_nhits = dedx_values.size();
                    float trk_trunk_dEdx = CalculateTrackTrunkdEdxByHits(dedx_values);
                    float trk_trunk_rr_dEdx = CalculateTrackTrunkdEdxByRange(dedx_values, rr);

                    float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);

                    if (plane == 0)
                    {
                        _trk_llr_pid_u_v.back() = llr_pid;
                        _trk_calo_energy_u_v.back() = calo_energy;
                        _trk_nhits_u_v.back() = trk_nhits;
                        _trk_trunk_dEdx_u_v.back() = trk_trunk_dEdx;
                        _trk_trunk_rr_dEdx_u_v.back() = trk_trunk_rr_dEdx;
                    }
                    else if (plane == 1)
                    {
                        _trk_llr_pid_v_v.back() = llr_pid;
                        _trk_calo_energy_v_v.back() = calo_energy;
                        _trk_nhits_v_v.back() = trk_nhits;
                        _trk_trunk_dEdx_v_v.back() = trk_trunk_dEdx;
                        _trk_trunk_rr_dEdx_v_v.back() = trk_trunk_rr_dEdx;
                    }
                    else if (plane == 2)
                    {
                        _trk_llr_pid_y_v.back() = llr_pid;
                        _trk_calo_energy_y_v.back() = calo_energy;
                        _trk_nhits_y_v.back() = trk_nhits;
                        _trk_trunk_dEdx_y_v.back() = trk_trunk_dEdx;
                        _trk_trunk_rr_dEdx_y_v.back() = trk_trunk_rr_dEdx;
                    }
                    _trk_llr_pid_v.back() += llr_pid;
                }
                _trk_llr_pid_score_v.back() = atan(_trk_llr_pid_v.back() / 100.) * 2 / 3.14159266;

                CalculateTrackDeflections(trk, _trk_avg_deflection_mean_v, _trk_avg_deflection_stdev_v, _trk_avg_deflection_separation_mean_v);

                int nPoints = 0;
                float distSquared = fEndSpacepointDistance*fEndSpacepointDistance;
                TVector3 trkEnd(_trk_end_sce[0], _trk_end_sce[1], _trk_end_sce[2]);
                for (auto &sp : spacePointCollection) {
                    float _sp_sce[3];
                    common::ApplySCECorrectionXYZ(sp->XYZ()[0], sp->XYZ()[1], sp->XYZ()[2], _sp_sce);
                    TVector3 spacePoint(_sp_sce[0], _sp_sce[1], _sp_sce[2]);
                    if ((trkEnd - spacePoint).Mag2() < distSquared) nPoints++;
                }
                _trk_end_spacepoints_v.push_back(nPoints);
            }
            else
            {
                fillDefault();
            }
        }

        std::cout << "finish track slice" << std::endl;
    }

    void TrackAnalysis::fillDefault()
    {
        _trk_pfp_id_v.push_back(std::numeric_limits<int>::lowest());

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

        _trk_llr_pid_u_v.push_back(std::numeric_limits<float>::lowest());
        _trk_llr_pid_v_v.push_back(std::numeric_limits<float>::lowest());
        _trk_llr_pid_y_v.push_back(std::numeric_limits<float>::lowest());
        _trk_llr_pid_v.push_back(std::numeric_limits<float>::lowest());
        _trk_llr_pid_score_v.push_back(std::numeric_limits<float>::lowest());

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
        _tree->Branch("trk_bragg_p_v", "std::vector< float >", &_trk_bragg_p_v);
        _tree->Branch("trk_bragg_mu_v", "std::vector< float >", &_trk_bragg_mu_v);
        _tree->Branch("trk_bragg_pion_v", "std::vector< float >", &_trk_bragg_pion_v);
        _tree->Branch("trk_bragg_mip_v", "std::vector< float >", &_trk_bragg_mip_v);
        _tree->Branch("trk_bragg_p_alt_dir_v", "std::vector< float >", &_trk_bragg_p_alt_dir_v);
        _tree->Branch("trk_bragg_mu_alt_dir_v", "std::vector< float >", &_trk_bragg_mu_alt_dir_v);
        _tree->Branch("trk_bragg_pion_alt_dir_v", "std::vector< float >", &_trk_bragg_pion_alt_dir_v);
        _tree->Branch("trk_bragg_p_fwd_preferred_v", "std::vector< bool >", &_trk_bragg_p_fwd_preferred_v);
        _tree->Branch("trk_bragg_mu_fwd_preferred_v", "std::vector< bool >", &_trk_bragg_mu_fwd_preferred_v);
        _tree->Branch("trk_bragg_pion_fwd_preferred_v", "std::vector< bool >", &_trk_bragg_pion_fwd_preferred_v);
        _tree->Branch("trk_pida_v", "std::vector< float >", &_trk_pida_v);
        _tree->Branch("trk_pid_chipr_v", "std::vector< float >", &_trk_pid_chipr_v);
        _tree->Branch("trk_pid_chipi_v", "std::vector< float >", &_trk_pid_chipi_v);
        _tree->Branch("trk_pid_chika_v", "std::vector< float >", &_trk_pid_chika_v);
        _tree->Branch("trk_pid_chimu_v", "std::vector< float >", &_trk_pid_chimu_v);

        _tree->Branch("trk_bragg_p_u_v", "std::vector< float >", &_trk_bragg_p_u_v);
        _tree->Branch("trk_bragg_mu_u_v", "std::vector< float >", &_trk_bragg_mu_u_v);
        _tree->Branch("trk_bragg_pion_u_v", "std::vector< float >", &_trk_bragg_pion_u_v);
        _tree->Branch("trk_bragg_mip_u_v", "std::vector< float >", &_trk_bragg_mip_u_v);
        _tree->Branch("trk_bragg_p_alt_dir_u_v", "std::vector< float >", &_trk_bragg_p_alt_dir_u_v);
        _tree->Branch("trk_bragg_mu_alt_dir_u_v", "std::vector< float >", &_trk_bragg_mu_alt_dir_u_v);
        _tree->Branch("trk_bragg_pion_alt_dir_u_v", "std::vector< float >", &_trk_bragg_pion_alt_dir_u_v);
        _tree->Branch("trk_bragg_p_fwd_preferred_u_v", "std::vector< bool >", &_trk_bragg_p_fwd_preferred_u_v);
        _tree->Branch("trk_bragg_mu_fwd_preferred_u_v", "std::vector< bool >", &_trk_bragg_mu_fwd_preferred_u_v);
        _tree->Branch("trk_bragg_pion_fwd_preferred_u_v", "std::vector< bool >", &_trk_bragg_pion_fwd_preferred_u_v);
        _tree->Branch("trk_pida_u_v", "std::vector< float >", &_trk_pida_u_v);
        _tree->Branch("trk_pid_chipr_u_v", "std::vector< float >", &_trk_pid_chipr_u_v);
        _tree->Branch("trk_pid_chipi_u_v", "std::vector< float >", &_trk_pid_chipi_u_v);
        _tree->Branch("trk_pid_chika_u_v", "std::vector< float >", &_trk_pid_chika_u_v);
        _tree->Branch("trk_pid_chimu_u_v", "std::vector< float >", &_trk_pid_chimu_u_v);

        _tree->Branch("trk_bragg_p_v_v", "std::vector< float >", &_trk_bragg_p_v_v);
        _tree->Branch("trk_bragg_mu_v_v", "std::vector< float >", &_trk_bragg_mu_v_v);
        _tree->Branch("trk_bragg_pion_v_v", "std::vector< float >", &_trk_bragg_pion_v_v);
        _tree->Branch("trk_bragg_mip_v_v", "std::vector< float >", &_trk_bragg_mip_v_v);
        _tree->Branch("trk_bragg_p_alt_dir_v_v", "std::vector< float >", &_trk_bragg_p_alt_dir_v_v);
        _tree->Branch("trk_bragg_mu_alt_dir_v_v", "std::vector< float >", &_trk_bragg_mu_alt_dir_v_v);
        _tree->Branch("trk_bragg_pion_alt_dir_v_v", "std::vector< float >", &_trk_bragg_pion_alt_dir_v_v);
        _tree->Branch("trk_bragg_p_fwd_preferred_v_v", "std::vector< bool >", &_trk_bragg_p_fwd_preferred_v_v);
        _tree->Branch("trk_bragg_mu_fwd_preferred_v_v", "std::vector< bool >", &_trk_bragg_mu_fwd_preferred_v_v);
        _tree->Branch("trk_bragg_pion_fwd_preferred_v_v", "std::vector< bool >", &_trk_bragg_pion_fwd_preferred_v_v);
        _tree->Branch("trk_pida_v_v", "std::vector< float >", &_trk_pida_v_v);
        _tree->Branch("trk_pid_chipr_v_v", "std::vector< float >", &_trk_pid_chipr_v_v);
        _tree->Branch("trk_pid_chipi_v_v", "std::vector< float >", &_trk_pid_chipi_v_v);
        _tree->Branch("trk_pid_chika_v_v", "std::vector< float >", &_trk_pid_chika_v_v);
        _tree->Branch("trk_pid_chimu_v_v", "std::vector< float >", &_trk_pid_chimu_v_v);

        _tree->Branch("trk_pfp_id_v", "std::vector< size_t >", &_trk_pfp_id_v);
        _tree->Branch("trk_dir_x_v", "std::vector< float >", &_trk_dir_x_v);
        _tree->Branch("trk_dir_y_v", "std::vector< float >", &_trk_dir_y_v);
        _tree->Branch("trk_dir_z_v", "std::vector< float >", &_trk_dir_z_v);

        _tree->Branch("trk_start_x_v", "std::vector< float >", &_trk_start_x_v);
        _tree->Branch("trk_start_y_v", "std::vector< float >", &_trk_start_y_v);
        _tree->Branch("trk_start_z_v", "std::vector< float >", &_trk_start_z_v);

        _tree->Branch("trk_sce_start_x_v", "std::vector< float >", &_trk_sce_start_x_v);
        _tree->Branch("trk_sce_start_y_v", "std::vector< float >", &_trk_sce_start_y_v);
        _tree->Branch("trk_sce_start_z_v", "std::vector< float >", &_trk_sce_start_z_v);

        _tree->Branch("trk_end_x_v", "std::vector< float >", &_trk_end_x_v);
        _tree->Branch("trk_end_y_v", "std::vector< float >", &_trk_end_y_v);
        _tree->Branch("trk_end_z_v", "std::vector< float >", &_trk_end_z_v);

        _tree->Branch("trk_sce_end_x_v", "std::vector< float >", &_trk_sce_end_x_v);
        _tree->Branch("trk_sce_end_y_v", "std::vector< float >", &_trk_sce_end_y_v);
        _tree->Branch("trk_sce_end_z_v", "std::vector< float >", &_trk_sce_end_z_v);

        _tree->Branch("trk_distance_v", "std::vector< float >", &_trk_distance_v);
        _tree->Branch("trk_theta_v", "std::vector< float >", &_trk_theta_v);
        _tree->Branch("trk_phi_v", "std::vector< float >", &_trk_phi_v);

        _tree->Branch("trk_len_v", "std::vector< float >", &_trk_len_v);
        _tree->Branch("trk_mcs_muon_mom_v", "std::vector< float >", &_trk_mcs_muon_mom_v);
        _tree->Branch("trk_range_muon_mom_v", "std::vector< float >", &_trk_range_muon_mom_v);
        _tree->Branch("trk_energy_proton_v", "std::vector< float >", &_trk_energy_proton_v);
        _tree->Branch("trk_energy_muon_v", "std::vector< float >", &_trk_energy_muon_v);
        _tree->Branch("trk_calo_energy_u_v", "std::vector< float >", &_trk_calo_energy_u_v);
        _tree->Branch("trk_calo_energy_v_v", "std::vector< float >", &_trk_calo_energy_v_v);
        _tree->Branch("trk_calo_energy_y_v", "std::vector< float >", &_trk_calo_energy_y_v);

        _tree->Branch("trk_llr_pid_u_v", "std::vector<float>", &_trk_llr_pid_u_v);
        _tree->Branch("trk_llr_pid_v_v", "std::vector<float>", &_trk_llr_pid_v_v);
        _tree->Branch("trk_llr_pid_y_v", "std::vector<float>", &_trk_llr_pid_y_v);
        _tree->Branch("trk_llr_pid_v", "std::vector<float>", &_trk_llr_pid_v);
        _tree->Branch("trk_llr_pid_score_v", "std::vector<float>", &_trk_llr_pid_score_v);

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

        _trk_llr_pid_u_v.clear();
        _trk_llr_pid_v_v.clear();
        _trk_llr_pid_y_v.clear();
        _trk_llr_pid_v.clear();
        _trk_llr_pid_score_v.clear();

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

        int firstHitIdx = trk_nhits - 3 - 1;

        int lastHitIdx = trk_nhits - (int)(trk_nhits/3) - 1;

        if (firstHitIdx - lastHitIdx < 5) {
            return std::numeric_limits<float>::lowest();
        }
        else {
            std::vector<float> trk_trunk_dEdx_values;
            trk_trunk_dEdx_values.reserve(firstHitIdx - lastHitIdx);

            for (int i = trk_nhits - 1; i >= 0; i--) {
                if (i > firstHitIdx) continue;

                trk_trunk_dEdx_values.push_back(dEdx_values[i]);

                if (i < lastHitIdx) break;
            }

            float median;
            std::sort(trk_trunk_dEdx_values.begin(), trk_trunk_dEdx_values.end());
            if (trk_trunk_dEdx_values.size() % 2 == 0) median = 0.5 * (trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2 - 1] + trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2]);
            else median = trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2];

            double sum = std::accumulate(std::begin(trk_trunk_dEdx_values), std::end(trk_trunk_dEdx_values), 0.0);
            double m = sum / trk_trunk_dEdx_values.size();

            double accum = 0.0;
            std::for_each(std::begin(trk_trunk_dEdx_values), std::end(trk_trunk_dEdx_values), [&](const double d) {accum += (d - m) * (d - m);});
            double stdev = sqrt(accum / (trk_trunk_dEdx_values.size()-1));

            std::vector<float> trk_trunk_dEdx_values_trimmed;
            trk_trunk_dEdx_values_trimmed.reserve(firstHitIdx - lastHitIdx);
            for (unsigned int i = 0; i < trk_trunk_dEdx_values.size(); i++) {
                if (trk_trunk_dEdx_values[i] <= median + stdev) trk_trunk_dEdx_values_trimmed.push_back(trk_trunk_dEdx_values[i]);
            }

            double sum_trimmed = std::accumulate(std::begin(trk_trunk_dEdx_values_trimmed), std::end(trk_trunk_dEdx_values_trimmed), 0.0);
            double trk_dEdx_trunk = sum_trimmed / trk_trunk_dEdx_values_trimmed.size();

            return trk_dEdx_trunk;
        }
    }

    float TrackAnalysis::CalculateTrackTrunkdEdxByRange(const std::vector<float> &dedxPerHit, const std::vector<float> &residualRangePerHit)
    {
        const auto nHitsToSkip = 3u;
        const auto lengthFraction = 1.f/3;

        if (residualRangePerHit.size() <= nHitsToSkip)
            return -std::numeric_limits<float>::max();

        std::vector<std::pair<float, unsigned int> > residualRangeIndices;
        float maxResidualRange = -std::numeric_limits<float>::max();
        for (unsigned int i = 0; i < residualRangePerHit.size(); ++i)
        {
            const auto residualRange = residualRangePerHit.at(i);
            maxResidualRange = std::max(maxResidualRange, residualRange);
            residualRangeIndices.emplace_back(residualRange, i);
        }

        const auto residualRangeCutoff = maxResidualRange * lengthFraction;

        std::sort(residualRangeIndices.begin(), residualRangeIndices.end(), [](auto &a, auto &b) {
            return a.first > b.first;
        });

        std::vector<float> dedxPerHitAtStart;
        for (unsigned int i = nHitsToSkip; i < residualRangeIndices.size(); ++i)
        {
            const auto entry = residualRangeIndices.at(i);
            const auto residualRange = entry.first;
            const auto hitIndex = entry.second;

            if (residualRange < residualRangeCutoff)
                continue;

            dedxPerHitAtStart.push_back(dedxPerHit.at(hitIndex));
        }

        const auto nHits = dedxPerHitAtStart.size();
        if (nHits == 0)
            return -std::numeric_limits<float>::max();

        std::sort(dedxPerHitAtStart.begin(), dedxPerHitAtStart.end());
        const auto median = dedxPerHitAtStart.at(nHits / 2);

        float total = 0.f;
        for (const auto &dEdx : dedxPerHitAtStart)
            total += dEdx;
        const auto mean = total / static_cast<float>(nHits);

        float squareSum = 0.f;
        for (const auto &dEdx : dedxPerHitAtStart)
            squareSum += std::pow(dEdx - mean, 2);
        const auto variance = squareSum / static_cast<float>(nHits);

        float truncatedTotal = 0.f;
        unsigned int nTruncatedHits = 0;
        for (const auto &dEdx : dedxPerHitAtStart)
        {
            if (std::pow(dEdx - median, 2) > variance)
                continue;

            truncatedTotal += dEdx;
            nTruncatedHits++;
        }

        if (nTruncatedHits == 0)
            return -std::numeric_limits<float>::max();

        return truncatedTotal / static_cast<float>(nTruncatedHits);
    }

    void TrackAnalysis::CalculateTrackDeflections(const art::Ptr<recob::Track> &trk, std::vector<float> &mean_v, std::vector<float> &stdev_v, std::vector<float> &separation_mean_v)
    {
        std::vector<size_t> validPoints;
        auto firstValidPoint = trk->FirstValidPoint();
        validPoints.push_back(firstValidPoint);
        auto nextValidPoint = trk->NextValidPoint(firstValidPoint + 1);
        while (nextValidPoint != recob::TrackTrajectory::InvalidIndex)
        {
            validPoints.push_back(nextValidPoint);
            nextValidPoint = trk->NextValidPoint(nextValidPoint + 1);
        }
        if (validPoints.size() < 3)
        {
            mean_v.push_back(0);
            stdev_v.push_back(0);
            separation_mean_v.push_back(0);
        }
        else
        {
            std::vector<float> thetaVector;
            float thetaSum = 0.f;
            float separationSum = 0.f;
            for (unsigned int i = 1; i < validPoints.size(); ++i) {
                const auto dir = trk->DirectionAtPoint(validPoints.at(i));
                const auto dirPrev = trk->DirectionAtPoint(validPoints.at(i - 1));

                const auto cosTheta = std::min(1.f, std::max(-1.f, static_cast<float>(dir.Dot(dirPrev))));
                const auto theta = std::acos(cosTheta);

                thetaSum += theta;
                thetaVector.push_back(theta);

                const TVector3 point(trk->LocationAtPoint(validPoints.at(i)).X(), trk->LocationAtPoint(validPoints.at(i)).Y(), trk->LocationAtPoint(validPoints.at(i)).Z());
                const TVector3 pointPrev(trk->LocationAtPoint(validPoints.at(i - 1)).X(), trk->LocationAtPoint(validPoints.at(i - 1)).Y(), trk->LocationAtPoint(validPoints.at(i - 1)).Z());
                const TVector3 separation = point - pointPrev;
                separationSum += separation.Mag();
            }

            float thetaMean = thetaSum / static_cast<float>(thetaVector.size());
            float separationMean = separationSum / static_cast<float>(thetaVector.size());

            float thetaDiffSum = 0.f;
            for (const auto &theta : thetaVector)
            {
                thetaDiffSum += std::pow(theta - thetaMean, 2);
            }

            const auto variance = thetaDiffSum / static_cast<float>(thetaVector.size() - 1);

            mean_v.push_back(thetaMean);
            stdev_v.push_back(std::sqrt(variance));
            separation_mean_v.push_back(separationMean);
        }
    }

    DEFINE_ART_CLASS_TOOL(TrackAnalysis)
} // namespace analysis

#endif