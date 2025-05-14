#ifndef ANALYSIS_RECONSTRUCTIONANALYSIS_CXX
#define ANALYSIS_RECONSTRUCTIONANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "simdata/SimulationBase/MCShower.h"
#include "simdata/SimulationBase/MCTrack.h"


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

#include "TTree.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <limits>
#include <cmath>
#include <numeric>

namespace analysis 
{
    class ReconstructionAnalysis : public AnalysisToolBase {
    public:
        explicit ReconstructionAnalysis(fhicl::ParameterSet const& p);
        virtual ~ReconstructionAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;
        void analyseEvent(art::Event const& e, bool is_data) override;
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fHproducer;
        art::InputTag fBacktrackTag;
        art::InputTag fMCRproducer;
        art::InputTag fCALOproducer;
        art::InputTag fPIDproducer;
        art::InputTag fTRKproducer;
        bool fRecalibrateHits;
        float fEnergyThresholdForMCHits;
        std::vector<float> fADCtoE;
        float fEndSpacepointDistance;

        const trkf::TrackMomentumCalculator _trkmom;
        trkf::TrajectoryMCSFitter _mcsfitter;
        TParticlePDG *proton_pdg;
        TParticlePDG *muon_pdg;
        common::LLRPID llr_pid_calculator;
        common::ProtonMuonLookUpParameters protonmuon_parameters;
        common::CorrectionLookUpParameters correction_parameters;

        int run, sub, evt;
        float reco_nu_vtx_x, reco_nu_vtx_y, reco_nu_vtx_z;
        float topological_score;

        std::vector<size_t> trk_pfp_id_v;
        std::vector<float> trk_start_x_v, trk_start_y_v, trk_start_z_v;
        std::vector<float> trk_sce_start_x_v, trk_sce_start_y_v, trk_sce_start_z_v;
        std::vector<float> trk_end_x_v, trk_end_y_v, trk_end_z_v;
        std::vector<float> trk_sce_end_x_v, trk_sce_end_y_v, trk_sce_end_z_v;
        std::vector<float> trk_dir_x_v, trk_dir_y_v, trk_dir_z_v;
        std::vector<float> trk_len_v;
        std::vector<float> trk_distance_v;
        std::vector<float> trk_theta_v;
        std::vector<float> trk_phi_v;
        std::vector<float> trk_bragg_p_v;
        std::vector<float> trk_bragg_mu_v;
        std::vector<float> trk_bragg_pion_v;
        std::vector<float> trk_bragg_mip_v;
        std::vector<float> trk_pid_chipr_v;
        std::vector<float> trk_pid_chika_v;
        std::vector<float> trk_pid_chipi_v;
        std::vector<float> trk_pid_chimu_v;
        std::vector<float> trk_pida_v;
        std::vector<float> trk_llr_pid_u_v;
        std::vector<float> trk_llr_pid_v_v;
        std::vector<float> trk_llr_pid_y_v;
        std::vector<float> trk_llr_pid_v;
        std::vector<float> trk_llr_pid_score_v;
        std::vector<float> trk_mcs_muon_mom_v;
        std::vector<float> trk_range_muon_mom_v;
        std::vector<float> trk_energy_proton_v;
        std::vector<float> trk_energy_muon_v;
        std::vector<float> trk_calo_energy_u_v;
        std::vector<float> trk_calo_energy_v_v;
        std::vector<float> trk_calo_energy_y_v;
        std::vector<float> trk_trunk_dEdx_u_v;
        std::vector<float> trk_trunk_dEdx_v_v;
        std::vector<float> trk_trunk_dEdx_y_v;
        std::vector<float> trk_trunk_rr_dEdx_u_v;
        std::vector<float> trk_trunk_rr_dEdx_v_v;
        std::vector<float> trk_trunk_rr_dEdx_y_v;
        std::vector<int> trk_nhits_u_v;
        std::vector<int> trk_nhits_v_v;
        std::vector<int> trk_nhits_y_v;
        std::vector<float> trk_avg_deflection_mean_v;
        std::vector<float> trk_avg_deflection_stdev_v;
        std::vector<float> trk_avg_deflection_separation_mean_v;
        std::vector<int> trk_end_spacepoints_v;
        std::vector<int> trk_bt_pdg_v;

        float CalculateTrackTrunkdEdxByHits(const std::vector<float>& dedxPerHit);
        float CalculateTrackTrunkdEdxByRange(const std::vector<float>& dedxPerHit, const std::vector<float>& residualRangePerHit);
        void CalculateTrackDeflections(const art::Ptr<recob::Track>& trk, std::vector<float>& mean_vec, std::vector<float>& stdev_vec, std::vector<float>& separation_mean_vec);
        void fillDefaultTrackVariables();
    };

    ReconstructionAnalysis::ReconstructionAnalysis(const fhicl::ParameterSet& pset) :
        _mcsfitter(fhicl::Table<trkf::TrajectoryMCSFitter::Config>(pset.get<fhicl::ParameterSet>("mcsfitmu"))) 
    {
        this->configure(pset);
        proton_pdg = TDatabasePDG::Instance()->GetParticle(2212);
        muon_pdg = TDatabasePDG::Instance()->GetParticle(13);

        llr_pid_calculator.set_dedx_binning(0, protonmuon_parameters.dedx_edges_pl_0);
        llr_pid_calculator.set_par_binning(0, protonmuon_parameters.parameters_edges_pl_0);
        llr_pid_calculator.set_lookup_tables(0, protonmuon_parameters.dedx_pdf_pl_0);
        llr_pid_calculator.set_dedx_binning(1, protonmuon_parameters.dedx_edges_pl_1);
        llr_pid_calculator.set_par_binning(1, protonmuon_parameters.parameters_edges_pl_1);
        llr_pid_calculator.set_lookup_tables(1, protonmuon_parameters.dedx_pdf_pl_1);
        llr_pid_calculator.set_dedx_binning(2, protonmuon_parameters.dedx_edges_pl_2);
        llr_pid_calculator.set_par_binning(2, protonmuon_parameters.parameters_edges_pl_2);
        llr_pid_calculator.set_lookup_tables(2, protonmuon_parameters.dedx_pdf_pl_2);
        if (fRecalibrateHits) {
            llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
            llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);
            llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
            llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);
            llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
            llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
        }
    }

    void ReconstructionAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
        fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
        fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
        fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
        fMCRproducer = p.get<art::InputTag>("MCRproducer", "mcreco");
        fCALOproducer = p.get<art::InputTag>("CALOproducer", "pandoracalo");
        fPIDproducer = p.get<art::InputTag>("PIDproducer", "pandorapid");
        fTRKproducer = p.get<art::InputTag>("TRKproducer", "pandora");
        fRecalibrateHits = p.get<bool>("RecalibrateHits", false);
        fEnergyThresholdForMCHits = p.get<float>("EnergyThresholdForMCHits", 0.1);
        fADCtoE = p.get<std::vector<float>>("ADCtoE");
        fEndSpacepointDistance = p.get<float>("EndSpacepointDistance", 5.0);
    }

    void ReconstructionAnalysis::analyseEvent(art::Event const& e, bool is_data) {
        run = e.run();
        sub = e.subRun();
        evt = e.event();
    }

    void ReconstructionAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
        reco_nu_vtx_x = std::numeric_limits<float>::lowest();
        reco_nu_vtx_y = std::numeric_limits<float>::lowest();
        reco_nu_vtx_z = std::numeric_limits<float>::lowest();
        topological_score = std::numeric_limits<float>::lowest();

        TVector3 nuvtx_internal(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());

        for (const auto& pfp : slice_pfp_v) {
            if (pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14 || std::abs(pfp->PdgCode()) == 16)) {
                auto vtx_v = pfp.get<recob::Vertex>();
                if (!vtx_v.empty()) {
                    double xyz[3] = {0.0, 0.0, 0.0};
                    vtx_v.at(0)->XYZ(xyz);
                    reco_nu_vtx_x = static_cast<float>(xyz[0]);
                    reco_nu_vtx_y = static_cast<float>(xyz[1]);
                    reco_nu_vtx_z = static_cast<float>(xyz[2]);
                    nuvtx_internal.SetXYZ(xyz[0], xyz[1], xyz[2]);
                }

                auto metadata_pxy_v = pfp.get<larpandoraobj::PFParticleMetadata>();
                if (!metadata_pxy_v.empty()) {
                    const art::Ptr<larpandoraobj::PFParticleMetadata>& pfParticleMetadata = metadata_pxy_v.at(0);
                    auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
                    if (pfParticlePropertiesMap.count("NuScore")) {
                        topological_score = pfParticlePropertiesMap.at("NuScore");
                    }
                }
                break; 
            }
        }
        
        common::ProxyCaloColl_t const* calo_proxy_ptr = nullptr;
        common::ProxyPIDColl_t const* pid_proxy_ptr = nullptr;
        try {
             calo_proxy_ptr = &proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALOproducer));
             pid_proxy_ptr = &proxy::getCollection<std::vector<recob::Track>>(e, fTRKproducer, proxy::withAssociated<anab::ParticleID>(fPIDproducer));
        } catch (const cet::exception& ex) {
            std::cerr << "WARNING: Could not get Track/Calo/PID collections: " << ex.what() << std::endl;
        }
        
        std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
        std::vector<common::BtPart> btparts_v;
        if (!is_data) {
            art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
            assocMCPart = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(inputHits, e, fBacktrackTag);
            
            art::Handle<std::vector<sim::MCShower>> inputMCShowerHandle;
            art::Handle<std::vector<sim::MCTrack>> inputMCTrackHandle;
            if (e.getByLabel(fMCRproducer, inputMCShowerHandle) && e.getByLabel(fMCRproducer, inputMCTrackHandle)) {
                 btparts_v = common::initBacktrackingParticleVec(*inputMCShowerHandle, *inputMCTrackHandle, *inputHits, assocMCPart);
            } else {
                std::cerr << "WARNING: Could not get MCShower/MCTrack collections for backtracking." << std::endl;
            }
        }

        art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
        std::vector<art::Ptr<recob::SpacePoint>> spacePointCollection;
        if (e.getByLabel(fCLSproducer, spacePointHandle)) {
            for (size_t i_sp = 0; i_sp < spacePointHandle->size(); ++i_sp) {
                spacePointCollection.emplace_back(spacePointHandle, i_sp);
            }
        }

        for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); ++i_pfp) {
            const auto& pfp = slice_pfp_v[i_pfp];
            if (pfp->IsPrimary()) continue;

            auto trk_v = pfp.get<recob::Track>();
            if (trk_v.size() == 1 && calo_proxy_ptr && pid_proxy_ptr) {
                const auto& trk = trk_v.at(0);
                trk_pfp_id_v.push_back(pfp->Self());

                auto pid_assns = (*pid_proxy_ptr)[trk.key()].get<anab::ParticleID>();
                if (!pid_assns.empty()) {
                    const auto& pid = pid_assns.at(0);
                    trk_bragg_p_v.push_back(std::max(common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 2212, 2), common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 2212, 2)));
                    trk_bragg_mu_v.push_back(std::max(common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 13, 2), common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 13, 2)));
                    trk_bragg_pion_v.push_back(std::max(common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 211, 2), common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kBackward, 211, 2)));
                    trk_bragg_mip_v.push_back(common::PID(pid, "BraggPeakLLH", anab::kLikelihood, anab::kForward, 0, 2));
                    trk_pid_chipr_v.push_back(common::PID(pid, "Chi2", anab::kGOF, anab::kForward, 2212, 2));
                    trk_pid_chika_v.push_back(common::PID(pid, "Chi2", anab::kGOF, anab::kForward, 321, 2));
                    trk_pid_chipi_v.push_back(common::PID(pid, "Chi2", anab::kGOF, anab::kForward, 211, 2));
                    trk_pid_chimu_v.push_back(common::PID(pid, "Chi2", anab::kGOF, anab::kForward, 13, 2));
                    trk_pida_v.push_back(common::PID(pid, "PIDA_mean", anab::kPIDA, anab::kForward, 0, 2));
                } else {
                    trk_bragg_p_v.push_back(std::numeric_limits<float>::lowest());
                    trk_bragg_mu_v.push_back(std::numeric_limits<float>::lowest());
                    trk_bragg_pion_v.push_back(std::numeric_limits<float>::lowest());
                    trk_bragg_mip_v.push_back(std::numeric_limits<float>::lowest());
                    trk_pid_chipr_v.push_back(std::numeric_limits<float>::lowest());
                    trk_pid_chika_v.push_back(std::numeric_limits<float>::lowest());
                    trk_pid_chipi_v.push_back(std::numeric_limits<float>::lowest());
                    trk_pid_chimu_v.push_back(std::numeric_limits<float>::lowest());
                    trk_pida_v.push_back(std::numeric_limits<float>::lowest());
                }
                
                float mcs_mom = -1.f, range_mom_mu = -1.f, energy_p = -1.f, energy_mu = -1.f;
                if (muon_pdg && proton_pdg) {
                     mcs_mom = _mcsfitter.fitMcs(trk->Trajectory(), muon_pdg->PdgCode()).bestMomentum();
                     range_mom_mu = _trkmom.GetTrackMomentum(common::GetSCECorrTrackLength(trk), muon_pdg->PdgCode());
                     float range_mom_p = _trkmom.GetTrackMomentum(common::GetSCECorrTrackLength(trk), proton_pdg->PdgCode());
                     if (range_mom_p >=0) energy_p = std::sqrt(std::pow(range_mom_p, 2) + std::pow(proton_pdg->Mass(), 2)) - proton_pdg->Mass();
                     if (range_mom_mu >=0) energy_mu = std::sqrt(std::pow(range_mom_mu, 2) + std::pow(muon_pdg->Mass(), 2)) - muon_pdg->Mass();
                }
                trk_mcs_muon_mom_v.push_back(mcs_mom);
                trk_range_muon_mom_v.push_back(range_mom_mu);
                trk_energy_proton_v.push_back(energy_p);
                trk_energy_muon_v.push_back(energy_mu);

                trk_start_x_v.push_back(trk->Start().X());
                trk_start_y_v.push_back(trk->Start().Y());
                trk_start_z_v.push_back(trk->Start().Z());
                float start_sce[3]; common::ApplySCECorrectionXYZ(trk->Start().X(), trk->Start().Y(), trk->Start().Z(), start_sce);
                trk_sce_start_x_v.push_back(start_sce[0]);
                trk_sce_start_y_v.push_back(start_sce[1]);
                trk_sce_start_z_v.push_back(start_sce[2]);
                trk_end_x_v.push_back(trk->End().X());
                trk_end_y_v.push_back(trk->End().Y());
                trk_end_z_v.push_back(trk->End().Z());
                float end_sce[3]; common::ApplySCECorrectionXYZ(trk->End().X(), trk->End().Y(), trk->End().Z(), end_sce);
                trk_sce_end_x_v.push_back(end_sce[0]);
                trk_sce_end_y_v.push_back(end_sce[1]);
                trk_sce_end_z_v.push_back(end_sce[2]);
                trk_dir_x_v.push_back(trk->StartDirection().X());
                trk_dir_y_v.push_back(trk->StartDirection().Y());
                trk_dir_z_v.push_back(trk->StartDirection().Z());
                trk_len_v.push_back(common::GetSCECorrTrackLength(trk));
                TVector3 trk_start_v3(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
                if (!std::isnan(nuvtx_internal.X())) {
                     trk_distance_v.push_back((trk_start_v3 - nuvtx_internal).Mag());
                } else {
                     trk_distance_v.push_back(std::numeric_limits<float>::lowest());
                }
                trk_theta_v.push_back(trk->Theta());
                trk_phi_v.push_back(trk->Phi());

                float llr_pid_sum = 0.f;
                float current_llr_pid_u = std::numeric_limits<float>::lowest();
                float current_llr_pid_v_plane = std::numeric_limits<float>::lowest();
                float current_llr_pid_y = std::numeric_limits<float>::lowest();
                float current_calo_e_u = 0.f, current_calo_e_v = 0.f, current_calo_e_y = 0.f;
                int current_nhits_u = 0, current_nhits_v = 0, current_nhits_y = 0;
                float current_trunk_dedx_u = std::numeric_limits<float>::lowest();
                float current_trunk_dedx_v = std::numeric_limits<float>::lowest();
                float current_trunk_dedx_y = std::numeric_limits<float>::lowest();
                float current_trunk_rr_dedx_u = std::numeric_limits<float>::lowest();
                float current_trunk_rr_dedx_v = std::numeric_limits<float>::lowest();
                float current_trunk_rr_dedx_y = std::numeric_limits<float>::lowest();

                auto calo_assns = (*calo_proxy_ptr)[trk.key()].get<anab::Calorimetry>();
                for (const auto& calo : calo_assns) {
                    if (!calo) continue;
                    auto const& plane_id = calo->PlaneID().Plane;
                    std::vector<float> dqdx_corr = (is_data || !fRecalibrateHits) ? calo->dQdx() : llr_pid_calculator.correct_many_hits_one_plane(calo, trk.value(), assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, false);
                    
                    std::vector<float> dedx_corr;
                    float current_calo_e_plane = 0.f;
                    for (size_t i_hit = 0; i_hit < dqdx_corr.size(); ++i_hit) {
                        float dedx_val = common::ModBoxCorrection(dqdx_corr[i_hit] * fADCtoE[plane_id], calo->XYZ()[i_hit].X(), calo->XYZ()[i_hit].Y(), calo->XYZ()[i_hit].Z());
                        dedx_corr.push_back(dedx_val);
                        if (calo->TrkPitchVec().size() > i_hit) current_calo_e_plane += dedx_val * calo->TrkPitchVec()[i_hit];
                    }

                    std::vector<std::vector<float>> par_values_for_llr;
                    par_values_for_llr.push_back(calo->ResidualRange());
                    par_values_for_llr.push_back(calo->TrkPitchVec());
                    float llr_plane = llr_pid_calculator.LLR_many_hits_one_plane(dedx_corr, par_values_for_llr, plane_id);
                    llr_pid_sum += llr_plane;

                    if (plane_id == 0) {
                        current_llr_pid_u = llr_plane; current_calo_e_u = current_calo_e_plane; current_nhits_u = dedx_corr.size();
                        current_trunk_dedx_u = CalculateTrackTrunkdEdxByHits(dedx_corr);
                        current_trunk_rr_dedx_u = CalculateTrackTrunkdEdxByRange(dedx_corr, calo->ResidualRange());
                    } else if (plane_id == 1) {
                        current_llr_pid_v_plane = llr_plane; current_calo_e_v = current_calo_e_plane; current_nhits_v = dedx_corr.size();
                        current_trunk_dedx_v = CalculateTrackTrunkdEdxByHits(dedx_corr);
                        current_trunk_rr_dedx_v = CalculateTrackTrunkdEdxByRange(dedx_corr, calo->ResidualRange());
                    } else if (plane_id == 2) {
                        current_llr_pid_y = llr_plane; current_calo_e_y = current_calo_e_plane; current_nhits_y = dedx_corr.size();
                        current_trunk_dedx_y = CalculateTrackTrunkdEdxByHits(dedx_corr);
                        current_trunk_rr_dedx_y = CalculateTrackTrunkdEdxByRange(dedx_corr, calo->ResidualRange());
                    }
                }
                trk_llr_pid_u_v.push_back(current_llr_pid_u);
                trk_llr_pid_v_v.push_back(current_llr_pid_v_plane);
                trk_llr_pid_y_v.push_back(current_llr_pid_y);
                trk_llr_pid_v.push_back(llr_pid_sum);
                trk_llr_pid_score_v.push_back(std::atan(llr_pid_sum / 100.f) * 2.f / M_PI);
                trk_calo_energy_u_v.push_back(current_calo_e_u);
                trk_calo_energy_v_v.push_back(current_calo_e_v);
                trk_calo_energy_y_v.push_back(current_calo_e_y);
                trk_nhits_u_v.push_back(current_nhits_u);
                trk_nhits_v_v.push_back(current_nhits_v);
                trk_nhits_y_v.push_back(current_nhits_y);
                trk_trunk_dEdx_u_v.push_back(current_trunk_dedx_u);
                trk_trunk_dEdx_v_v.push_back(current_trunk_dedx_v);
                trk_trunk_dEdx_y_v.push_back(current_trunk_dedx_y);
                trk_trunk_rr_dEdx_u_v.push_back(current_trunk_rr_dedx_u);
                trk_trunk_rr_dEdx_v_v.push_back(current_trunk_rr_dedx_v);
                trk_trunk_rr_dEdx_y_v.push_back(current_trunk_rr_dedx_y);

                std::vector<float> avg_defl_mean, avg_defl_stdev, avg_defl_sep_mean;
                CalculateTrackDeflections(trk, avg_defl_mean, avg_defl_stdev, avg_defl_sep_mean);
                trk_avg_deflection_mean_v.insert(trk_avg_deflection_mean_v.end(), avg_defl_mean.begin(), avg_defl_mean.end());
                trk_avg_deflection_stdev_v.insert(trk_avg_deflection_stdev_v.end(), avg_defl_stdev.begin(), avg_defl_stdev.end());
                trk_avg_deflection_separation_mean_v.insert(trk_avg_deflection_separation_mean_v.end(), avg_defl_sep_mean.begin(), avg_defl_sep_mean.end());

                int n_end_sp = 0;
                if (spacePointHandle.isValid()) {
                    TVector3 trkEndSCE(end_sce[0], end_sce[1], end_sce[2]);
                    float distSqCut = fEndSpacepointDistance * fEndSpacepointDistance;
                    for (const auto& sp_ptr : spacePointCollection) {
                        float sp_xyz_sce[3];
                        common::ApplySCECorrectionXYZ(sp_ptr->XYZ()[0], sp_ptr->XYZ()[1], sp_ptr->XYZ()[2], sp_xyz_sce);
                        TVector3 sp_v3_sce(sp_xyz_sce[0], sp_xyz_sce[1], sp_xyz_sce[2]);
                        if ((trkEndSCE - sp_v3_sce).Mag2() < distSqCut) {
                            n_end_sp++;
                        }
                    }
                }
                trk_end_spacepoints_v.push_back(n_end_sp);
                
                int bt_pdg = 0;
                if (!is_data && assocMCPart) {
                    auto trkHits = (*calo_proxy_ptr)[trk.key()].get<recob::Hit>();
                    float purity = 0.0, completeness = 0.0, overlay_purity = 0.0;
                    int ibt = common::getAssocBtPart(trkHits, assocMCPart, btparts_v, purity, completeness, overlay_purity);
                    if (ibt >=0 && static_cast<size_t>(ibt) < btparts_v.size()) {
                        bt_pdg = btparts_v[ibt].part->PdgCode();
                    }
                }
                trk_bt_pdg_v.push_back(bt_pdg);

            } else {
                fillDefaultTrackVariables();
            }
        }
    }

    void ReconstructionAnalysis::fillDefaultTrackVariables() {
        trk_pfp_id_v.push_back(std::numeric_limits<size_t>::max());
        trk_start_x_v.push_back(std::numeric_limits<float>::lowest()); trk_start_y_v.push_back(std::numeric_limits<float>::lowest()); trk_start_z_v.push_back(std::numeric_limits<float>::lowest());
        trk_sce_start_x_v.push_back(std::numeric_limits<float>::lowest()); trk_sce_start_y_v.push_back(std::numeric_limits<float>::lowest()); trk_sce_start_z_v.push_back(std::numeric_limits<float>::lowest());
        trk_end_x_v.push_back(std::numeric_limits<float>::lowest()); trk_end_y_v.push_back(std::numeric_limits<float>::lowest()); trk_end_z_v.push_back(std::numeric_limits<float>::lowest());
        trk_sce_end_x_v.push_back(std::numeric_limits<float>::lowest()); trk_sce_end_y_v.push_back(std::numeric_limits<float>::lowest()); trk_sce_end_z_v.push_back(std::numeric_limits<float>::lowest());
        trk_dir_x_v.push_back(std::numeric_limits<float>::lowest()); trk_dir_y_v.push_back(std::numeric_limits<float>::lowest()); trk_dir_z_v.push_back(std::numeric_limits<float>::lowest());
        trk_len_v.push_back(std::numeric_limits<float>::lowest());
        trk_distance_v.push_back(std::numeric_limits<float>::lowest());
        trk_theta_v.push_back(std::numeric_limits<float>::lowest());
        trk_phi_v.push_back(std::numeric_limits<float>::lowest());
        trk_bragg_p_v.push_back(std::numeric_limits<float>::lowest()); trk_bragg_mu_v.push_back(std::numeric_limits<float>::lowest()); trk_bragg_pion_v.push_back(std::numeric_limits<float>::lowest()); trk_bragg_mip_v.push_back(std::numeric_limits<float>::lowest());
        trk_pid_chipr_v.push_back(std::numeric_limits<float>::lowest()); trk_pid_chika_v.push_back(std::numeric_limits<float>::lowest()); trk_pid_chipi_v.push_back(std::numeric_limits<float>::lowest()); trk_pid_chimu_v.push_back(std::numeric_limits<float>::lowest());
        trk_pida_v.push_back(std::numeric_limits<float>::lowest());
        trk_llr_pid_u_v.push_back(std::numeric_limits<float>::lowest()); trk_llr_pid_v_v.push_back(std::numeric_limits<float>::lowest()); trk_llr_pid_y_v.push_back(std::numeric_limits<float>::lowest());
        trk_llr_pid_v.push_back(std::numeric_limits<float>::lowest()); trk_llr_pid_score_v.push_back(std::numeric_limits<float>::lowest());
        trk_mcs_muon_mom_v.push_back(std::numeric_limits<float>::lowest()); trk_range_muon_mom_v.push_back(std::numeric_limits<float>::lowest());
        trk_energy_proton_v.push_back(std::numeric_limits<float>::lowest()); trk_energy_muon_v.push_back(std::numeric_limits<float>::lowest());
        trk_calo_energy_u_v.push_back(std::numeric_limits<float>::lowest()); trk_calo_energy_v_v.push_back(std::numeric_limits<float>::lowest()); trk_calo_energy_y_v.push_back(std::numeric_limits<float>::lowest());
        trk_trunk_dEdx_u_v.push_back(std::numeric_limits<float>::lowest()); trk_trunk_dEdx_v_v.push_back(std::numeric_limits<float>::lowest()); trk_trunk_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());
        trk_trunk_rr_dEdx_u_v.push_back(std::numeric_limits<float>::lowest()); trk_trunk_rr_dEdx_v_v.push_back(std::numeric_limits<float>::lowest()); trk_trunk_rr_dEdx_y_v.push_back(std::numeric_limits<float>::lowest());
        trk_nhits_u_v.push_back(std::numeric_limits<int>::lowest()); trk_nhits_v_v.push_back(std::numeric_limits<int>::lowest()); trk_nhits_y_v.push_back(std::numeric_limits<int>::lowest());
        trk_avg_deflection_mean_v.push_back(std::numeric_limits<float>::lowest()); trk_avg_deflection_stdev_v.push_back(std::numeric_limits<float>::lowest()); trk_avg_deflection_separation_mean_v.push_back(std::numeric_limits<float>::lowest());
        trk_end_spacepoints_v.push_back(std::numeric_limits<int>::lowest());
        trk_bt_pdg_v.push_back(std::numeric_limits<int>::lowest());
    }


    void ReconstructionAnalysis::setBranches(TTree* tree) {
        if (!tree) return;
        tree->Branch("run", &run, "run/I");
        tree->Branch("sub", &sub, "sub/I");
        tree->Branch("evt", &evt, "evt/I");
        tree->Branch("reco_nu_vtx_x", &reco_nu_vtx_x, "reco_nu_vtx_x/F");
        tree->Branch("reco_nu_vtx_y", &reco_nu_vtx_y, "reco_nu_vtx_y/F");
        tree->Branch("reco_nu_vtx_z", &reco_nu_vtx_z, "reco_nu_vtx_z/F");
        tree->Branch("topological_score", &topological_score, "topological_score/F");

        tree->Branch("trk_pfp_id_v", &trk_pfp_id_v);
        tree->Branch("trk_start_x_v", &trk_start_x_v); tree->Branch("trk_start_y_v", &trk_start_y_v); tree->Branch("trk_start_z_v", &trk_start_z_v);
        tree->Branch("trk_sce_start_x_v", &trk_sce_start_x_v); tree->Branch("trk_sce_start_y_v", &trk_sce_start_y_v); tree->Branch("trk_sce_start_z_v", &trk_sce_start_z_v);
        tree->Branch("trk_end_x_v", &trk_end_x_v); tree->Branch("trk_end_y_v", &trk_end_y_v); tree->Branch("trk_end_z_v", &trk_end_z_v);
        tree->Branch("trk_sce_end_x_v", &trk_sce_end_x_v); tree->Branch("trk_sce_end_y_v", &trk_sce_end_y_v); tree->Branch("trk_sce_end_z_v", &trk_sce_end_z_v);
        tree->Branch("trk_dir_x_v", &trk_dir_x_v); tree->Branch("trk_dir_y_v", &trk_dir_y_v); tree->Branch("trk_dir_z_v", &trk_dir_z_v);
        tree->Branch("trk_len_v", &trk_len_v);
        tree->Branch("trk_distance_v", &trk_distance_v);
        tree->Branch("trk_theta_v", &trk_theta_v);
        tree->Branch("trk_phi_v", &trk_phi_v);
        tree->Branch("trk_bragg_p_v", &trk_bragg_p_v); tree->Branch("trk_bragg_mu_v", &trk_bragg_mu_v); tree->Branch("trk_bragg_pion_v", &trk_bragg_pion_v); tree->Branch("trk_bragg_mip_v", &trk_bragg_mip_v);
        tree->Branch("trk_pid_chipr_v", &trk_pid_chipr_v); tree->Branch("trk_pid_chika_v", &trk_pid_chika_v); tree->Branch("trk_pid_chipi_v", &trk_pid_chipi_v); tree->Branch("trk_pid_chimu_v", &trk_pid_chimu_v);
        tree->Branch("trk_pida_v", &trk_pida_v);
        tree->Branch("trk_llr_pid_u_v", &trk_llr_pid_u_v); tree->Branch("trk_llr_pid_v_v", &trk_llr_pid_v_v); tree->Branch("trk_llr_pid_y_v", &trk_llr_pid_y_v);
        tree->Branch("trk_llr_pid_v", &trk_llr_pid_v); tree->Branch("trk_llr_pid_score_v", &trk_llr_pid_score_v);
        tree->Branch("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v); tree->Branch("trk_range_muon_mom_v", &trk_range_muon_mom_v);
        tree->Branch("trk_energy_proton_v", &trk_energy_proton_v); tree->Branch("trk_energy_muon_v", &trk_energy_muon_v);
        tree->Branch("trk_calo_energy_u_v", &trk_calo_energy_u_v); tree->Branch("trk_calo_energy_v_v", &trk_calo_energy_v_v); tree->Branch("trk_calo_energy_y_v", &trk_calo_energy_y_v);
        tree->Branch("trk_trunk_dEdx_u_v", &trk_trunk_dEdx_u_v); tree->Branch("trk_trunk_dEdx_v_v", &trk_trunk_dEdx_v_v); tree->Branch("trk_trunk_dEdx_y_v", &trk_trunk_dEdx_y_v);
        tree->Branch("trk_trunk_rr_dEdx_u_v", &trk_trunk_rr_dEdx_u_v); tree->Branch("trk_trunk_rr_dEdx_v_v", &trk_trunk_rr_dEdx_v_v); tree->Branch("trk_trunk_rr_dEdx_y_v", &trk_trunk_rr_dEdx_y_v);
        tree->Branch("trk_nhits_u_v", &trk_nhits_u_v); tree->Branch("trk_nhits_v_v", &trk_nhits_v_v); tree->Branch("trk_nhits_y_v", &trk_nhits_y_v);
        tree->Branch("trk_avg_deflection_mean_v", &trk_avg_deflection_mean_v); tree->Branch("trk_avg_deflection_stdev_v", &trk_avg_deflection_stdev_v); tree->Branch("trk_avg_deflection_separation_mean_v", &trk_avg_deflection_separation_mean_v);
        tree->Branch("trk_end_spacepoints_v", &trk_end_spacepoints_v);
        tree->Branch("trk_bt_pdg_v", &trk_bt_pdg_v);
    }

    void ReconstructionAnalysis::resetTTree(TTree* tree) {
        run = -1; sub = -1; evt = -1;
        reco_nu_vtx_x = std::numeric_limits<float>::lowest(); reco_nu_vtx_y = std::numeric_limits<float>::lowest(); reco_nu_vtx_z = std::numeric_limits<float>::lowest();
        topological_score = std::numeric_limits<float>::lowest();

        trk_pfp_id_v.clear();
        trk_start_x_v.clear(); trk_start_y_v.clear(); trk_start_z_v.clear();
        trk_sce_start_x_v.clear(); trk_sce_start_y_v.clear(); trk_sce_start_z_v.clear();
        trk_end_x_v.clear(); trk_end_y_v.clear(); trk_end_z_v.clear();
        trk_sce_end_x_v.clear(); trk_sce_end_y_v.clear(); trk_sce_end_z_v.clear();
        trk_dir_x_v.clear(); trk_dir_y_v.clear(); trk_dir_z_v.clear();
        trk_len_v.clear();
        trk_distance_v.clear();
        trk_theta_v.clear();
        trk_phi_v.clear();
        trk_bragg_p_v.clear(); trk_bragg_mu_v.clear(); trk_bragg_pion_v.clear(); trk_bragg_mip_v.clear();
        trk_pid_chipr_v.clear(); trk_pid_chika_v.clear(); trk_pid_chipi_v.clear(); trk_pid_chimu_v.clear();
        trk_pida_v.clear();
        trk_llr_pid_u_v.clear(); trk_llr_pid_v_v.clear(); trk_llr_pid_y_v.clear();
        trk_llr_pid_v.clear(); trk_llr_pid_score_v.clear();
        trk_mcs_muon_mom_v.clear(); trk_range_muon_mom_v.clear();
        trk_energy_proton_v.clear(); trk_energy_muon_v.clear();
        trk_calo_energy_u_v.clear(); trk_calo_energy_v_v.clear(); trk_calo_energy_y_v.clear();
        trk_trunk_dEdx_u_v.clear(); trk_trunk_dEdx_v_v.clear(); trk_trunk_dEdx_y_v.clear();
        trk_trunk_rr_dEdx_u_v.clear(); trk_trunk_rr_dEdx_v_v.clear(); trk_trunk_rr_dEdx_y_v.clear();
        trk_nhits_u_v.clear(); trk_nhits_v_v.clear(); trk_nhits_y_v.clear();
        trk_avg_deflection_mean_v.clear(); trk_avg_deflection_stdev_v.clear(); trk_avg_deflection_separation_mean_v.clear();
        trk_end_spacepoints_v.clear();
        trk_bt_pdg_v.clear();
    }

    float ReconstructionAnalysis::CalculateTrackTrunkdEdxByHits(const std::vector<float> &dEdx_values) {
        if (dEdx_values.empty()) return std::numeric_limits<float>::lowest();
        unsigned int trk_nhits = dEdx_values.size();
        if (trk_nhits <= 3) return (trk_nhits > 0) ? std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0f) / trk_nhits : std::numeric_limits<float>::lowest();

        int firstHitIdx = trk_nhits - 3 - 1; 
        int lastHitIdx = trk_nhits - (int)(trk_nhits/3.0) - 1; 
        
        if (firstHitIdx < 0) firstHitIdx = 0;
        if (lastHitIdx < 0 || lastHitIdx < firstHitIdx) { 
            if (trk_nhits > 0) return std::accumulate(dEdx_values.begin(), dEdx_values.end(), 0.0f) / trk_nhits;
            return std::numeric_limits<float>::lowest();
        }


        std::vector<float> trk_trunk_dEdx_values;
        for (int i = trk_nhits - 1; i >= 0; --i) {
            if (i > firstHitIdx) continue; 
            if (i < lastHitIdx) break;    
            trk_trunk_dEdx_values.push_back(dEdx_values[i]);
        }

        if (trk_trunk_dEdx_values.size() < 5) { 
             if (!trk_trunk_dEdx_values.empty()) return std::accumulate(trk_trunk_dEdx_values.begin(), trk_trunk_dEdx_values.end(), 0.0f) / trk_trunk_dEdx_values.size();
             return std::numeric_limits<float>::lowest();
        }
        
        std::sort(trk_trunk_dEdx_values.begin(), trk_trunk_dEdx_values.end());
        float median = (trk_trunk_dEdx_values.size() % 2 != 0) ?
                       trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2] :
                       0.5f * (trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2 - 1] + trk_trunk_dEdx_values[trk_trunk_dEdx_values.size()/2]);
        
        double sum = std::accumulate(trk_trunk_dEdx_values.begin(), trk_trunk_dEdx_values.end(), 0.0);
        double mean_val = sum / trk_trunk_dEdx_values.size();
        double accum_sq_diff = 0.0;
        for(const double d : trk_trunk_dEdx_values) { accum_sq_diff += (d - mean_val) * (d - mean_val); }
        double stdev = (trk_trunk_dEdx_values.size() > 1) ? std::sqrt(accum_sq_diff / (trk_trunk_dEdx_values.size()-1)) : 0.0;

        std::vector<float> trk_trunk_dEdx_values_trimmed;
        for (float val : trk_trunk_dEdx_values) {
            if (std::abs(val - median) <= stdev) { 
                trk_trunk_dEdx_values_trimmed.push_back(val);
            }
        }

        if (trk_trunk_dEdx_values_trimmed.empty()) return std::numeric_limits<float>::lowest();
        double sum_trimmed = std::accumulate(trk_trunk_dEdx_values_trimmed.begin(), trk_trunk_dEdx_values_trimmed.end(), 0.0);
        return static_cast<float>(sum_trimmed / trk_trunk_dEdx_values_trimmed.size());
    }

    float ReconstructionAnalysis::CalculateTrackTrunkdEdxByRange(const std::vector<float> &dedxPerHit, const std::vector<float> &residualRangePerHit) {
        const auto nHitsToSkip = 3u;
        const float lengthFraction = 1.f/3.f;

        if (residualRangePerHit.size() <= nHitsToSkip || dedxPerHit.size() != residualRangePerHit.size())
            return std::numeric_limits<float>::lowest();

        std::vector<std::pair<float, unsigned int>> residualRangeIndices;
        float maxResidualRange = std::numeric_limits<float>::lowest();
        for (unsigned int i = 0; i < residualRangePerHit.size(); ++i) {
            const auto rr = residualRangePerHit.at(i);
            if (rr > maxResidualRange) maxResidualRange = rr;
            residualRangeIndices.emplace_back(rr, i);
        }

        if (maxResidualRange <= 0.f) return std::numeric_limits<float>::lowest();

        const auto residualRangeCutoff = maxResidualRange * lengthFraction;
        std::sort(residualRangeIndices.begin(), residualRangeIndices.end(), [](const auto &a, const auto &b) {
            return a.first > b.first;
        });

        std::vector<float> dedxPerHitAtStart;
        for (unsigned int i = nHitsToSkip; i < residualRangeIndices.size(); ++i) {
            const auto& entry = residualRangeIndices.at(i);
            if (entry.first < residualRangeCutoff) continue;
            dedxPerHitAtStart.push_back(dedxPerHit.at(entry.second));
        }

        if (dedxPerHitAtStart.empty()) return std::numeric_limits<float>::lowest();
        
        std::sort(dedxPerHitAtStart.begin(), dedxPerHitAtStart.end());
        const float median = dedxPerHitAtStart.at(dedxPerHitAtStart.size() / 2);
        float total = 0.f;
        for (const auto &dEdx : dedxPerHitAtStart) total += dEdx;
        const float mean = total / static_cast<float>(dedxPerHitAtStart.size());
        
        float squareSum = 0.f;
        for (const auto &dEdx : dedxPerHitAtStart) squareSum += std::pow(dEdx - mean, 2);
        const float variance = (dedxPerHitAtStart.size() > 1) ? (squareSum / static_cast<float>(dedxPerHitAtStart.size() -1)) : 0.f;
        
        float truncatedTotal = 0.f;
        unsigned int nTruncatedHits = 0;
        for (const auto &dEdx : dedxPerHitAtStart) {
            if (std::pow(dEdx - median, 2) > variance) continue;
            truncatedTotal += dEdx;
            nTruncatedHits++;
        }

        if (nTruncatedHits == 0) return std::numeric_limits<float>::lowest();
        return truncatedTotal / static_cast<float>(nTruncatedHits);
    }


    void ReconstructionAnalysis::CalculateTrackDeflections(const art::Ptr<recob::Track> &trk, std::vector<float> &mean_vec, std::vector<float> &stdev_vec, std::vector<float> &separation_mean_vec) {
        std::vector<size_t> validPoints;
        auto firstValidPoint = trk->FirstValidPoint();
        if (firstValidPoint == recob::TrackTrajectory::InvalidIndex) {
            mean_vec.push_back(0); stdev_vec.push_back(0); separation_mean_vec.push_back(0);
            return;
        }
        validPoints.push_back(firstValidPoint);
        auto nextValidPoint = trk->NextValidPoint(firstValidPoint + 1);
        while (nextValidPoint != recob::TrackTrajectory::InvalidIndex) {
            validPoints.push_back(nextValidPoint);
            nextValidPoint = trk->NextValidPoint(nextValidPoint + 1);
        }

        if (validPoints.size() < 2) {
            mean_vec.push_back(0); stdev_vec.push_back(0); separation_mean_vec.push_back(0);
            return;
        }
        
        if (validPoints.size() < 3) { 
             mean_vec.push_back(0); stdev_vec.push_back(0);
             const TVector3 p0(trk->LocationAtPoint(validPoints.at(0)).X(), trk->LocationAtPoint(validPoints.at(0)).Y(), trk->LocationAtPoint(validPoints.at(0)).Z());
             const TVector3 p1(trk->LocationAtPoint(validPoints.at(1)).X(), trk->LocationAtPoint(validPoints.at(1)).Y(), trk->LocationAtPoint(validPoints.at(1)).Z());
             separation_mean_vec.push_back((p1-p0).Mag());
             return;
        }


        std::vector<float> thetaVector;
        float thetaSum = 0.f;
        float separationSum = 0.f;
        for (unsigned int i = 1; i < validPoints.size() -1; ++i) { 
            const auto dirCurrent = trk->DirectionAtPoint(validPoints.at(i));
            const auto dirNext = trk->DirectionAtPoint(validPoints.at(i+1)); 
            
            if(dirCurrent.Mag2() < 1e-6 || dirNext.Mag2() < 1e-6) continue;

            const double cosTheta = dirCurrent.Dot(dirNext) / (dirCurrent.Mag() * dirNext.Mag());
            const float theta = std::acos(std::min(1.0, std::max(-1.0, cosTheta))); 
            thetaSum += theta;
            thetaVector.push_back(theta);
            
            const TVector3 pointCurrent(trk->LocationAtPoint(validPoints.at(i)).X(), trk->LocationAtPoint(validPoints.at(i)).Y(), trk->LocationAtPoint(validPoints.at(i)).Z());
            const TVector3 pointPrev(trk->LocationAtPoint(validPoints.at(i-1)).X(), trk->LocationAtPoint(validPoints.at(i-1)).Y(), trk->LocationAtPoint(validPoints.at(i-1)).Z());
            separationSum += (pointCurrent - pointPrev).Mag();
        }
        if (validPoints.size() >=2) {
             const TVector3 plast(trk->LocationAtPoint(validPoints.at(validPoints.size()-1)).X(), trk->LocationAtPoint(validPoints.at(validPoints.size()-1)).Y(), trk->LocationAtPoint(validPoints.at(validPoints.size()-1)).Z());
             const TVector3 psecondlast(trk->LocationAtPoint(validPoints.at(validPoints.size()-2)).X(), trk->LocationAtPoint(validPoints.at(validPoints.size()-2)).Y(), trk->LocationAtPoint(validPoints.at(validPoints.size()-2)).Z());
             separationSum += (plast-psecondlast).Mag();
        }


        if (thetaVector.empty()) {
            mean_vec.push_back(0);
            stdev_vec.push_back(0);
        } else {
            float thetaMean = thetaSum / static_cast<float>(thetaVector.size());
            float thetaDiffSum = 0.f;
            for (const auto &theta_val : thetaVector) {
                thetaDiffSum += std::pow(theta_val - thetaMean, 2);
            }
            const auto variance = (thetaVector.size() > 1) ? (thetaDiffSum / static_cast<float>(thetaVector.size() - 1)) : 0.f;
            mean_vec.push_back(thetaMean);
            stdev_vec.push_back(std::sqrt(variance));
        }
        float separationMean = (validPoints.size() > 1) ? (separationSum / static_cast<float>(validPoints.size() -1)) : 0.f;
        separation_mean_vec.push_back(separationMean);
    }


    DEFINE_ART_CLASS_TOOL(ReconstructionAnalysis)
}
#endif
