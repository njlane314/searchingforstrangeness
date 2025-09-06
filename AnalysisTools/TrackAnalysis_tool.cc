#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"
#include "Common/EnergyCalibration.h"
#include "Common/GeometryUtils.h"
#include "Common/LLRProtonMuonLookup.h"
#include "Common/LLRCorrectionLookup.h"
#include "Common/LLRPIDCalculator.h"
#include "Common/ProxyTypes.h"
#include "Common/SpaceChargeCorrections.h"
#include "Common/TrackShowerScore.h"
#include "Common/Quantile.h"

#include <cmath>  // for std::isfinite
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>

namespace analysis {

class TrackAnalysis : public AnalysisToolBase {
public:
    TrackAnalysis(const fhicl::ParameterSet& parameter_set);
    ~TrackAnalysis() {}
    void configure(const fhicl::ParameterSet& parameter_set) override;
    void analyseEvent(const art::Event& event, bool is_data) override;
    void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec, bool is_data, bool is_selected) override;
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    void fill_default();
    float calculate_track_trunk_dedx_by_hits(const std::vector<float>& dedx_values);
    float calculate_track_trunk_dedx_by_range(const std::vector<float>& dedx_values, const std::vector<float>& residual_range_values);
    float calculate_track_trunk_rr_qtrim(const std::vector<float>& dedx,
                                         const std::vector<float>& rr,
                                         float alpha=0.25f, float beta=0.75f,
                                         unsigned skip_each_end=3,
                                         double qlow=0.10, double qhigh=0.80,
                                         unsigned min_hits_in_band=8);
    void calculate_track_deflections(const art::Ptr<recob::Track>& track, std::vector<float>& mean_deflections,
                                   std::vector<float>& stdev_deflections, std::vector<float>& mean_separation);
    void calculate_track_deflections_qtrim(const art::Ptr<recob::Track>& track,
                                           std::vector<float>& mean_qtrim_out,
                                           std::vector<float>& rms_qtrim_out,
                                           std::vector<float>& q90_out,
                                           std::vector<int>& njoints_mid_out,
                                           float alpha=0.25f, float beta=0.75f,
                                           unsigned skip_each_end=2,
                                           double qlow=0.10, double qhigh=0.90,
                                           float min_seg_cm=2.0f, float max_seg_cm=15.0f,
                                           unsigned min_joints=5);

    art::InputTag fCALproducer;
    art::InputTag fTRKproducer;
    art::InputTag fHITproducer;
    art::InputTag fCLSproducer;
    art::InputTag fBKTproducer;

    std::vector<float> fADCtoE;
    float fEndSpacepointDistance;
    bool fRecalibrateHits;
    float fEnergyThresholdForMCHits;

    common::LLRPID llr_pid_calculator;

    int _run;
    int _sbr;
    int _evt;

    std::vector<size_t> _track_pfp_ids;
    std::vector<float> _track_start_x;
    std::vector<float> _track_start_y;
    std::vector<float> _track_start_z;
    std::vector<float> _track_sce_start_x;
    std::vector<float> _track_sce_start_y;
    std::vector<float> _track_sce_start_z;
    std::vector<float> _track_end_x;
    std::vector<float> _track_end_y;
    std::vector<float> _track_end_z;
    std::vector<float> _track_sce_end_x;
    std::vector<float> _track_sce_end_y;
    std::vector<float> _track_sce_end_z;
    std::vector<float> _track_direction_x;
    std::vector<float> _track_direction_y;
    std::vector<float> _track_direction_z;
    std::vector<float> _track_distance_to_vertex;
    std::vector<float> _track_theta;
    std::vector<float> _track_phi;
    std::vector<float> _track_length;
    std::vector<float> _track_calo_energy_u;
    std::vector<float> _track_calo_energy_v;
    std::vector<float> _track_calo_energy_y;
    std::vector<float> _track_trunk_dedx_u;
    std::vector<float> _track_trunk_dedx_v;
    std::vector<float> _track_trunk_dedx_y;
    std::vector<float> _track_trunk_rr_dedx_u;
    std::vector<float> _track_trunk_rr_dedx_v;
    std::vector<float> _track_trunk_rr_dedx_y;
    // --- Track-level dE/dx (for downstream quantile fit) ---
    std::vector<float> _track_tm_dedx;
    std::vector<float> _track_tm_dedx_plane_spread;
    std::vector<int>   _track_tm_dedx_nplanes;
    // --- Scatter (for downstream standardisation/fit) ---
    std::vector<int>   _track_mcs_njoints;
    std::vector<float> _track_deflection_rms;
    std::vector<float> _track_deflection_max;
    std::vector<int> _track_nhits_u;
    std::vector<int> _track_nhits_v;
    std::vector<int> _track_nhits_y;
    std::vector<float> _track_avg_deflection_mean;
    std::vector<float> _track_avg_deflection_stdev;
    std::vector<float> _track_avg_deflection_separation_mean;
    // Deflection (angles) with middle-window + quantile trimming
    std::vector<float> _track_deflection_mean_qtrim;
    std::vector<float> _track_deflection_rms_qtrim;
    std::vector<float> _track_deflection_q90;
    std::vector<int>   _track_deflection_njoints_mid;
    std::vector<int> _track_end_spacepoints;

    std::vector<float> _trk_llr_pid_u_v;
    std::vector<float> _trk_llr_pid_v_v;
    std::vector<float> _trk_llr_pid_y_v;
    std::vector<float> _trk_llr_pid_v;
    std::vector<float> _trk_llr_pid_score_v;
};

TrackAnalysis::TrackAnalysis(const fhicl::ParameterSet& parameter_set) {
    fHITproducer = parameter_set.get<art::InputTag>("HITproducer", "gaushit");
    fCLSproducer = parameter_set.get<art::InputTag>("CLSproducer", "pandora");
    fTRKproducer = parameter_set.get<art::InputTag>("TRKproducer", "pandoraTrack");
    fCALproducer = parameter_set.get<art::InputTag>("CALproducer", "pandoraTrackcali");
    fBKTproducer = parameter_set.get<art::InputTag>("BKTproducer", "gaushitTruthMatch");

    fADCtoE = parameter_set.get<std::vector<float>>("ADCtoE");
    fEndSpacepointDistance = parameter_set.get<float>("EndSpacepointDistance", 5.0);
    fRecalibrateHits = parameter_set.get<bool>("RecalibrateHits", true);
    fEnergyThresholdForMCHits = parameter_set.get<float>("EnergyThresholdForMCHits", 0.1);

    common::ProtonMuonLookUpParameters protonmuon_parameters;
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
        common::CorrectionLookUpParameters correction_parameters;
        llr_pid_calculator.set_corr_par_binning(0, correction_parameters.parameter_correction_edges_pl_0);
        llr_pid_calculator.set_correction_tables(0, correction_parameters.correction_table_pl_0);
        llr_pid_calculator.set_corr_par_binning(1, correction_parameters.parameter_correction_edges_pl_1);
        llr_pid_calculator.set_correction_tables(1, correction_parameters.correction_table_pl_1);
        llr_pid_calculator.set_corr_par_binning(2, correction_parameters.parameter_correction_edges_pl_2);
        llr_pid_calculator.set_correction_tables(2, correction_parameters.correction_table_pl_2);
    }
}

void TrackAnalysis::configure(const fhicl::ParameterSet& parameter_set) {}

void TrackAnalysis::analyseEvent(const art::Event& event, bool is_data) {
    _evt = event.event();
    _sbr = event.subRun();
    _run = event.run();
}

void TrackAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                                  bool is_data, bool is_selected) {
    const auto& calorimetry_proxy = proxy::getCollection<std::vector<recob::Track>>(
        event, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALproducer));

    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;
    if (!is_data) {
        auto inputHits = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
        assocMCPart = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(inputHits, event, fBKTproducer);
    }

    TVector3 neutrino_vertex;
    for (const auto& pfp : slice_pfp_vec) {
        if (pfp->IsPrimary()) {
            double xyz[3] = {};
            auto vertices = pfp.get<recob::Vertex>();
            if (vertices.size() != 1) {
                continue;
            } else {
                vertices.at(0)->XYZ(xyz);
                neutrino_vertex.SetXYZ(xyz[0], xyz[1], xyz[2]);
            }
            break;
        }
    }

    auto space_point_handle = event.getValidHandle<std::vector<recob::SpacePoint>>(fCLSproducer);
    std::vector<art::Ptr<recob::SpacePoint>> space_point_collection;
    for (size_t i = 0; i < space_point_handle->size(); ++i) {
        space_point_collection.emplace_back(space_point_handle, i);
    }

    for (size_t pfp_index = 0; pfp_index < slice_pfp_vec.size(); ++pfp_index) {
        const auto& pfp = slice_pfp_vec[pfp_index];
        if (pfp->IsPrimary()) continue;

        auto tracks = pfp.get<recob::Track>();
        if (tracks.size() == 1) {
            const auto& track = tracks.at(0);

            _track_pfp_ids.push_back(pfp_index);
            _track_start_x.push_back(track->Start().X());
            _track_start_y.push_back(track->Start().Y());
            _track_start_z.push_back(track->Start().Z());

            float sce_start[3];
            common::ApplySCECorrectionXYZ(track->Start().X(), track->Start().Y(), track->Start().Z(), sce_start);
            _track_sce_start_x.push_back(sce_start[0]);
            _track_sce_start_y.push_back(sce_start[1]);
            _track_sce_start_z.push_back(sce_start[2]);

            _track_end_x.push_back(track->End().X());
            _track_end_y.push_back(track->End().Y());
            _track_end_z.push_back(track->End().Z());

            float sce_end[3];
            common::ApplySCECorrectionXYZ(track->End().X(), track->End().Y(), track->End().Z(), sce_end);
            _track_sce_end_x.push_back(sce_end[0]);
            _track_sce_end_y.push_back(sce_end[1]);
            _track_sce_end_z.push_back(sce_end[2]);

            _track_direction_x.push_back(track->StartDirection().X());
            _track_direction_y.push_back(track->StartDirection().Y());
            _track_direction_z.push_back(track->StartDirection().Z());

            _track_theta.push_back(track->Theta());
            _track_phi.push_back(track->Phi());
            _track_length.push_back(common::GetSCECorrTrackLength(track));

            TVector3 track_vertex(track->Start().X(), track->Start().Y(), track->Start().Z());
            track_vertex -= neutrino_vertex;
            _track_distance_to_vertex.push_back(track_vertex.Mag());

            _trk_llr_pid_u_v.push_back(std::numeric_limits<float>::quiet_NaN());
            _trk_llr_pid_v_v.push_back(std::numeric_limits<float>::quiet_NaN());
            _trk_llr_pid_y_v.push_back(std::numeric_limits<float>::quiet_NaN());
            _trk_llr_pid_v.push_back(0.0f);
            _trk_llr_pid_score_v.push_back(std::numeric_limits<float>::quiet_NaN());

            auto calorimetry_objects = calorimetry_proxy[track.key()].get<anab::Calorimetry>();

            // Store trunk dE/dx estimates per plane to later compute a track-level metric
            float trunk_by_hits[3];
            float trunk_by_rr[3];
            for (int i = 0; i < 3; ++i) {
                trunk_by_hits[i] = std::numeric_limits<float>::lowest();
                trunk_by_rr[i]   = std::numeric_limits<float>::lowest();
            }

            for (const auto& calo : calorimetry_objects) {
                int plane = calo->PlaneID().Plane;
                if (plane < 0 || plane > 2) continue;

                const auto& dqdx_values = calo->dQdx();
                const auto& residual_range = calo->ResidualRange();
                const auto& pitch_values = calo->TrkPitchVec();
                const auto& xyz_values = calo->XYZ();

                std::vector<float> dqdx_values_corrected;
                if (is_data || !fRecalibrateHits) {
                    dqdx_values_corrected = dqdx_values;
                } else {
                    dqdx_values_corrected = llr_pid_calculator.correct_many_hits_one_plane(calo, *track, assocMCPart, fRecalibrateHits, fEnergyThresholdForMCHits, false);
                }

                float calo_energy = 0.0f;
                std::vector<float> dedx_values_corrected;
                for (size_t i = 0; i < dqdx_values_corrected.size(); ++i) {
                    float dedx = common::ModBoxCorrection(dqdx_values_corrected[i] * fADCtoE[plane],
                                                         xyz_values[i].X(), xyz_values[i].Y(), xyz_values[i].Z());
                    dedx_values_corrected.push_back(dedx);
                    calo_energy += dedx * pitch_values[i];
                }

                float trk_nhits = static_cast<float>(dedx_values_corrected.size());
                float trunk_dedx = calculate_track_trunk_dedx_by_hits(dedx_values_corrected);
                float trunk_rr_dedx = calculate_track_trunk_rr_qtrim(
                    dedx_values_corrected, residual_range,
                    0.25f, 0.75f, 3, 0.10, 0.80, 8);

                std::vector<std::vector<float>> par_values;
                par_values.push_back(residual_range);
                par_values.push_back(pitch_values);
                float llr_pid = llr_pid_calculator.LLR_many_hits_one_plane(dedx_values_corrected, par_values, plane);

                if (plane == 0) {
                    _trk_llr_pid_u_v.back() = llr_pid;
                    _track_calo_energy_u.push_back(calo_energy);
                    _track_nhits_u.push_back(trk_nhits);
                    _track_trunk_dedx_u.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_u.push_back(trunk_rr_dedx);
                    trunk_by_hits[0] = trunk_dedx;
                    trunk_by_rr[0]   = trunk_rr_dedx;
                } else if (plane == 1) {
                    _trk_llr_pid_v_v.back() = llr_pid;
                    _track_calo_energy_v.push_back(calo_energy);
                    _track_nhits_v.push_back(trk_nhits);
                    _track_trunk_dedx_v.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_v.push_back(trunk_rr_dedx);
                    trunk_by_hits[1] = trunk_dedx;
                    trunk_by_rr[1]   = trunk_rr_dedx;
                } else if (plane == 2) {
                    _trk_llr_pid_y_v.back() = llr_pid;
                    _track_calo_energy_y.push_back(calo_energy);
                    _track_nhits_y.push_back(trk_nhits);
                    _track_trunk_dedx_y.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_y.push_back(trunk_rr_dedx);
                    trunk_by_hits[2] = trunk_dedx;
                    trunk_by_rr[2]   = trunk_rr_dedx;
                }
                _trk_llr_pid_v.back() += llr_pid;
            }

            auto valid = [](float x){
              return std::isfinite(x) && x != std::numeric_limits<float>::lowest();
            };

            std::vector<float> candidates;
            candidates.reserve(3);
            for (int i=0;i<3;++i) {
              if (valid(trunk_by_rr[i]))   candidates.push_back(trunk_by_rr[i]);
              else if (valid(trunk_by_hits[i])) candidates.push_back(trunk_by_hits[i]);
            }

            if (candidates.size() >= 2) {
              std::sort(candidates.begin(), candidates.end());
              float median_tm = (candidates.size()==2) ? 0.5f*(candidates[0]+candidates[1])
                                                      : candidates[candidates.size()/2];
              float plane_spread = candidates.back() - candidates.front();
              _track_tm_dedx.push_back(median_tm);
              _track_tm_dedx_plane_spread.push_back(plane_spread);
              _track_tm_dedx_nplanes.push_back(static_cast<int>(candidates.size()));
            } else {
              _track_tm_dedx.push_back(std::numeric_limits<float>::lowest());
              _track_tm_dedx_plane_spread.push_back(std::numeric_limits<float>::lowest());
              _track_tm_dedx_nplanes.push_back(std::numeric_limits<int>::lowest());
            }

            _trk_llr_pid_score_v.back() = atan(_trk_llr_pid_v.back() / 100.f) * 2 / 3.14159266f;

            calculate_track_deflections(track, _track_avg_deflection_mean, _track_avg_deflection_stdev,
                                       _track_avg_deflection_separation_mean);
            calculate_track_deflections_qtrim(track,
                                             _track_deflection_mean_qtrim,
                                             _track_deflection_rms_qtrim,
                                             _track_deflection_q90,
                                             _track_deflection_njoints_mid);

            int end_spacepoint_count = 0;
            float distance_squared = fEndSpacepointDistance * fEndSpacepointDistance;
            TVector3 track_end(_track_sce_end_x.back(), _track_sce_end_y.back(), _track_sce_end_z.back());
            for (const auto& space_point : space_point_collection) {
                float sce_space_point[3];
                common::ApplySCECorrectionXYZ(space_point->XYZ()[0], space_point->XYZ()[1], space_point->XYZ()[2], sce_space_point);
                TVector3 space_point_vector(sce_space_point[0], sce_space_point[1], sce_space_point[2]);
                if ((track_end - space_point_vector).Mag2() < distance_squared) {
                    ++end_spacepoint_count;
                }
            }
            _track_end_spacepoints.push_back(end_spacepoint_count);
        } else {
            fill_default();
        }
    }
}

void TrackAnalysis::fill_default() {
    _track_pfp_ids.push_back(std::numeric_limits<size_t>::max());
    _track_start_x.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_start_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_start_z.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_sce_start_x.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_sce_start_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_sce_start_z.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_end_x.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_end_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_end_z.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_sce_end_x.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_sce_end_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_sce_end_z.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_direction_x.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_direction_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_direction_z.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_distance_to_vertex.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_theta.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_phi.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_length.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_calo_energy_u.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_calo_energy_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_calo_energy_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_trunk_dedx_u.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_trunk_dedx_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_trunk_dedx_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_trunk_rr_dedx_u.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_trunk_rr_dedx_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_trunk_rr_dedx_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_nhits_u.push_back(-1);
    _track_nhits_v.push_back(-1);
    _track_nhits_y.push_back(-1);
    _track_avg_deflection_mean.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_avg_deflection_stdev.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_avg_deflection_separation_mean.push_back(std::numeric_limits<float>::quiet_NaN());
    _track_end_spacepoints.push_back(-1);
    _trk_llr_pid_u_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _trk_llr_pid_v_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _trk_llr_pid_y_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _trk_llr_pid_v.push_back(std::numeric_limits<float>::quiet_NaN());
    _trk_llr_pid_score_v.push_back(std::numeric_limits<float>::quiet_NaN());
}

void TrackAnalysis::setBranches(TTree* tree) {
    tree->Branch("track_pfp_ids", "std::vector<size_t>", &_track_pfp_ids);
    tree->Branch("track_start_x", "std::vector<float>", &_track_start_x);
    tree->Branch("track_start_y", "std::vector<float>", &_track_start_y);
    tree->Branch("track_start_z", "std::vector<float>", &_track_start_z);
    tree->Branch("track_sce_start_x", "std::vector<float>", &_track_sce_start_x);
    tree->Branch("track_sce_start_y", "std::vector<float>", &_track_sce_start_y);
    tree->Branch("track_sce_start_z", "std::vector<float>", &_track_sce_start_z);
    tree->Branch("track_end_x", "std::vector<float>", &_track_end_x);
    tree->Branch("track_end_y", "std::vector<float>", &_track_end_y);
    tree->Branch("track_end_z", "std::vector<float>", &_track_end_z);
    tree->Branch("track_sce_end_x", "std::vector<float>", &_track_sce_end_x);
    tree->Branch("track_sce_end_y", "std::vector<float>", &_track_sce_end_y);
    tree->Branch("track_sce_end_z", "std::vector<float>", &_track_sce_end_z);
    tree->Branch("track_direction_x", "std::vector<float>", &_track_direction_x);
    tree->Branch("track_direction_y", "std::vector<float>", &_track_direction_y);
    tree->Branch("track_direction_z", "std::vector<float>", &_track_direction_z);
    tree->Branch("track_distance_to_vertex", "std::vector<float>", &_track_distance_to_vertex);
    tree->Branch("track_theta", "std::vector<float>", &_track_theta);
    tree->Branch("track_phi", "std::vector<float>", &_track_phi);
    tree->Branch("track_length", "std::vector<float>", &_track_length);
    tree->Branch("track_calo_energy_u", "std::vector<float>", &_track_calo_energy_u);
    tree->Branch("track_calo_energy_v", "std::vector<float>", &_track_calo_energy_v);
    tree->Branch("track_calo_energy_y", "std::vector<float>", &_track_calo_energy_y);
    tree->Branch("track_trunk_dedx_u", "std::vector<float>", &_track_trunk_dedx_u);
    tree->Branch("track_trunk_dedx_v", "std::vector<float>", &_track_trunk_dedx_v);
    tree->Branch("track_trunk_dedx_y", "std::vector<float>", &_track_trunk_dedx_y);
    tree->Branch("track_trunk_rr_dedx_u", "std::vector<float>", &_track_trunk_rr_dedx_u);
    tree->Branch("track_trunk_rr_dedx_v", "std::vector<float>", &_track_trunk_rr_dedx_v);
    tree->Branch("track_trunk_rr_dedx_y", "std::vector<float>", &_track_trunk_rr_dedx_y);
    tree->Branch("track_tm_dedx", "std::vector<float>", &_track_tm_dedx);
    tree->Branch("track_tm_dedx_plane_spread", "std::vector<float>", &_track_tm_dedx_plane_spread);
    tree->Branch("track_tm_dedx_nplanes", "std::vector<int>", &_track_tm_dedx_nplanes);
    tree->Branch("track_mcs_njoints", "std::vector<int>", &_track_mcs_njoints);
    tree->Branch("track_deflection_rms", "std::vector<float>", &_track_deflection_rms);
    tree->Branch("track_deflection_max", "std::vector<float>", &_track_deflection_max);
    tree->Branch("track_nhits_u", "std::vector<int>", &_track_nhits_u);
    tree->Branch("track_nhits_v", "std::vector<int>", &_track_nhits_v);
    tree->Branch("track_nhits_y", "std::vector<int>", &_track_nhits_y);
    tree->Branch("track_avg_deflection_mean", "std::vector<float>", &_track_avg_deflection_mean);
    tree->Branch("track_avg_deflection_stdev", "std::vector<float>", &_track_avg_deflection_stdev);
    tree->Branch("track_avg_deflection_separation_mean", "std::vector<float>", &_track_avg_deflection_separation_mean);
    tree->Branch("track_deflection_mean_qtrim", "std::vector<float>", &_track_deflection_mean_qtrim);
    tree->Branch("track_deflection_rms_qtrim", "std::vector<float>", &_track_deflection_rms_qtrim);
    tree->Branch("track_deflection_q90", "std::vector<float>", &_track_deflection_q90);
    tree->Branch("track_deflection_njoints_mid", "std::vector<int>", &_track_deflection_njoints_mid);
    tree->Branch("track_end_spacepoints", "std::vector<int>", &_track_end_spacepoints);
    tree->Branch("trk_llr_pid_u", "std::vector<float>", &_trk_llr_pid_u_v);
    tree->Branch("trk_llr_pid_v", "std::vector<float>", &_trk_llr_pid_v_v);
    tree->Branch("trk_llr_pid_y", "std::vector<float>", &_trk_llr_pid_y_v);
    tree->Branch("trk_llr_pid", "std::vector<float>", &_trk_llr_pid_v);
    tree->Branch("trk_llr_pid_score", "std::vector<float>", &_trk_llr_pid_score_v);
}

void TrackAnalysis::resetTTree(TTree* tree) {
    _track_pfp_ids.clear();
    _track_start_x.clear();
    _track_start_y.clear();
    _track_start_z.clear();
    _track_sce_start_x.clear();
    _track_sce_start_y.clear();
    _track_sce_start_z.clear();
    _track_end_x.clear();
    _track_end_y.clear();
    _track_end_z.clear();
    _track_sce_end_x.clear();
    _track_sce_end_y.clear();
    _track_sce_end_z.clear();
    _track_direction_x.clear();
    _track_direction_y.clear();
    _track_direction_z.clear();
    _track_distance_to_vertex.clear();
    _track_theta.clear();
    _track_phi.clear();
    _track_length.clear();
    _track_calo_energy_u.clear();
    _track_calo_energy_v.clear();
    _track_calo_energy_y.clear();
    _track_trunk_dedx_u.clear();
    _track_trunk_dedx_v.clear();
    _track_trunk_dedx_y.clear();
    _track_trunk_rr_dedx_u.clear();
    _track_trunk_rr_dedx_v.clear();
    _track_trunk_rr_dedx_y.clear();
    _track_tm_dedx.clear();
    _track_tm_dedx_plane_spread.clear();
    _track_tm_dedx_nplanes.clear();
    _track_mcs_njoints.clear();
    _track_deflection_rms.clear();
    _track_deflection_max.clear();
    _track_nhits_u.clear();
    _track_nhits_v.clear();
    _track_nhits_y.clear();
    _track_avg_deflection_mean.clear();
    _track_avg_deflection_stdev.clear();
    _track_avg_deflection_separation_mean.clear();
    _track_deflection_mean_qtrim.clear();
    _track_deflection_rms_qtrim.clear();
    _track_deflection_q90.clear();
    _track_deflection_njoints_mid.clear();
    _track_end_spacepoints.clear();
    _trk_llr_pid_u_v.clear();
    _trk_llr_pid_v_v.clear();
    _trk_llr_pid_y_v.clear();
    _trk_llr_pid_v.clear();
    _trk_llr_pid_score_v.clear();
}

float TrackAnalysis::calculate_track_trunk_dedx_by_hits(const std::vector<float>& dedx_values) {
    size_t hit_count = dedx_values.size();
    int first_hit_index = static_cast<int>(hit_count) - 4; 
    int last_hit_index = static_cast<int>(hit_count - (hit_count / 3)) - 1;

    if (first_hit_index - last_hit_index < 5) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::vector<float> trunk_dedx_values;
    trunk_dedx_values.reserve(first_hit_index - last_hit_index + 1);
    for (int i = hit_count - 1; i >= 0 && i >= last_hit_index; --i) {
        if (i > first_hit_index) continue;
        trunk_dedx_values.push_back(dedx_values[i]);
    }

    std::sort(trunk_dedx_values.begin(), trunk_dedx_values.end());
    float median = (trunk_dedx_values.size() % 2 == 0) 
        ? 0.5f * (trunk_dedx_values[trunk_dedx_values.size() / 2 - 1] + trunk_dedx_values[trunk_dedx_values.size() / 2])
        : trunk_dedx_values[trunk_dedx_values.size() / 2];

    double sum = std::accumulate(trunk_dedx_values.begin(), trunk_dedx_values.end(), 0.0);
    double mean = sum / trunk_dedx_values.size();
    double variance_accum = 0.0;
    for (const auto& value : trunk_dedx_values) {
        variance_accum += (value - mean) * (value - mean);
    }
    double stdev = std::sqrt(variance_accum / (trunk_dedx_values.size() - 1));

    std::vector<float> trimmed_dedx_values;
    trimmed_dedx_values.reserve(trunk_dedx_values.size());
    for (const auto& value : trunk_dedx_values) {
        if (value <= median + stdev) {
            trimmed_dedx_values.push_back(value);
        }
    }

    double trimmed_sum = std::accumulate(trimmed_dedx_values.begin(), trimmed_dedx_values.end(), 0.0);
    return static_cast<float>(trimmed_sum / trimmed_dedx_values.size());
}

float TrackAnalysis::calculate_track_trunk_dedx_by_range(const std::vector<float>& dedx_values,
                                                   const std::vector<float>& residual_range_values) {
    const unsigned int hits_to_skip = 3;
    const float length_fraction = 1.0f / 3.0f;

    if (residual_range_values.size() <= hits_to_skip) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::vector<std::pair<float, unsigned int>> residual_range_indices;
    float max_residual_range = std::numeric_limits<float>::quiet_NaN();
    for (unsigned int i = 0; i < residual_range_values.size(); ++i) {
        max_residual_range = std::max(max_residual_range, residual_range_values[i]);
        residual_range_indices.emplace_back(residual_range_values[i], i);
    }

    const float residual_range_cutoff = max_residual_range * length_fraction;
    std::sort(residual_range_indices.begin(), residual_range_indices.end(), 
              [](const auto& a, const auto& b) { return a.first > b.first; });

    std::vector<float> dedx_at_start;
    for (unsigned int i = hits_to_skip; i < residual_range_indices.size(); ++i) {
        if (residual_range_indices[i].first < residual_range_cutoff) continue;
        dedx_at_start.push_back(dedx_values[residual_range_indices[i].second]);
    }

    if (dedx_at_start.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::sort(dedx_at_start.begin(), dedx_at_start.end());
    float median = dedx_at_start[dedx_at_start.size() / 2];
    float total = std::accumulate(dedx_at_start.begin(), dedx_at_start.end(), 0.0f);
    float mean = total / dedx_at_start.size();

    float square_sum = 0.0f;
    for (const auto& dedx : dedx_at_start) {
        square_sum += (dedx - mean) * (dedx - mean);
    }
    float variance = square_sum / dedx_at_start.size();

    float truncated_total = 0.0f;
    unsigned int truncated_hit_count = 0;
    for (const auto& dedx : dedx_at_start) {
        if ((dedx - median) * (dedx - median) <= variance) {
            truncated_total += dedx;
            ++truncated_hit_count;
        }
    }

    return (truncated_hit_count == 0) ? std::numeric_limits<float>::quiet_NaN() 
                                   : truncated_total / static_cast<float>(truncated_hit_count);
}

float TrackAnalysis::calculate_track_trunk_rr_qtrim(
    const std::vector<float>& dedx,
    const std::vector<float>& rr,
    float alpha, float beta,
    unsigned skip_each_end,
    double qlow, double qhigh,
    unsigned min_hits_in_band) {
    if (dedx.size() != rr.size() || dedx.empty())
        return std::numeric_limits<float>::lowest();

    std::vector<size_t> idx(rr.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
        [&](size_t a, size_t b){ return rr[a] > rr[b]; });

    float rrmax = *std::max_element(rr.begin(), rr.end());
    if (!(rrmax > 0)) return std::numeric_limits<float>::lowest();
    float lo_cut = alpha * rrmax;
    float hi_cut = beta  * rrmax;

    std::vector<float> vals;
    vals.reserve(dedx.size());
    for (size_t k = 0; k < idx.size(); ++k) {
        if (k < skip_each_end) continue;
        if (k >= idx.size() - skip_each_end) continue;
        size_t i = idx[k];
        float r = rr[i];
        if (r < lo_cut || r > hi_cut) continue;
        float v = dedx[i];
        if (std::isfinite(v)) vals.push_back(v);
    }

    if (vals.size() < min_hits_in_band)
        return std::numeric_limits<float>::lowest();

    float q_lo = common::quantile_linear(vals, qlow);
    float q_hi = common::quantile_linear(vals, qhigh);
    if (!(q_lo <= q_hi)) return std::numeric_limits<float>::lowest();

    double sum = 0.0; unsigned n = 0;
    for (float v : vals) {
        if (v < q_lo || v > q_hi) continue;
        sum += v; ++n;
    }
    if (n == 0) return std::numeric_limits<float>::lowest();
    return static_cast<float>(sum / n);
}

void TrackAnalysis::calculate_track_deflections(const art::Ptr<recob::Track>& track, 
                                              std::vector<float>& mean_deflections, 
                                              std::vector<float>& stdev_deflections, 
                                              std::vector<float>& mean_separation) {
    std::vector<size_t> valid_points;
    size_t first_valid_point = track->FirstValidPoint();
    valid_points.push_back(first_valid_point);
    size_t next_valid_point = track->NextValidPoint(first_valid_point + 1);
    while (next_valid_point != recob::TrackTrajectory::InvalidIndex) {
        valid_points.push_back(next_valid_point);
        next_valid_point = track->NextValidPoint(next_valid_point + 1);
    }

    if (valid_points.size() < 3) {
        mean_deflections.push_back(std::numeric_limits<float>::quiet_NaN());
        stdev_deflections.push_back(std::numeric_limits<float>::quiet_NaN());
        mean_separation.push_back(std::numeric_limits<float>::quiet_NaN());
        return;
    }

    std::vector<float> deflection_angles;
    float angle_sum = 0.0f;
    float separation_sum = 0.0f;
    for (size_t i = 1; i < valid_points.size(); ++i) {
        const auto current_direction = track->DirectionAtPoint(valid_points[i]);
        const auto previous_direction = track->DirectionAtPoint(valid_points[i - 1]);
        float cos_theta = std::min(1.0f, std::max(-1.0f, static_cast<float>(current_direction.Dot(previous_direction))));
        float theta = std::acos(cos_theta);
        angle_sum += theta;
        deflection_angles.push_back(theta);

        TVector3 current_point(track->LocationAtPoint(valid_points[i]).X(), 
                              track->LocationAtPoint(valid_points[i]).Y(), 
                              track->LocationAtPoint(valid_points[i]).Z());
        TVector3 previous_point(track->LocationAtPoint(valid_points[i - 1]).X(), 
                               track->LocationAtPoint(valid_points[i - 1]).Y(), 
                               track->LocationAtPoint(valid_points[i - 1]).Z());
        separation_sum += (current_point - previous_point).Mag();
    }

    float mean_angle = angle_sum / deflection_angles.size();
    float mean_sep = separation_sum / deflection_angles.size();
    float variance_sum = 0.0f;
    for (const auto& angle : deflection_angles) {
        variance_sum += (angle - mean_angle) * (angle - mean_angle);
    }
    float variance = variance_sum / (deflection_angles.size() - 1);
    const int J = static_cast<int>(deflection_angles.size());
    float sum_theta2 = 0.0f;
    float theta_max  = 0.0f;
    for (const auto& th : deflection_angles) {
        sum_theta2 += th*th;
        theta_max   = std::max(theta_max, th);
    }
    float theta_rms = (J>0) ? std::sqrt(sum_theta2 / J) : std::numeric_limits<float>::lowest();

    mean_deflections.push_back(mean_angle);
    stdev_deflections.push_back(std::sqrt(variance));
    mean_separation.push_back(mean_sep);
    _track_mcs_njoints.push_back((J>0) ? J : std::numeric_limits<int>::lowest());
    _track_deflection_rms.push_back(theta_rms);
    _track_deflection_max.push_back((J>0) ? theta_max : std::numeric_limits<float>::lowest());
}

void TrackAnalysis::calculate_track_deflections_qtrim(
    const art::Ptr<recob::Track>& track,
    std::vector<float>& mean_qtrim_out,
    std::vector<float>& rms_qtrim_out,
    std::vector<float>& q90_out,
    std::vector<int>& njoints_mid_out,
    float alpha, float beta,
    unsigned skip_each_end,
    double qlow, double qhigh,
    float min_seg_cm, float max_seg_cm,
    unsigned min_joints) {
    std::vector<size_t> vtx;
    size_t i0 = track->FirstValidPoint();
    if (i0 == recob::TrackTrajectory::InvalidIndex) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }
    vtx.push_back(i0);
    size_t inext = track->NextValidPoint(i0 + 1);
    while (inext != recob::TrackTrajectory::InvalidIndex) {
        vtx.push_back(inext);
        inext = track->NextValidPoint(inext + 1);
    }
    if (vtx.size() < 3) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }

    std::vector<float> s(vtx.size(), 0.0f);
    for (size_t k=1; k<vtx.size(); ++k) {
        TVector3 a(track->LocationAtPoint(vtx[k-1]).X(),
                   track->LocationAtPoint(vtx[k-1]).Y(),
                   track->LocationAtPoint(vtx[k-1]).Z());
        TVector3 b(track->LocationAtPoint(vtx[k]).X(),
                   track->LocationAtPoint(vtx[k]).Y(),
                   track->LocationAtPoint(vtx[k]).Z());
        s[k] = s[k-1] + (b-a).Mag();
    }
    const float Ltot = s.back();
    if (!(Ltot > 0)) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }

    const float slo = alpha * Ltot;
    const float shi = beta  * Ltot;

    struct Joint { float mid_s; float theta; };
    std::vector<Joint> joints;
    joints.reserve(vtx.size());
    for (size_t k=1; k<vtx.size(); ++k) {
        const auto u = track->DirectionAtPoint(vtx[k-1]).Unit();
        const auto v = track->DirectionAtPoint(vtx[k]).Unit();
        float cosang = std::min(1.0f, std::max(-1.0f, static_cast<float>(u.Dot(v))));
        float theta = std::acos(cosang);

        TVector3 a(track->LocationAtPoint(vtx[k-1]).X(),
                   track->LocationAtPoint(vtx[k-1]).Y(),
                   track->LocationAtPoint(vtx[k-1]).Z());
        TVector3 b(track->LocationAtPoint(vtx[k]).X(),
                   track->LocationAtPoint(vtx[k]).Y(),
                   track->LocationAtPoint(vtx[k]).Z());
        float ell = (b-a).Mag();
        if (ell < min_seg_cm || ell > max_seg_cm) continue;

        float smid = 0.5f * (s[k-1] + s[k]);
        if (smid < slo || smid > shi) continue;

        joints.push_back({smid, theta});
    }

    if (joints.size() < min_joints) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }

    std::sort(joints.begin(), joints.end(),
              [](const Joint& A, const Joint& B){ return A.mid_s < B.mid_s; });
    size_t start = std::min<size_t>(skip_each_end, joints.size());
    size_t stop  = (joints.size() > skip_each_end) ? joints.size() - skip_each_end : start;
    if (stop <= start) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }

    std::vector<float> th;
    th.reserve(stop - start);
    for (size_t k=start; k<stop; ++k) th.push_back(joints[k].theta);

    float q_lo = common::quantile_linear(th, qlow);
    float q_hi = common::quantile_linear(th, qhigh);
    if (!(q_lo <= q_hi)) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }

    double sum=0.0, sum2=0.0; unsigned n=0;
    std::vector<float> kept; kept.reserve(th.size());
    for (float t : th) {
        if (t < q_lo || t > q_hi) continue;
        kept.push_back(t);
        sum  += t;
        sum2 += t*t;
        ++n;
    }
    if (n==0) {
        mean_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        rms_qtrim_out.push_back(std::numeric_limits<float>::lowest());
        q90_out.push_back(std::numeric_limits<float>::lowest());
        njoints_mid_out.push_back(std::numeric_limits<int>::lowest());
        return;
    }

    float mean_qtrim = static_cast<float>(sum / n);
    float rms_qtrim  = static_cast<float>(std::sqrt(sum2 / n));
    float q90        = common::quantile_linear(kept, 0.90);

    mean_qtrim_out.push_back(mean_qtrim);
    rms_qtrim_out.push_back(rms_qtrim);
    q90_out.push_back(q90);
    njoints_mid_out.push_back(static_cast<int>(n));
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)

}

#endif
