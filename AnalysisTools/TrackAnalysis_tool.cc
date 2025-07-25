#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include <iostream>
#include <limits>
#include "AnalysisToolBase.h"
#include "../CommonDefs/Types.h"
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/CalibrationFuncs.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace analysis {

class TrackAnalysis : public AnalysisToolBase {
public:
    TrackAnalysis(const fhicl::ParameterSet& parameter_set);
    ~TrackAnalysis() {}
    void configure(const fhicl::ParameterSet& parameter_set) override;
    void analyseEvent(const art::Event& event, bool is_data) override;
    void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vector, bool is_data, bool is_selected) override;
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    void fill_default();
    float calculate_track_trunk_dedx_by_hits(const std::vector<float>& dedx_values);
    float calculate_track_trunk_dedx_by_range(const std::vector<float>& dedx_values, const std::vector<float>& residual_range_values);
    void calculate_track_deflections(const art::Ptr<recob::Track>& track, std::vector<float>& mean_deflections, 
                                   std::vector<float>& stdev_deflections, std::vector<float>& mean_separation);

    art::InputTag fCALproducer;
    art::InputTag fTRKproducer;
    art::InputTag fHITproducer;
    art::InputTag fCLSproducer;
    art::InputTag fBKTproducer;

    std::vector<float> fADCtoE;
    float fEndSpacepointDistance;

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
    std::vector<int> _track_nhits_u;
    std::vector<int> _track_nhits_v;
    std::vector<int> _track_nhits_y;
    std::vector<float> _track_avg_deflection_mean;
    std::vector<float> _track_avg_deflection_stdev;
    std::vector<float> _track_avg_deflection_separation_mean;
    std::vector<int> _track_end_spacepoints;
};

TrackAnalysis::TrackAnalysis(const fhicl::ParameterSet& parameter_set) {
    fHITproducer = parameter_set.get<art::InputTag>("HITproducer", "gaushit");
    fCLSproducer = parameter_set.get<art::InputTag>("CLSproducer", "pandora");
    fTRKproducer = parameter_set.get<art::InputTag>("TRKproducer", "pandoraTrack");
    fCALproducer = parameter_set.get<art::InputTag>("CALproducer", "pandoraTrackcali");
    fBKTproducer = parameter_set.get<art::InputTag>("BKTproducer", "gaushitTruthMatch");

    fADCtoE = parameter_set.get<std::vector<float>>("ADCtoE");
    fEndSpacepointDistance = parameter_set.get<float>("EndSpacepointDistance", 5.0);
}

void TrackAnalysis::configure(const fhicl::ParameterSet& parameter_set) {}

void TrackAnalysis::analyseEvent(const art::Event& event, bool is_data) {
    _evt = event.event();
    _sbr = event.subRun();
    _run = event.run();
}

void TrackAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vector, 
                                 bool is_data, bool is_selected) {
    const auto& calorimetry_proxy = proxy::getCollection<std::vector<recob::Track>>(
        event, fTRKproducer, proxy::withAssociated<anab::Calorimetry>(fCALproducer));

    TVector3 neutrino_vertex;
    for (const auto& pfp : slice_pfp_vector) {
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

    for (size_t pfp_index = 0; pfp_index < slice_pfp_vector.size(); ++pfp_index) {
        const auto& pfp = slice_pfp_vector[pfp_index];
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

            auto calorimetry_objects = calorimetry_proxy[track.key()].get<anab::Calorimetry>();
            for (const auto& calo : calorimetry_objects) {
                int plane = calo->PlaneID().Plane;
                if (plane < 0 || plane > 2) continue;

                const auto& dqdx_values = calo->dQdx();
                const auto& residual_range = calo->ResidualRange();
                const auto& pitch_values = calo->TrkPitchVec();
                const auto& xyz_values = calo->XYZ();

                float calo_energy = 0.0f;
                std::vector<float> corrected_dedx_values;
                for (size_t i = 0; i < dqdx_values.size(); ++i) {
                    float dedx = common::ModBoxCorrection(dqdx_values[i] * fADCtoE[plane], 
                                                         xyz_values[i].X(), xyz_values[i].Y(), xyz_values[i].Z());
                    corrected_dedx_values.push_back(dedx);
                    calo_energy += dedx * pitch_values[i];
                }

                int hit_count = dqdx_values.size();
                float trunk_dedx = calculate_track_trunk_dedx_by_hits(corrected_dedx_values);
                float trunk_rr_dedx = calculate_track_trunk_dedx_by_range(corrected_dedx_values, residual_range);

                if (plane == 0) {
                    _track_calo_energy_u.push_back(calo_energy);
                    _track_nhits_u.push_back(hit_count);
                    _track_trunk_dedx_u.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_u.push_back(trunk_rr_dedx);
                } else if (plane == 1) {
                    _track_calo_energy_v.push_back(calo_energy);
                    _track_nhits_v.push_back(hit_count);
                    _track_trunk_dedx_v.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_v.push_back(trunk_rr_dedx);
                } else if (plane == 2) {
                    _track_calo_energy_y.push_back(calo_energy);
                    _track_nhits_y.push_back(hit_count);
                    _track_trunk_dedx_y.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_y.push_back(trunk_rr_dedx);
                }
            }

            calculate_track_deflections(track, _track_avg_deflection_mean, _track_avg_deflection_stdev, 
                                       _track_avg_deflection_separation_mean);

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
    _track_start_x.push_back(std::numeric_limits<float>::lowest());
    _track_start_y.push_back(std::numeric_limits<float>::lowest());
    _track_start_z.push_back(std::numeric_limits<float>::lowest());
    _track_sce_start_x.push_back(std::numeric_limits<float>::lowest());
    _track_sce_start_y.push_back(std::numeric_limits<float>::lowest());
    _track_sce_start_z.push_back(std::numeric_limits<float>::lowest());
    _track_end_x.push_back(std::numeric_limits<float>::lowest());
    _track_end_y.push_back(std::numeric_limits<float>::lowest());
    _track_end_z.push_back(std::numeric_limits<float>::lowest());
    _track_sce_end_x.push_back(std::numeric_limits<float>::lowest());
    _track_sce_end_y.push_back(std::numeric_limits<float>::lowest());
    _track_sce_end_z.push_back(std::numeric_limits<float>::lowest());
    _track_direction_x.push_back(std::numeric_limits<float>::lowest());
    _track_direction_y.push_back(std::numeric_limits<float>::lowest());
    _track_direction_z.push_back(std::numeric_limits<float>::lowest());
    _track_distance_to_vertex.push_back(std::numeric_limits<float>::lowest());
    _track_theta.push_back(std::numeric_limits<float>::lowest());
    _track_phi.push_back(std::numeric_limits<float>::lowest());
    _track_length.push_back(std::numeric_limits<float>::lowest());
    _track_calo_energy_u.push_back(std::numeric_limits<float>::lowest());
    _track_calo_energy_v.push_back(std::numeric_limits<float>::lowest());
    _track_calo_energy_y.push_back(std::numeric_limits<float>::lowest());
    _track_trunk_dedx_u.push_back(std::numeric_limits<float>::lowest());
    _track_trunk_dedx_v.push_back(std::numeric_limits<float>::lowest());
    _track_trunk_dedx_y.push_back(std::numeric_limits<float>::lowest());
    _track_trunk_rr_dedx_u.push_back(std::numeric_limits<float>::lowest());
    _track_trunk_rr_dedx_v.push_back(std::numeric_limits<float>::lowest());
    _track_trunk_rr_dedx_y.push_back(std::numeric_limits<float>::lowest());
    _track_nhits_u.push_back(std::numeric_limits<int>::lowest());
    _track_nhits_v.push_back(std::numeric_limits<int>::lowest());
    _track_nhits_y.push_back(std::numeric_limits<int>::lowest());
    _track_avg_deflection_mean.push_back(std::numeric_limits<float>::lowest());
    _track_avg_deflection_stdev.push_back(std::numeric_limits<float>::lowest());
    _track_avg_deflection_separation_mean.push_back(std::numeric_limits<float>::lowest());
    _track_end_spacepoints.push_back(std::numeric_limits<int>::lowest());
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
    tree->Branch("track_nhits_u", "std::vector<int>", &_track_nhits_u);
    tree->Branch("track_nhits_v", "std::vector<int>", &_track_nhits_v);
    tree->Branch("track_nhits_y", "std::vector<int>", &_track_nhits_y);
    tree->Branch("track_avg_deflection_mean", "std::vector<float>", &_track_avg_deflection_mean);
    tree->Branch("track_avg_deflection_stdev", "std::vector<float>", &_track_avg_deflection_stdev);
    tree->Branch("track_avg_deflection_separation_mean", "std::vector<float>", &_track_avg_deflection_separation_mean);
    tree->Branch("track_end_spacepoints", "std::vector<int>", &_track_end_spacepoints);
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
    _track_nhits_u.clear();
    _track_nhits_v.clear();
    _track_nhits_y.clear();
    _track_avg_deflection_mean.clear();
    _track_avg_deflection_stdev.clear();
    _track_avg_deflection_separation_mean.clear();
    _track_end_spacepoints.clear();
}

float TrackAnalysis::calculate_track_trunk_dedx_by_hits(const std::vector<float>& dedx_values) {
    size_t hit_count = dedx_values.size();
    int first_hit_index = static_cast<int>(hit_count) - 4; 
    int last_hit_index = static_cast<int>(hit_count - (hit_count / 3)) - 1;

    if (first_hit_index - last_hit_index < 5) {
        return std::numeric_limits<float>::lowest();
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
        return std::numeric_limits<float>::lowest();
    }

    std::vector<std::pair<float, unsigned int>> residual_range_indices;
    float max_residual_range = std::numeric_limits<float>::lowest();
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
        return std::numeric_limits<float>::lowest();
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

    return (truncated_hit_count == 0) ? std::numeric_limits<float>::lowest() 
                                   : truncated_total / static_cast<float>(truncated_hit_count);
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
        mean_deflections.push_back(std::numeric_limits<float>::lowest());
        stdev_deflections.push_back(std::numeric_limits<float>::lowest());
        mean_separation.push_back(std::numeric_limits<float>::lowest());
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

    mean_deflections.push_back(mean_angle);
    stdev_deflections.push_back(std::sqrt(variance));
    mean_separation.push_back(mean_sep);
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)

}

#endif