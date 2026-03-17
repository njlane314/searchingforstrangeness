#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "Common/EnergyCalibration.h"
#include "Common/GeometryUtils.h"
#include "Common/ProxyTypes.h"
#include "Common/SpaceChargeCorrections.h"
#include "Common/StatisticsQuantile.h"

#include "TVector3.h"

#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include <utility>

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
    bool compute_selected_median_iqr(std::vector<float> values,
                                     float& median,
                                     float& iqr,
                                     unsigned min_hits);

    bool compute_vertex_trunk_mip_stats(const std::vector<float>& dedx,
                                        const std::vector<float>& rr,
                                        bool vertex_at_track_start,
                                        float& median,
                                        float& iqr,
                                        float d_lo = 5.0f,
                                        float d_hi = 20.0f,
                                        unsigned min_hits = 6);

    bool compute_bulk_mip_stats(const std::vector<float>& dedx,
                                const std::vector<float>& rr,
                                float& median,
                                float& iqr,
                                float rr_lo = 0.20f,
                                float rr_hi = 0.80f,
                                unsigned min_hits = 8);

    float compute_bragg_ratio(const std::vector<float>& dedx,
                              const std::vector<float>& rr,
                              float end_rr_frac = 0.10f,
                              float bulk_rr_lo = 0.20f,
                              float bulk_rr_hi = 0.80f,
                              unsigned min_hits_end = 4,
                              unsigned min_hits_bulk = 8);
    art::InputTag fTRKproducer;
    art::InputTag fCALproducer;

    std::vector<float> fADCtoE;

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

    std::vector<float> _track_vertex_trunk_mip_median_u;
    std::vector<float> _track_vertex_trunk_mip_median_v;
    std::vector<float> _track_vertex_trunk_mip_median_y;
    std::vector<float> _track_vertex_trunk_mip_iqr_u;
    std::vector<float> _track_vertex_trunk_mip_iqr_v;
    std::vector<float> _track_vertex_trunk_mip_iqr_y;

    std::vector<float> _track_bulk_mip_median_u;
    std::vector<float> _track_bulk_mip_median_v;
    std::vector<float> _track_bulk_mip_median_y;
    std::vector<float> _track_bulk_mip_iqr_u;
    std::vector<float> _track_bulk_mip_iqr_v;
    std::vector<float> _track_bulk_mip_iqr_y;

    std::vector<float> _track_bragg_ratio_u;
    std::vector<float> _track_bragg_ratio_v;
    std::vector<float> _track_bragg_ratio_y;
};

TrackAnalysis::TrackAnalysis(const fhicl::ParameterSet& parameter_set) {
    fTRKproducer = parameter_set.get<art::InputTag>("TRKproducer", "pandoraTrack");
    fCALproducer = parameter_set.get<art::InputTag>("CALproducer", "pandoraTrackcali");

    fADCtoE = parameter_set.get<std::vector<float>>("ADCtoE");
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

    TVector3 neutrino_vertex(0.0, 0.0, 0.0);
    bool have_neutrino_vertex = false;
    for (const auto& pfp : slice_pfp_vec) {
        if (pfp->IsPrimary()) {
            double xyz[3] = {};
            auto vertices = pfp.get<recob::Vertex>();
            if (vertices.size() != 1) {
                continue;
            } else {
                vertices.at(0)->XYZ(xyz);
                neutrino_vertex.SetXYZ(xyz[0], xyz[1], xyz[2]);
                have_neutrino_vertex = true;
            }
            break;
        }
    }

    const float nan = std::numeric_limits<float>::quiet_NaN();

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

            bool vertex_at_track_start = false;
            if (have_neutrino_vertex) {
                TVector3 start_to_vertex(track->Start().X(), track->Start().Y(), track->Start().Z());
                start_to_vertex -= neutrino_vertex;
                TVector3 end_to_vertex(track->End().X(), track->End().Y(), track->End().Z());
                end_to_vertex -= neutrino_vertex;

                const float d_start = start_to_vertex.Mag();
                const float d_end = end_to_vertex.Mag();
                vertex_at_track_start = (d_start <= d_end);
                _track_distance_to_vertex.push_back(d_start);
            } else {
                _track_distance_to_vertex.push_back(nan);
            }

            std::array<float, 3> calo_energy;
            std::array<float, 3> vertex_trunk_median;
            std::array<float, 3> vertex_trunk_iqr;
            std::array<float, 3> bulk_median;
            std::array<float, 3> bulk_iqr;
            std::array<float, 3> bragg_ratio;

            calo_energy.fill(nan);
            vertex_trunk_median.fill(nan);
            vertex_trunk_iqr.fill(nan);
            bulk_median.fill(nan);
            bulk_iqr.fill(nan);
            bragg_ratio.fill(nan);

            auto calorimetry_objects = calorimetry_proxy[track.key()].get<anab::Calorimetry>();

            for (const auto& calo : calorimetry_objects) {
                const int plane = calo->PlaneID().Plane;
                if (plane < 0 || plane > 2) continue;

                const auto& dqdx_values = calo->dQdx();
                const auto& residual_range = calo->ResidualRange();
                const auto& pitch_values = calo->TrkPitchVec();
                const auto& xyz_values = calo->XYZ();

                const size_t n_hits = std::min({
                    dqdx_values.size(),
                    residual_range.size(),
                    pitch_values.size(),
                    xyz_values.size()
                });

                if (n_hits == 0) continue;

                float plane_calo_energy = 0.0f;
                std::vector<float> dedx_values;
                dedx_values.reserve(n_hits);

                std::vector<float> rr_values;
                rr_values.reserve(n_hits);

                for (size_t i = 0; i < n_hits; ++i) {
                    const float dedx = common::ModBoxCorrection(dqdx_values[i] * fADCtoE[plane],
                                                                xyz_values[i].X(), xyz_values[i].Y(), xyz_values[i].Z());
                    dedx_values.push_back(dedx);
                    rr_values.push_back(residual_range[i]);

                    if (std::isfinite(dedx) && std::isfinite(pitch_values[i])) {
                        plane_calo_energy += dedx * pitch_values[i];
                    }
                }

                calo_energy[plane] = plane_calo_energy;

                if (have_neutrino_vertex) {
                    compute_vertex_trunk_mip_stats(dedx_values, rr_values,
                                                   vertex_at_track_start,
                                                   vertex_trunk_median[plane],
                                                   vertex_trunk_iqr[plane],
                                                   5.0f, 20.0f, 6);
                }

                compute_bulk_mip_stats(dedx_values, rr_values,
                                       bulk_median[plane],
                                       bulk_iqr[plane],
                                       0.20f, 0.80f, 8);

                bragg_ratio[plane] = compute_bragg_ratio(dedx_values, rr_values,
                                                         0.10f, 0.20f, 0.80f, 4, 8);
            }

            _track_calo_energy_u.push_back(calo_energy[0]);
            _track_calo_energy_v.push_back(calo_energy[1]);
            _track_calo_energy_y.push_back(calo_energy[2]);

            _track_vertex_trunk_mip_median_u.push_back(vertex_trunk_median[0]);
            _track_vertex_trunk_mip_median_v.push_back(vertex_trunk_median[1]);
            _track_vertex_trunk_mip_median_y.push_back(vertex_trunk_median[2]);
            _track_vertex_trunk_mip_iqr_u.push_back(vertex_trunk_iqr[0]);
            _track_vertex_trunk_mip_iqr_v.push_back(vertex_trunk_iqr[1]);
            _track_vertex_trunk_mip_iqr_y.push_back(vertex_trunk_iqr[2]);

            _track_bulk_mip_median_u.push_back(bulk_median[0]);
            _track_bulk_mip_median_v.push_back(bulk_median[1]);
            _track_bulk_mip_median_y.push_back(bulk_median[2]);
            _track_bulk_mip_iqr_u.push_back(bulk_iqr[0]);
            _track_bulk_mip_iqr_v.push_back(bulk_iqr[1]);
            _track_bulk_mip_iqr_y.push_back(bulk_iqr[2]);

            _track_bragg_ratio_u.push_back(bragg_ratio[0]);
            _track_bragg_ratio_v.push_back(bragg_ratio[1]);
            _track_bragg_ratio_y.push_back(bragg_ratio[2]);
        } else {
            fill_default();
        }
    }
}

void TrackAnalysis::fill_default() {
    const float nan = std::numeric_limits<float>::quiet_NaN();

    _track_pfp_ids.push_back(std::numeric_limits<size_t>::max());
    _track_start_x.push_back(nan);
    _track_start_y.push_back(nan);
    _track_start_z.push_back(nan);
    _track_sce_start_x.push_back(nan);
    _track_sce_start_y.push_back(nan);
    _track_sce_start_z.push_back(nan);
    _track_end_x.push_back(nan);
    _track_end_y.push_back(nan);
    _track_end_z.push_back(nan);
    _track_sce_end_x.push_back(nan);
    _track_sce_end_y.push_back(nan);
    _track_sce_end_z.push_back(nan);
    _track_direction_x.push_back(nan);
    _track_direction_y.push_back(nan);
    _track_direction_z.push_back(nan);
    _track_distance_to_vertex.push_back(nan);
    _track_theta.push_back(nan);
    _track_phi.push_back(nan);
    _track_length.push_back(nan);

    _track_calo_energy_u.push_back(nan);
    _track_calo_energy_v.push_back(nan);
    _track_calo_energy_y.push_back(nan);

    _track_vertex_trunk_mip_median_u.push_back(nan);
    _track_vertex_trunk_mip_median_v.push_back(nan);
    _track_vertex_trunk_mip_median_y.push_back(nan);
    _track_vertex_trunk_mip_iqr_u.push_back(nan);
    _track_vertex_trunk_mip_iqr_v.push_back(nan);
    _track_vertex_trunk_mip_iqr_y.push_back(nan);

    _track_bulk_mip_median_u.push_back(nan);
    _track_bulk_mip_median_v.push_back(nan);
    _track_bulk_mip_median_y.push_back(nan);
    _track_bulk_mip_iqr_u.push_back(nan);
    _track_bulk_mip_iqr_v.push_back(nan);
    _track_bulk_mip_iqr_y.push_back(nan);

    _track_bragg_ratio_u.push_back(nan);
    _track_bragg_ratio_v.push_back(nan);
    _track_bragg_ratio_y.push_back(nan);
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

    tree->Branch("track_vertex_trunk_mip_median_u", "std::vector<float>", &_track_vertex_trunk_mip_median_u);
    tree->Branch("track_vertex_trunk_mip_median_v", "std::vector<float>", &_track_vertex_trunk_mip_median_v);
    tree->Branch("track_vertex_trunk_mip_median_y", "std::vector<float>", &_track_vertex_trunk_mip_median_y);
    tree->Branch("track_vertex_trunk_mip_iqr_u", "std::vector<float>", &_track_vertex_trunk_mip_iqr_u);
    tree->Branch("track_vertex_trunk_mip_iqr_v", "std::vector<float>", &_track_vertex_trunk_mip_iqr_v);
    tree->Branch("track_vertex_trunk_mip_iqr_y", "std::vector<float>", &_track_vertex_trunk_mip_iqr_y);

    tree->Branch("track_bulk_mip_median_u", "std::vector<float>", &_track_bulk_mip_median_u);
    tree->Branch("track_bulk_mip_median_v", "std::vector<float>", &_track_bulk_mip_median_v);
    tree->Branch("track_bulk_mip_median_y", "std::vector<float>", &_track_bulk_mip_median_y);
    tree->Branch("track_bulk_mip_iqr_u", "std::vector<float>", &_track_bulk_mip_iqr_u);
    tree->Branch("track_bulk_mip_iqr_v", "std::vector<float>", &_track_bulk_mip_iqr_v);
    tree->Branch("track_bulk_mip_iqr_y", "std::vector<float>", &_track_bulk_mip_iqr_y);

    tree->Branch("track_bragg_ratio_u", "std::vector<float>", &_track_bragg_ratio_u);
    tree->Branch("track_bragg_ratio_v", "std::vector<float>", &_track_bragg_ratio_v);
    tree->Branch("track_bragg_ratio_y", "std::vector<float>", &_track_bragg_ratio_y);
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

    _track_vertex_trunk_mip_median_u.clear();
    _track_vertex_trunk_mip_median_v.clear();
    _track_vertex_trunk_mip_median_y.clear();
    _track_vertex_trunk_mip_iqr_u.clear();
    _track_vertex_trunk_mip_iqr_v.clear();
    _track_vertex_trunk_mip_iqr_y.clear();

    _track_bulk_mip_median_u.clear();
    _track_bulk_mip_median_v.clear();
    _track_bulk_mip_median_y.clear();
    _track_bulk_mip_iqr_u.clear();
    _track_bulk_mip_iqr_v.clear();
    _track_bulk_mip_iqr_y.clear();

    _track_bragg_ratio_u.clear();
    _track_bragg_ratio_v.clear();
    _track_bragg_ratio_y.clear();
}

bool TrackAnalysis::compute_selected_median_iqr(std::vector<float> values,
                                                float& median,
                                                float& iqr,
                                                unsigned min_hits) {
    median = std::numeric_limits<float>::quiet_NaN();
    iqr = std::numeric_limits<float>::quiet_NaN();

    if (values.size() < min_hits) {
        return false;
    }

    std::sort(values.begin(), values.end());

    const float q25 = common::quantile_linear(values, 0.25);
    median = common::quantile_linear(values, 0.50);
    const float q75 = common::quantile_linear(values, 0.75);
    iqr = q75 - q25;

    return std::isfinite(median) && std::isfinite(iqr);
}

bool TrackAnalysis::compute_vertex_trunk_mip_stats(const std::vector<float>& dedx,
                                                   const std::vector<float>& rr,
                                                   bool vertex_at_track_start,
                                                   float& median,
                                                   float& iqr,
                                                   float d_lo,
                                                   float d_hi,
                                                   unsigned min_hits) {
    median = std::numeric_limits<float>::quiet_NaN();
    iqr = std::numeric_limits<float>::quiet_NaN();

    if (dedx.size() != rr.size() || dedx.empty()) {
        return false;
    }

    float rrmax = -1.0f;
    for (float r : rr) {
        if (std::isfinite(r)) {
            rrmax = std::max(rrmax, r);
        }
    }
    if (!(rrmax > 0.0f)) {
        return false;
    }

    std::vector<float> values;
    values.reserve(dedx.size());

    for (size_t i = 0; i < dedx.size(); ++i) {
        const float v = dedx[i];
        const float r = rr[i];

        if (!std::isfinite(v) || v <= 0.0f) continue;
        if (!std::isfinite(r) || r < 0.0f) continue;

        const float distance_from_vertex = vertex_at_track_start ? (rrmax - r) : r;
        if (distance_from_vertex < d_lo || distance_from_vertex > d_hi) continue;

        values.push_back(v);
    }

    return compute_selected_median_iqr(std::move(values), median, iqr, min_hits);
}

bool TrackAnalysis::compute_bulk_mip_stats(const std::vector<float>& dedx,
                                           const std::vector<float>& rr,
                                           float& median,
                                           float& iqr,
                                           float rr_lo,
                                           float rr_hi,
                                           unsigned min_hits) {
    median = std::numeric_limits<float>::quiet_NaN();
    iqr = std::numeric_limits<float>::quiet_NaN();

    if (dedx.size() != rr.size() || dedx.empty()) {
        return false;
    }

    float rrmax = -1.0f;
    for (float r : rr) {
        if (std::isfinite(r)) {
            rrmax = std::max(rrmax, r);
        }
    }
    if (!(rrmax > 0.0f)) {
        return false;
    }

    std::vector<float> values;
    values.reserve(dedx.size());

    for (size_t i = 0; i < dedx.size(); ++i) {
        const float v = dedx[i];
        const float r = rr[i];

        if (!std::isfinite(v) || v <= 0.0f) continue;
        if (!std::isfinite(r) || r < 0.0f) continue;

        const float u = r / rrmax;
        if (u < rr_lo || u > rr_hi) continue;

        values.push_back(v);
    }

    return compute_selected_median_iqr(std::move(values), median, iqr, min_hits);
}

float TrackAnalysis::compute_bragg_ratio(const std::vector<float>& dedx,
                                         const std::vector<float>& rr,
                                         float end_rr_frac,
                                         float bulk_rr_lo,
                                         float bulk_rr_hi,
                                         unsigned min_hits_end,
                                         unsigned min_hits_bulk) {
    if (dedx.size() != rr.size() || dedx.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    float rrmax = -1.0f;
    for (float r : rr) {
        if (std::isfinite(r)) {
            rrmax = std::max(rrmax, r);
        }
    }
    if (!(rrmax > 0.0f)) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::vector<float> end_values;
    std::vector<float> bulk_values;
    end_values.reserve(dedx.size());
    bulk_values.reserve(dedx.size());

    for (size_t i = 0; i < dedx.size(); ++i) {
        const float v = dedx[i];
        const float r = rr[i];

        if (!std::isfinite(v) || v <= 0.0f) continue;
        if (!std::isfinite(r) || r < 0.0f) continue;

        const float u = r / rrmax;
        if (u <= end_rr_frac) {
            end_values.push_back(v);
        }
        if (u >= bulk_rr_lo && u <= bulk_rr_hi) {
            bulk_values.push_back(v);
        }
    }

    float end_median = std::numeric_limits<float>::quiet_NaN();
    float end_iqr = std::numeric_limits<float>::quiet_NaN();
    float bulk_median = std::numeric_limits<float>::quiet_NaN();
    float bulk_iqr = std::numeric_limits<float>::quiet_NaN();

    if (!compute_selected_median_iqr(std::move(end_values), end_median, end_iqr, min_hits_end)) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    if (!compute_selected_median_iqr(std::move(bulk_values), bulk_median, bulk_iqr, min_hits_bulk)) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    if (!(bulk_median > 0.0f)) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    return end_median / bulk_median;
}

DEFINE_ART_CLASS_TOOL(TrackAnalysis)

}

#endif
