#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "Common/EnergyCalibration.h"
#include "Common/GeometryUtils.h"
#include "Common/ProxyTypes.h"
#include "Common/SpaceChargeCorrections.h"
#include "Common/Quantile.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <limits>
#include "TVector3.h"

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
    float calculate_track_trunk_rr_qtrim(const std::vector<float>& dedx,
                                         const std::vector<float>& rr,
                                         float alpha=0.25f, float beta=0.75f,
                                         unsigned skip_each_end=3,
                                         double qlow=0.10, double qhigh=0.80,
                                         unsigned min_hits_in_band=8);
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
    std::vector<float> _track_trunk_dedx_u;
    std::vector<float> _track_trunk_dedx_v;
    std::vector<float> _track_trunk_dedx_y;
    std::vector<float> _track_trunk_rr_dedx_u;
    std::vector<float> _track_trunk_rr_dedx_v;
    std::vector<float> _track_trunk_rr_dedx_y;
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

            auto calorimetry_objects = calorimetry_proxy[track.key()].get<anab::Calorimetry>();

            for (const auto& calo : calorimetry_objects) {
                int plane = calo->PlaneID().Plane;
                if (plane < 0 || plane > 2) continue;

                const auto& dqdx_values = calo->dQdx();
                const auto& residual_range = calo->ResidualRange();
                const auto& pitch_values = calo->TrkPitchVec();
                const auto& xyz_values = calo->XYZ();

                float calo_energy = 0.0f;
                std::vector<float> dedx_values;
                dedx_values.reserve(dqdx_values.size());
                for (size_t i = 0; i < dqdx_values.size(); ++i) {
                    float dedx = common::ModBoxCorrection(dqdx_values[i] * fADCtoE[plane],
                                                         xyz_values[i].X(), xyz_values[i].Y(), xyz_values[i].Z());
                    dedx_values.push_back(dedx);
                    calo_energy += dedx * pitch_values[i];
                }

                float trunk_dedx = calculate_track_trunk_dedx_by_hits(dedx_values);
                float trunk_rr_dedx = calculate_track_trunk_rr_qtrim(
                    dedx_values, residual_range,
                    0.25f, 0.75f, 3, 0.10, 0.80, 8);

                if (plane == 0) {
                    _track_calo_energy_u.push_back(calo_energy);
                    _track_trunk_dedx_u.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_u.push_back(trunk_rr_dedx);
                } else if (plane == 1) {
                    _track_calo_energy_v.push_back(calo_energy);
                    _track_trunk_dedx_v.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_v.push_back(trunk_rr_dedx);
                } else if (plane == 2) {
                    _track_calo_energy_y.push_back(calo_energy);
                    _track_trunk_dedx_y.push_back(trunk_dedx);
                    _track_trunk_rr_dedx_y.push_back(trunk_rr_dedx);
                }
            }
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
DEFINE_ART_CLASS_TOOL(TrackAnalysis)

}

#endif
