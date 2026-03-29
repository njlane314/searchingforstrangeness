#ifndef ANALYSIS_COMPONENTGRAPHANALYSIS_CXX
#define ANALYSIS_COMPONENTGRAPHANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "Common/ProxyTypes.h"
#include "Common/GraphEngine.h"

#include "TTree.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace analysis {

class GraphAnalysis : public AnalysisToolBase {
public:
    explicit GraphAnalysis(const fhicl::ParameterSet& parameter_set);
    ~GraphAnalysis() override = default;

    void configure(const fhicl::ParameterSet& parameter_set) override;
    void analyseEvent(const art::Event& event, bool is_data) override;
    void analyseSlice(const art::Event& event,
                      std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                      bool is_data,
                      bool is_selected) override;
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    struct SlicePV {
        bool valid = false;
        recob::Vertex::Point_t pos;
    };

    struct PlaneKey {
        int cryo = 0;
        int tpc = 0;
        int plane = 0;

        bool operator<(const PlaneKey& other) const {
            if (cryo != other.cryo) return cryo < other.cryo;
            if (tpc != other.tpc) return tpc < other.tpc;
            return plane < other.plane;
        }
    };

    void load_bad_channels(const std::string& filename);
    void build_bad_wire_ranges(const geo::GeometryCore* geom);
    void read_config(const fhicl::ParameterSet& parameter_set);

    cg::PlaneGeometry make_plane_geometry(const geo::PlaneID& plane_id,
                                          const detinfo::DetectorProperties* det_prop,
                                          const geo::GeometryCore* geom) const;

    void collect_slice_hits(const art::Event& event,
                            int slice_index,
                            std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                            const detinfo::DetectorProperties* det_prop,
                            const geo::GeometryCore* geom,
                            std::map<PlaneKey, std::vector<cg::HitInput>>& hits_by_plane) const;

    void write_plane_replay(int plane_uid,
                            const cg::PlaneInput& input);

    void write_nominal_output(int plane_uid,
                              const cg::Result& result,
                              bool have_plane_pv,
                              float pv_wire_cm,
                              float pv_drift_cm);

    bool find_slice_primary_vertex(std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                                   SlicePV& out) const;

    bool project_pv_to_plane(recob::Vertex::Point_t const& pos,
                             geo::PlaneID const& plane_id,
                             geo::GeometryCore const* geom,
                             float& wire_cm,
                             float& drift_cm) const;

    art::InputTag fHITproducer;
    art::InputTag fCLSproducer;
    bool fOnlySelectedSlices = true;
    std::string fBadChannelFile;
    std::set<unsigned int> fBadChannels;
    std::map<PlaneKey, std::vector<std::pair<int, int>>> fBadWireRanges;
    bool fBadWireRangesBuilt = false;

    cg::Params fNominalParams;
    cg::GraphEngine fEngine;

    int _run = -1;
    int _sbr = -1;
    int _evt = -1;
    int _da_slice_counter = 0;
    int _da_next_plane_uid = 0;
    int _da_next_component_id = 0;

    // immutable replay state per plane
    std::vector<int> _da_plane_uid;
    std::vector<int> _da_plane_slice_index;
    std::vector<int> _da_plane_cryo;
    std::vector<int> _da_plane_tpc;
    std::vector<int> _da_plane_plane;
    std::vector<float> _da_plane_wire_pitch_cm;
    std::vector<float> _da_plane_drift_cm_per_tick;
    std::vector<float> _da_plane_wire_min_cm;
    std::vector<float> _da_plane_wire_max_cm;
    std::vector<float> _da_plane_drift_min_cm;
    std::vector<float> _da_plane_drift_max_cm;
    std::vector<int> _da_plane_has_primary_vtx;
    std::vector<float> _da_plane_primary_vtx_wire_cm;
    std::vector<float> _da_plane_primary_vtx_drift_cm;

    // dead-region replay overlay: one drawable rectangle per masked wire band
    std::vector<int> _da_dead_plane_uid;
    std::vector<int> _da_dead_slice_index;
    std::vector<int> _da_dead_cryo;
    std::vector<int> _da_dead_tpc;
    std::vector<int> _da_dead_plane;
    std::vector<float> _da_dead_wire_min_cm;
    std::vector<float> _da_dead_wire_max_cm;
    std::vector<float> _da_dead_drift_min_cm;
    std::vector<float> _da_dead_drift_max_cm;

    // hit branches (raw replay inputs + nominal local geometry / labels)
    std::vector<int> _da_hit_plane_uid;
    std::vector<int> _da_hit_slice_index;
    std::vector<std::size_t> _da_hit_key;
    std::vector<int> _da_hit_cryo;
    std::vector<int> _da_hit_tpc;
    std::vector<int> _da_hit_plane;
    std::vector<int> _da_hit_wire;
    std::vector<float> _da_hit_tick;
    std::vector<float> _da_hit_wire_cm;
    std::vector<float> _da_hit_drift_cm;
    std::vector<float> _da_hit_charge;
    std::vector<float> _da_hit_peakamp;
    std::vector<float> _da_hit_rms;
    std::vector<float> _da_hit_dir_wire;
    std::vector<float> _da_hit_dir_drift;
    std::vector<float> _da_hit_anis;
    std::vector<float> _da_hit_scale;
    std::vector<int> _da_hit_component_id;
    std::vector<int> _da_hit_activity_id;
    std::vector<int> _da_hit_is_boundary;
    std::vector<int> _da_hit_is_start;
    std::vector<int> _da_hit_label;

    // component branches (nominal interpretation)
    std::vector<int> _da_comp_plane_uid;
    std::vector<int> _da_comp_slice_index;
    std::vector<int> _da_comp_id;
    std::vector<int> _da_comp_activity_id;
    std::vector<int> _da_comp_cryo;
    std::vector<int> _da_comp_tpc;
    std::vector<int> _da_comp_plane;
    std::vector<int> _da_comp_n_hits;
    std::vector<int> _da_comp_label;
    std::vector<int> _da_comp_parent_id;
    std::vector<float> _da_comp_qsum;
    std::vector<float> _da_comp_start_wire;
    std::vector<float> _da_comp_start_tick;
    std::vector<float> _da_comp_start_wire_cm;
    std::vector<float> _da_comp_start_drift_cm;
    std::vector<float> _da_comp_dir_wire;
    std::vector<float> _da_comp_dir_drift;
    std::vector<float> _da_comp_start_asym;
    std::vector<float> _da_comp_bbox_min_wire_cm;
    std::vector<float> _da_comp_bbox_max_wire_cm;
    std::vector<float> _da_comp_bbox_min_drift_cm;
    std::vector<float> _da_comp_bbox_max_drift_cm;
    std::vector<float> _da_comp_gap_cm;
    std::vector<float> _da_comp_emptiness;
    std::vector<float> _da_comp_visibility;
    std::vector<float> _da_comp_pointing;
    std::vector<float> _da_comp_path_break;
    std::vector<float> _da_comp_score;
    std::vector<float> _da_comp_activity_vtx_wire_cm;
    std::vector<float> _da_comp_activity_vtx_drift_cm;
    std::vector<float> _da_comp_activity_vtx_residual;
    std::vector<float> _da_comp_pv_to_activity_dist2d_cm;

    // edge branches (nominal interpretation)
    std::vector<int> _da_edge_plane_uid;
    std::vector<int> _da_edge_slice_index;
    std::vector<int> _da_edge_plane;
    std::vector<int> _da_edge_src_component;
    std::vector<int> _da_edge_dst_component;
    std::vector<int> _da_edge_kind;
    std::vector<float> _da_edge_src_wire_cm;
    std::vector<float> _da_edge_src_drift_cm;
    std::vector<float> _da_edge_dst_wire_cm;
    std::vector<float> _da_edge_dst_drift_cm;
    std::vector<float> _da_edge_gap_cm;
    std::vector<float> _da_edge_emptiness;
    std::vector<float> _da_edge_visibility;
    std::vector<float> _da_edge_pointing;
    std::vector<float> _da_edge_path_break;
    std::vector<float> _da_edge_score;
};

GraphAnalysis::GraphAnalysis(const fhicl::ParameterSet& parameter_set)
    : fEngine(cg::Params {}) {
    this->read_config(parameter_set);
}

void GraphAnalysis::load_bad_channels(const std::string& filename) {
    fBadChannels.clear();

    std::ifstream in(filename);
    if (!in.is_open()) {
        throw art::Exception(art::errors::Configuration)
            << "Cannot open bad channel file: " << filename;
    }

    std::string line;
    while (std::getline(in, line)) {
        const auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line.erase(comment_pos);
        }

        std::stringstream ss(line);
        unsigned first = 0u;
        unsigned second = 0u;
        if (!(ss >> first)) continue;

        if (ss >> second) {
            if (second < first) std::swap(first, second);
            for (unsigned ch = first; ch <= second; ++ch) {
                fBadChannels.insert(ch);
            }
        } else {
            fBadChannels.insert(first);
        }
    }
}

void GraphAnalysis::build_bad_wire_ranges(const geo::GeometryCore* geom) {
    fBadWireRanges.clear();

    if (!geom || fBadChannels.empty()) {
        fBadWireRangesBuilt = true;
        return;
    }

    std::map<PlaneKey, std::vector<int>> bad_wires_by_plane;
    for (const unsigned int channel : fBadChannels) {
        if (channel >= static_cast<unsigned int>(geom->Nchannels())) continue;

        const auto wire_ids = geom->ChannelToWire(channel);
        for (const auto& wire_id : wire_ids) {
            bad_wires_by_plane[PlaneKey {static_cast<int>(wire_id.Cryostat),
                                         static_cast<int>(wire_id.TPC),
                                         static_cast<int>(wire_id.Plane)}]
                .push_back(static_cast<int>(wire_id.Wire));
        }
    }

    for (auto& item : bad_wires_by_plane) {
        auto& wires = item.second;
        if (wires.empty()) continue;

        std::sort(wires.begin(), wires.end());
        wires.erase(std::unique(wires.begin(), wires.end()), wires.end());

        int first = wires.front();
        int last = wires.front();
        auto& ranges = fBadWireRanges[item.first];
        for (std::size_t i = 1u; i < wires.size(); ++i) {
            if (wires[i] <= last + 1) {
                last = wires[i];
                continue;
            }
            ranges.emplace_back(first, last);
            first = last = wires[i];
        }
        ranges.emplace_back(first, last);
    }

    fBadWireRangesBuilt = true;
}

void GraphAnalysis::read_config(const fhicl::ParameterSet& p) {
    fHITproducer = p.get<art::InputTag>("HITproducer", "gaushit");
    fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
    fOnlySelectedSlices = p.get<bool>("OnlySelectedSlices", true);
    fBadChannelFile = p.get<std::string>("BadChannelFile", "");
    fBadChannels.clear();
    fBadWireRanges.clear();
    fBadWireRangesBuilt = false;
    if (!fBadChannelFile.empty()) this->load_bad_channels(fBadChannelFile);

    fNominalParams.save_all_component_graph_edges =
        p.get<bool>("SaveAllComponentGraphEdges", false);
    fNominalParams.merge_detached_activities =
        p.get<bool>("MergeDetachedActivities", true);

    fNominalParams.min_component_hits =
        p.get<unsigned>("MinComponentHits", 6u);
    fNominalParams.boundary_candidates_per_side =
        p.get<unsigned>("BoundaryCandidatesPerSide", 3u);
    fNominalParams.max_geom_neighbors =
        p.get<unsigned>("MaxGeomNeighbors", 16u);
    fNominalParams.max_candidate_neighbors =
        p.get<unsigned>("MaxCandidateNeighbors", 12u);

    fNominalParams.cell_size_scale =
        p.get<float>("CellSizeScale", 1.0f);
    fNominalParams.candidate_radius_scale =
        p.get<float>("CandidateRadiusScale", 2.5f);
    fNominalParams.geom_radius_scale =
        p.get<float>("GeomRadiusScale", 3.5f);
    fNominalParams.capsule_radius_scale =
        p.get<float>("CapsuleRadiusScale", 1.0f);
    fNominalParams.tube_radius_scale =
        p.get<float>("TubeRadiusScale", 1.0f);
    fNominalParams.start_radius_scale =
        p.get<float>("StartRadiusScale", 4.0f);
    fNominalParams.min_gap_scale =
        p.get<float>("MinGapScale", 2.0f);
    fNominalParams.max_gap_scale =
        p.get<float>("MaxGapScale", 30.0f);
    fNominalParams.parent_search_radius_scale =
        p.get<float>("ParentSearchRadiusScale", 30.0f);

    fNominalParams.alpha_dist = p.get<float>("AlphaDist", 1.0f);
    fNominalParams.beta_angle = p.get<float>("BetaAngle", 1.0f);
    fNominalParams.gamma_support = p.get<float>("GammaSupport", 1.0f);
    fNominalParams.delta_charge = p.get<float>("DeltaCharge", 0.25f);
    fNominalParams.tau_ch = p.get<float>("TauCh", 2.0f);

    fNominalParams.cone_angle_deg = p.get<float>("ConeAngleDeg", 35.0f);
    fNominalParams.start_kappa = p.get<float>("StartKappa", 1.5f);
    fNominalParams.min_start_asym = p.get<float>("MinStartAsym", 0.70f);

    fNominalParams.min_pointing = p.get<float>("MinPointing", 0.70f);
    fNominalParams.min_emptiness = p.get<float>("MinEmptiness", 0.85f);
    fNominalParams.min_visibility = p.get<float>("MinVisibility", 0.80f);
    fNominalParams.min_path_break = p.get<float>("MinPathBreak", 0.70f);

    fNominalParams.capsule_interior_frac =
        cg::clamp_value(p.get<float>("CapsuleInteriorFrac", 0.10f), 0.0f, 0.49f);
    fNominalParams.corridor_end_exclude_frac =
        cg::clamp_value(p.get<float>("CorridorEndExcludeFrac", 0.15f), 0.0f, 0.49f);

    fNominalParams.max_secondary_vertex_residual =
        p.get<float>("MaxSecondaryVertexResidual", 2.0f);
    fNominalParams.max_secondary_vertex_back_projection =
        p.get<float>("MaxSecondaryVertexBackProjection", 2.0f);

    fEngine.set_params(fNominalParams);
}

void GraphAnalysis::configure(const fhicl::ParameterSet& parameter_set) {
    this->read_config(parameter_set);
}

void GraphAnalysis::analyseEvent(const art::Event& event, bool /*is_data*/) {
    _evt = event.event();
    _sbr = event.subRun();
    _run = event.run();
    _da_slice_counter = 0;
    _da_next_plane_uid = 0;
    _da_next_component_id = 0;
}

cg::PlaneGeometry GraphAnalysis::make_plane_geometry(
    const geo::PlaneID& plane_id,
    const detinfo::DetectorProperties* det_prop,
    const geo::GeometryCore* geom) const {

    cg::PlaneGeometry out;
    out.wire_pitch = static_cast<float>(geom->WirePitch(plane_id));

    const double x0 = det_prop->ConvertTicksToX(0.0, plane_id.Plane, plane_id.TPC, plane_id.Cryostat);
    const double x1 = det_prop->ConvertTicksToX(1.0, plane_id.Plane, plane_id.TPC, plane_id.Cryostat);
    out.drift_cm_per_tick = static_cast<float>(std::fabs(x1 - x0));

    const auto& tpc_geo = geom->TPC(plane_id.asTPCID());
    const auto& plane_geo = tpc_geo.Plane(plane_id);

    out.wire_min = -0.5f * out.wire_pitch;
    out.wire_max = static_cast<float>((plane_geo.Nwires() - 0.5) * out.wire_pitch);

    const unsigned nticks = det_prop->NumberTimeSamples();
    const double xd0 = det_prop->ConvertTicksToX(0.0, plane_id.Plane, plane_id.TPC, plane_id.Cryostat);
    const double xd1 = det_prop->ConvertTicksToX(static_cast<double>(nticks - 1u),
                                                plane_id.Plane,
                                                plane_id.TPC,
                                                plane_id.Cryostat);

    out.drift_min = static_cast<float>(std::min(xd0, xd1));
    out.drift_max = static_cast<float>(std::max(xd0, xd1));

    const PlaneKey key {static_cast<int>(plane_id.Cryostat),
                        static_cast<int>(plane_id.TPC),
                        static_cast<int>(plane_id.Plane)};
    const auto range_iter = fBadWireRanges.find(key);
    if (range_iter != fBadWireRanges.end()) {
        out.dead_wire_intervals.reserve(range_iter->second.size());
        for (const auto& range : range_iter->second) {
            const float dead_min = std::max((static_cast<float>(range.first) - 0.5f) * out.wire_pitch,
                                            out.wire_min);
            const float dead_max = std::min((static_cast<float>(range.second) + 0.5f) * out.wire_pitch,
                                            out.wire_max);
            if (dead_max > dead_min) {
                out.dead_wire_intervals.emplace_back(dead_min, dead_max);
            }
        }
    }

    return out;
}

void GraphAnalysis::collect_slice_hits(
    const art::Event& event,
    const int slice_index,
    std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
    const detinfo::DetectorProperties* det_prop,
    const geo::GeometryCore* geom,
    std::map<PlaneKey, std::vector<cg::HitInput>>& hits_by_plane) const {

    std::unordered_set<std::size_t> used_hit_keys;

    auto const& cluster_handle = event.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
    art::FindManyP<recob::Hit> cluster_hits(cluster_handle, event, fCLSproducer);

    for (const auto& pfp : slice_pfp_vec) {
        const auto clusters = pfp.get<recob::Cluster>();
        for (const auto& cluster : clusters) {
            const auto& hits = cluster_hits.at(cluster.key());
            for (const auto& hit : hits) {
                if (!hit) continue;
                if (!used_hit_keys.insert(hit.key()).second) continue;

                const auto& wire_id = hit->WireID();
                const auto plane_id = wire_id.planeID();

                cg::HitInput node;
                node.hit_key = hit.key();
                node.slice_index = slice_index;
                node.cryo = static_cast<int>(wire_id.Cryostat);
                node.tpc = static_cast<int>(wire_id.TPC);
                node.plane = static_cast<int>(wire_id.Plane);
                node.wire = static_cast<int>(wire_id.Wire);
                node.tick = hit->PeakTime();
                node.charge = hit->Integral();
                node.peakamp = hit->PeakAmplitude();
                node.rms = hit->RMS();

                const float wire_pitch = static_cast<float>(geom->WirePitch(plane_id));
                node.wire_cm = static_cast<float>(wire_id.Wire) * wire_pitch;
                node.drift_cm = static_cast<float>(
                    det_prop->ConvertTicksToX(hit->PeakTime(),
                                              plane_id.Plane,
                                              plane_id.TPC,
                                              plane_id.Cryostat));

                hits_by_plane[PlaneKey {node.cryo, node.tpc, node.plane}].push_back(node);
            }
        }
    }
}

void GraphAnalysis::write_plane_replay(const int plane_uid,
                                                  const cg::PlaneInput& input) {
    _da_plane_uid.push_back(plane_uid);
    _da_plane_slice_index.push_back(input.slice_index);
    _da_plane_cryo.push_back(input.cryo);
    _da_plane_tpc.push_back(input.tpc);
    _da_plane_plane.push_back(input.plane);
    _da_plane_wire_pitch_cm.push_back(input.geometry.wire_pitch);
    _da_plane_drift_cm_per_tick.push_back(input.geometry.drift_cm_per_tick);
    _da_plane_wire_min_cm.push_back(input.geometry.wire_min);
    _da_plane_wire_max_cm.push_back(input.geometry.wire_max);
    _da_plane_drift_min_cm.push_back(input.geometry.drift_min);
    _da_plane_drift_max_cm.push_back(input.geometry.drift_max);

    for (const auto& interval : input.geometry.dead_wire_intervals) {
        _da_dead_plane_uid.push_back(plane_uid);
        _da_dead_slice_index.push_back(input.slice_index);
        _da_dead_cryo.push_back(input.cryo);
        _da_dead_tpc.push_back(input.tpc);
        _da_dead_plane.push_back(input.plane);
        _da_dead_wire_min_cm.push_back(interval.first);
        _da_dead_wire_max_cm.push_back(interval.second);
        _da_dead_drift_min_cm.push_back(input.geometry.drift_min);
        _da_dead_drift_max_cm.push_back(input.geometry.drift_max);
    }
}

bool GraphAnalysis::find_slice_primary_vertex(
    std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
    SlicePV& out) const {

    for (const auto& pfp_proxy : slice_pfp_vec) {
        const auto pfps = pfp_proxy.get<recob::PFParticle>();
        if (pfps.size() == 0u || !pfps.front()) continue;

        const auto& pfp = pfps.front();
        const int apdg = std::abs(pfp->PdgCode());

        if (!pfp->IsPrimary()) continue;
        if (apdg != 12 && apdg != 14 && apdg != 16) continue;

        const auto vertices = pfp_proxy.get<recob::Vertex>();
        if (vertices.size() == 0u || !vertices.front() || !vertices.front()->isValid()) continue;

        out.pos = vertices.front()->position();
        out.valid = true;
        return true;
    }

    return false;
}

bool GraphAnalysis::project_pv_to_plane(recob::Vertex::Point_t const& pos,
                                        geo::PlaneID const& plane_id,
                                        geo::GeometryCore const* geom,
                                        float& wire_cm,
                                        float& drift_cm) const {
    wire_cm = std::numeric_limits<float>::quiet_NaN();
    drift_cm = std::numeric_limits<float>::quiet_NaN();

    try {
        const auto tpc_id = geom->FindTPCAtPosition(pos);
        if (tpc_id.Cryostat != plane_id.Cryostat || tpc_id.TPC != plane_id.TPC) return false;
    }
    catch (...) {
        return false;
    }

    const double w = geom->WireCoordinate(pos, plane_id);
    wire_cm = static_cast<float>(w * geom->WirePitch(plane_id));
    drift_cm = static_cast<float>(pos.X());
    return true;
}

void GraphAnalysis::write_nominal_output(const int plane_uid,
                                                    const cg::Result& result,
                                                    const bool have_plane_pv,
                                                    const float pv_wire_cm,
                                                    const float pv_drift_cm) {
    std::vector<int> global_comp_id(result.components.size(), -1);
    for (const auto& comp : result.components) {
        if (comp.id < 0 || comp.id >= static_cast<int>(global_comp_id.size())) continue;
        global_comp_id[comp.id] = _da_next_component_id++;
    }

    std::vector<int> global_activity_id(global_comp_id.size(), -1);
    for (const auto& comp : result.components) {
        if (comp.activity_id < 0 || comp.activity_id >= static_cast<int>(global_comp_id.size())) continue;
        const int rep_global = global_comp_id[comp.activity_id];
        if (rep_global < 0) continue;
        global_activity_id[comp.activity_id] = rep_global;
    }

    auto resolve_global_comp = [&](const int local_id) -> int {
        if (local_id < 0 || local_id >= static_cast<int>(global_comp_id.size())) return -1;
        return global_comp_id[local_id];
    };

    auto resolve_global_activity = [&](const int local_activity_id) -> int {
        if (local_activity_id < 0 || local_activity_id >= static_cast<int>(global_activity_id.size())) return -1;
        return global_activity_id[local_activity_id];
    };

    for (const auto& hit : result.hits) {
        _da_hit_plane_uid.push_back(plane_uid);
        _da_hit_slice_index.push_back(hit.slice_index);
        _da_hit_key.push_back(hit.hit_key);
        _da_hit_cryo.push_back(hit.cryo);
        _da_hit_tpc.push_back(hit.tpc);
        _da_hit_plane.push_back(hit.plane);
        _da_hit_wire.push_back(hit.wire);
        _da_hit_tick.push_back(hit.tick);
        _da_hit_wire_cm.push_back(hit.wire_cm);
        _da_hit_drift_cm.push_back(hit.drift_cm);
        _da_hit_charge.push_back(hit.charge);
        _da_hit_peakamp.push_back(hit.peakamp);
        _da_hit_rms.push_back(hit.rms);
        _da_hit_dir_wire.push_back(hit.dir.x);
        _da_hit_dir_drift.push_back(hit.dir.y);
        _da_hit_anis.push_back(hit.anis);
        _da_hit_scale.push_back(hit.local_scale);
        _da_hit_component_id.push_back(resolve_global_comp(hit.component_id));
        _da_hit_activity_id.push_back(resolve_global_activity(hit.activity_id));
        _da_hit_is_boundary.push_back(hit.is_boundary ? 1 : 0);
        _da_hit_is_start.push_back(hit.is_start ? 1 : 0);
        _da_hit_label.push_back(hit.label);
    }

    for (const auto& comp : result.components) {
        float pv_to_activity = std::numeric_limits<float>::quiet_NaN();
        const bool is_activity_representative = (comp.id == comp.activity_id);
        if (have_plane_pv && comp.label == cg::kDetached && is_activity_representative) {
            const float dw = comp.activity_vertex.x - pv_wire_cm;
            const float dx = comp.activity_vertex.y - pv_drift_cm;
            pv_to_activity = std::hypot(dw, dx);
        }

        float start_wire = -1.0f;
        float start_tick = -1.0f;
        if (comp.start_hit >= 0 && comp.start_hit < static_cast<int>(result.hits.size())) {
            start_wire = static_cast<float>(result.hits[comp.start_hit].wire);
            start_tick = result.hits[comp.start_hit].tick;
        }

        _da_comp_plane_uid.push_back(plane_uid);
        _da_comp_slice_index.push_back(result.input.slice_index);
        _da_comp_id.push_back(resolve_global_comp(comp.id));
        _da_comp_activity_id.push_back(resolve_global_activity(comp.activity_id));
        _da_comp_cryo.push_back(comp.cryo);
        _da_comp_tpc.push_back(comp.tpc);
        _da_comp_plane.push_back(comp.plane);
        _da_comp_n_hits.push_back(static_cast<int>(comp.hit_indices.size()));
        _da_comp_label.push_back(comp.label);
        _da_comp_parent_id.push_back(resolve_global_comp(comp.parent_id));
        _da_comp_qsum.push_back(comp.qsum);
        _da_comp_start_wire.push_back(start_wire);
        _da_comp_start_tick.push_back(start_tick);
        _da_comp_start_wire_cm.push_back(comp.start_pos.x);
        _da_comp_start_drift_cm.push_back(comp.start_pos.y);
        _da_comp_dir_wire.push_back(comp.start_dir.x);
        _da_comp_dir_drift.push_back(comp.start_dir.y);
        _da_comp_start_asym.push_back(comp.start_asym);
        _da_comp_bbox_min_wire_cm.push_back(comp.bbox_min_wire);
        _da_comp_bbox_max_wire_cm.push_back(comp.bbox_max_wire);
        _da_comp_bbox_min_drift_cm.push_back(comp.bbox_min_drift);
        _da_comp_bbox_max_drift_cm.push_back(comp.bbox_max_drift);
        _da_comp_gap_cm.push_back(comp.gap);
        _da_comp_emptiness.push_back(comp.emptiness);
        _da_comp_visibility.push_back(comp.visibility);
        _da_comp_pointing.push_back(comp.pointing);
        _da_comp_path_break.push_back(comp.path_break);
        _da_comp_score.push_back(comp.score);
        _da_comp_activity_vtx_wire_cm.push_back(comp.activity_vertex.x);
        _da_comp_activity_vtx_drift_cm.push_back(comp.activity_vertex.y);
        _da_comp_activity_vtx_residual.push_back(comp.activity_residual);
        _da_comp_pv_to_activity_dist2d_cm.push_back(pv_to_activity);
    }

    for (const auto& edge : result.edges) {
        _da_edge_plane_uid.push_back(plane_uid);
        _da_edge_slice_index.push_back(result.input.slice_index);
        _da_edge_plane.push_back(edge.plane);
        _da_edge_src_component.push_back(resolve_global_comp(edge.src_component));
        _da_edge_dst_component.push_back(resolve_global_comp(edge.dst_component));
        _da_edge_kind.push_back(edge.kind);
        _da_edge_src_wire_cm.push_back(edge.src_pos.x);
        _da_edge_src_drift_cm.push_back(edge.src_pos.y);
        _da_edge_dst_wire_cm.push_back(edge.dst_pos.x);
        _da_edge_dst_drift_cm.push_back(edge.dst_pos.y);
        _da_edge_gap_cm.push_back(edge.gap);
        _da_edge_emptiness.push_back(edge.emptiness);
        _da_edge_visibility.push_back(edge.visibility);
        _da_edge_pointing.push_back(edge.pointing);
        _da_edge_path_break.push_back(edge.path_break);
        _da_edge_score.push_back(edge.score);
    }
}

void GraphAnalysis::analyseSlice(const art::Event& event,
                                            std::vector<common::ProxyPfpElem_t>& slice_pfp_vec,
                                            bool /*is_data*/,
                                            bool is_selected) {
    if (fOnlySelectedSlices && !is_selected) {
        ++_da_slice_counter;
        return;
    }

    auto const* det_prop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const* geom = lar::providerFrom<geo::Geometry>();
    SlicePV slice_pv;
    this->find_slice_primary_vertex(slice_pfp_vec, slice_pv);

    if (!fBadChannelFile.empty() && !fBadWireRangesBuilt) {
        this->build_bad_wire_ranges(geom);
    }

    std::map<PlaneKey, std::vector<cg::HitInput>> hits_by_plane;
    this->collect_slice_hits(event, _da_slice_counter, slice_pfp_vec, det_prop, geom, hits_by_plane);

    for (auto& plane_hits_pair : hits_by_plane) {
        auto& plane_hits = plane_hits_pair.second;
        if (plane_hits.empty()) continue;

        const auto& first_hit = plane_hits.front();
        const geo::PlaneID plane_id(static_cast<unsigned int>(first_hit.cryo),
                                    static_cast<unsigned int>(first_hit.tpc),
                                    static_cast<unsigned int>(first_hit.plane));

        cg::PlaneInput input;
        input.slice_index = _da_slice_counter;
        input.cryo = first_hit.cryo;
        input.tpc = first_hit.tpc;
        input.plane = first_hit.plane;
        input.geometry = this->make_plane_geometry(plane_id, det_prop, geom);
        input.hits = plane_hits;

        const int plane_uid = _da_next_plane_uid++;
        float pv_wire_cm = std::numeric_limits<float>::quiet_NaN();
        float pv_drift_cm = std::numeric_limits<float>::quiet_NaN();
        const bool have_plane_pv =
            slice_pv.valid && this->project_pv_to_plane(slice_pv.pos, plane_id, geom, pv_wire_cm, pv_drift_cm);

        _da_plane_has_primary_vtx.push_back(have_plane_pv ? 1 : 0);
        _da_plane_primary_vtx_wire_cm.push_back(pv_wire_cm);
        _da_plane_primary_vtx_drift_cm.push_back(pv_drift_cm);

        this->write_plane_replay(plane_uid, input);

        const cg::Result nominal = fEngine.run(input);
        this->write_nominal_output(plane_uid, nominal, have_plane_pv, pv_wire_cm, pv_drift_cm);
    }

    ++_da_slice_counter;
}

void GraphAnalysis::setBranches(TTree* tree) {
    tree->Branch("da_plane_uid", "std::vector<int>", &_da_plane_uid);
    tree->Branch("da_plane_slice_index", "std::vector<int>", &_da_plane_slice_index);
    tree->Branch("da_plane_cryo", "std::vector<int>", &_da_plane_cryo);
    tree->Branch("da_plane_tpc", "std::vector<int>", &_da_plane_tpc);
    tree->Branch("da_plane_plane", "std::vector<int>", &_da_plane_plane);
    tree->Branch("da_plane_wire_pitch_cm", "std::vector<float>", &_da_plane_wire_pitch_cm);
    tree->Branch("da_plane_drift_cm_per_tick", "std::vector<float>", &_da_plane_drift_cm_per_tick);
    tree->Branch("da_plane_wire_min_cm", "std::vector<float>", &_da_plane_wire_min_cm);
    tree->Branch("da_plane_wire_max_cm", "std::vector<float>", &_da_plane_wire_max_cm);
    tree->Branch("da_plane_drift_min_cm", "std::vector<float>", &_da_plane_drift_min_cm);
    tree->Branch("da_plane_drift_max_cm", "std::vector<float>", &_da_plane_drift_max_cm);
    tree->Branch("da_plane_has_primary_vtx", "std::vector<int>", &_da_plane_has_primary_vtx);
    tree->Branch("da_plane_primary_vtx_wire_cm", "std::vector<float>", &_da_plane_primary_vtx_wire_cm);
    tree->Branch("da_plane_primary_vtx_drift_cm", "std::vector<float>", &_da_plane_primary_vtx_drift_cm);

    tree->Branch("da_dead_plane_uid", "std::vector<int>", &_da_dead_plane_uid);
    tree->Branch("da_dead_slice_index", "std::vector<int>", &_da_dead_slice_index);
    tree->Branch("da_dead_cryo", "std::vector<int>", &_da_dead_cryo);
    tree->Branch("da_dead_tpc", "std::vector<int>", &_da_dead_tpc);
    tree->Branch("da_dead_plane", "std::vector<int>", &_da_dead_plane);
    tree->Branch("da_dead_wire_min_cm", "std::vector<float>", &_da_dead_wire_min_cm);
    tree->Branch("da_dead_wire_max_cm", "std::vector<float>", &_da_dead_wire_max_cm);
    tree->Branch("da_dead_drift_min_cm", "std::vector<float>", &_da_dead_drift_min_cm);
    tree->Branch("da_dead_drift_max_cm", "std::vector<float>", &_da_dead_drift_max_cm);

    tree->Branch("da_hit_plane_uid", "std::vector<int>", &_da_hit_plane_uid);
    tree->Branch("da_hit_slice_index", "std::vector<int>", &_da_hit_slice_index);
    tree->Branch("da_hit_key", "std::vector<size_t>", &_da_hit_key);
    tree->Branch("da_hit_cryo", "std::vector<int>", &_da_hit_cryo);
    tree->Branch("da_hit_tpc", "std::vector<int>", &_da_hit_tpc);
    tree->Branch("da_hit_plane", "std::vector<int>", &_da_hit_plane);
    tree->Branch("da_hit_wire", "std::vector<int>", &_da_hit_wire);
    tree->Branch("da_hit_tick", "std::vector<float>", &_da_hit_tick);
    tree->Branch("da_hit_wire_cm", "std::vector<float>", &_da_hit_wire_cm);
    tree->Branch("da_hit_drift_cm", "std::vector<float>", &_da_hit_drift_cm);
    tree->Branch("da_hit_charge", "std::vector<float>", &_da_hit_charge);
    tree->Branch("da_hit_peakamp", "std::vector<float>", &_da_hit_peakamp);
    tree->Branch("da_hit_rms", "std::vector<float>", &_da_hit_rms);
    tree->Branch("da_hit_dir_wire", "std::vector<float>", &_da_hit_dir_wire);
    tree->Branch("da_hit_dir_drift", "std::vector<float>", &_da_hit_dir_drift);
    tree->Branch("da_hit_anis", "std::vector<float>", &_da_hit_anis);
    tree->Branch("da_hit_scale", "std::vector<float>", &_da_hit_scale);
    tree->Branch("da_hit_component_id", "std::vector<int>", &_da_hit_component_id);
    tree->Branch("da_hit_activity_id", "std::vector<int>", &_da_hit_activity_id);
    tree->Branch("da_hit_is_boundary", "std::vector<int>", &_da_hit_is_boundary);
    tree->Branch("da_hit_is_start", "std::vector<int>", &_da_hit_is_start);
    tree->Branch("da_hit_label", "std::vector<int>", &_da_hit_label);

    tree->Branch("da_comp_plane_uid", "std::vector<int>", &_da_comp_plane_uid);
    tree->Branch("da_comp_slice_index", "std::vector<int>", &_da_comp_slice_index);
    tree->Branch("da_comp_id", "std::vector<int>", &_da_comp_id);
    tree->Branch("da_comp_activity_id", "std::vector<int>", &_da_comp_activity_id);
    tree->Branch("da_comp_cryo", "std::vector<int>", &_da_comp_cryo);
    tree->Branch("da_comp_tpc", "std::vector<int>", &_da_comp_tpc);
    tree->Branch("da_comp_plane", "std::vector<int>", &_da_comp_plane);
    tree->Branch("da_comp_n_hits", "std::vector<int>", &_da_comp_n_hits);
    tree->Branch("da_comp_label", "std::vector<int>", &_da_comp_label);
    tree->Branch("da_comp_parent_id", "std::vector<int>", &_da_comp_parent_id);
    tree->Branch("da_comp_qsum", "std::vector<float>", &_da_comp_qsum);
    tree->Branch("da_comp_start_wire", "std::vector<float>", &_da_comp_start_wire);
    tree->Branch("da_comp_start_tick", "std::vector<float>", &_da_comp_start_tick);
    tree->Branch("da_comp_start_wire_cm", "std::vector<float>", &_da_comp_start_wire_cm);
    tree->Branch("da_comp_start_drift_cm", "std::vector<float>", &_da_comp_start_drift_cm);
    tree->Branch("da_comp_dir_wire", "std::vector<float>", &_da_comp_dir_wire);
    tree->Branch("da_comp_dir_drift", "std::vector<float>", &_da_comp_dir_drift);
    tree->Branch("da_comp_start_asym", "std::vector<float>", &_da_comp_start_asym);
    tree->Branch("da_comp_bbox_min_wire_cm", "std::vector<float>", &_da_comp_bbox_min_wire_cm);
    tree->Branch("da_comp_bbox_max_wire_cm", "std::vector<float>", &_da_comp_bbox_max_wire_cm);
    tree->Branch("da_comp_bbox_min_drift_cm", "std::vector<float>", &_da_comp_bbox_min_drift_cm);
    tree->Branch("da_comp_bbox_max_drift_cm", "std::vector<float>", &_da_comp_bbox_max_drift_cm);
    tree->Branch("da_comp_gap_cm", "std::vector<float>", &_da_comp_gap_cm);
    tree->Branch("da_comp_emptiness", "std::vector<float>", &_da_comp_emptiness);
    tree->Branch("da_comp_visibility", "std::vector<float>", &_da_comp_visibility);
    tree->Branch("da_comp_pointing", "std::vector<float>", &_da_comp_pointing);
    tree->Branch("da_comp_path_break", "std::vector<float>", &_da_comp_path_break);
    tree->Branch("da_comp_score", "std::vector<float>", &_da_comp_score);
    tree->Branch("da_comp_activity_vtx_wire_cm", "std::vector<float>", &_da_comp_activity_vtx_wire_cm);
    tree->Branch("da_comp_activity_vtx_drift_cm", "std::vector<float>", &_da_comp_activity_vtx_drift_cm);
    tree->Branch("da_comp_activity_vtx_residual", "std::vector<float>", &_da_comp_activity_vtx_residual);
    tree->Branch("da_comp_pv_to_activity_dist2d_cm", "std::vector<float>", &_da_comp_pv_to_activity_dist2d_cm);

    tree->Branch("da_edge_plane_uid", "std::vector<int>", &_da_edge_plane_uid);
    tree->Branch("da_edge_slice_index", "std::vector<int>", &_da_edge_slice_index);
    tree->Branch("da_edge_plane", "std::vector<int>", &_da_edge_plane);
    tree->Branch("da_edge_src_component", "std::vector<int>", &_da_edge_src_component);
    tree->Branch("da_edge_dst_component", "std::vector<int>", &_da_edge_dst_component);
    tree->Branch("da_edge_kind", "std::vector<int>", &_da_edge_kind);
    tree->Branch("da_edge_src_wire_cm", "std::vector<float>", &_da_edge_src_wire_cm);
    tree->Branch("da_edge_src_drift_cm", "std::vector<float>", &_da_edge_src_drift_cm);
    tree->Branch("da_edge_dst_wire_cm", "std::vector<float>", &_da_edge_dst_wire_cm);
    tree->Branch("da_edge_dst_drift_cm", "std::vector<float>", &_da_edge_dst_drift_cm);
    tree->Branch("da_edge_gap_cm", "std::vector<float>", &_da_edge_gap_cm);
    tree->Branch("da_edge_emptiness", "std::vector<float>", &_da_edge_emptiness);
    tree->Branch("da_edge_visibility", "std::vector<float>", &_da_edge_visibility);
    tree->Branch("da_edge_pointing", "std::vector<float>", &_da_edge_pointing);
    tree->Branch("da_edge_path_break", "std::vector<float>", &_da_edge_path_break);
    tree->Branch("da_edge_score", "std::vector<float>", &_da_edge_score);
}

void GraphAnalysis::resetTTree(TTree* /*tree*/) {
    _da_plane_uid.clear();
    _da_plane_slice_index.clear();
    _da_plane_cryo.clear();
    _da_plane_tpc.clear();
    _da_plane_plane.clear();
    _da_plane_wire_pitch_cm.clear();
    _da_plane_drift_cm_per_tick.clear();
    _da_plane_wire_min_cm.clear();
    _da_plane_wire_max_cm.clear();
    _da_plane_drift_min_cm.clear();
    _da_plane_drift_max_cm.clear();
    _da_plane_has_primary_vtx.clear();
    _da_plane_primary_vtx_wire_cm.clear();
    _da_plane_primary_vtx_drift_cm.clear();

    _da_dead_plane_uid.clear();
    _da_dead_slice_index.clear();
    _da_dead_cryo.clear();
    _da_dead_tpc.clear();
    _da_dead_plane.clear();
    _da_dead_wire_min_cm.clear();
    _da_dead_wire_max_cm.clear();
    _da_dead_drift_min_cm.clear();
    _da_dead_drift_max_cm.clear();

    _da_hit_plane_uid.clear();
    _da_hit_slice_index.clear();
    _da_hit_key.clear();
    _da_hit_cryo.clear();
    _da_hit_tpc.clear();
    _da_hit_plane.clear();
    _da_hit_wire.clear();
    _da_hit_tick.clear();
    _da_hit_wire_cm.clear();
    _da_hit_drift_cm.clear();
    _da_hit_charge.clear();
    _da_hit_peakamp.clear();
    _da_hit_rms.clear();
    _da_hit_dir_wire.clear();
    _da_hit_dir_drift.clear();
    _da_hit_anis.clear();
    _da_hit_scale.clear();
    _da_hit_component_id.clear();
    _da_hit_activity_id.clear();
    _da_hit_is_boundary.clear();
    _da_hit_is_start.clear();
    _da_hit_label.clear();

    _da_comp_plane_uid.clear();
    _da_comp_slice_index.clear();
    _da_comp_id.clear();
    _da_comp_activity_id.clear();
    _da_comp_cryo.clear();
    _da_comp_tpc.clear();
    _da_comp_plane.clear();
    _da_comp_n_hits.clear();
    _da_comp_label.clear();
    _da_comp_parent_id.clear();
    _da_comp_qsum.clear();
    _da_comp_start_wire.clear();
    _da_comp_start_tick.clear();
    _da_comp_start_wire_cm.clear();
    _da_comp_start_drift_cm.clear();
    _da_comp_dir_wire.clear();
    _da_comp_dir_drift.clear();
    _da_comp_start_asym.clear();
    _da_comp_bbox_min_wire_cm.clear();
    _da_comp_bbox_max_wire_cm.clear();
    _da_comp_bbox_min_drift_cm.clear();
    _da_comp_bbox_max_drift_cm.clear();
    _da_comp_gap_cm.clear();
    _da_comp_emptiness.clear();
    _da_comp_visibility.clear();
    _da_comp_pointing.clear();
    _da_comp_path_break.clear();
    _da_comp_score.clear();
    _da_comp_activity_vtx_wire_cm.clear();
    _da_comp_activity_vtx_drift_cm.clear();
    _da_comp_activity_vtx_residual.clear();
    _da_comp_pv_to_activity_dist2d_cm.clear();

    _da_edge_plane_uid.clear();
    _da_edge_slice_index.clear();
    _da_edge_plane.clear();
    _da_edge_src_component.clear();
    _da_edge_dst_component.clear();
    _da_edge_kind.clear();
    _da_edge_src_wire_cm.clear();
    _da_edge_src_drift_cm.clear();
    _da_edge_dst_wire_cm.clear();
    _da_edge_dst_drift_cm.clear();
    _da_edge_gap_cm.clear();
    _da_edge_emptiness.clear();
    _da_edge_visibility.clear();
    _da_edge_pointing.clear();
    _da_edge_path_break.clear();
    _da_edge_score.clear();

    _da_slice_counter = 0;
    _da_next_plane_uid = 0;
    _da_next_component_id = 0;
}

DEFINE_ART_CLASS_TOOL(GraphAnalysis)

} // namespace analysis

#endif
