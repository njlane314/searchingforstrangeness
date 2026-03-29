#ifndef COMMON_GRAPHENGINE_H
#define COMMON_GRAPHENGINE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace analysis {
namespace cg {

constexpr float kPi = 3.14159265358979323846f;
constexpr float kEps = 1.0e-5f;

inline float qnanf() {
    return std::numeric_limits<float>::quiet_NaN();
}

template <typename T>
inline T clamp_value(const T& x, const T& lo, const T& hi) {
    return std::max(lo, std::min(x, hi));
}

enum Label : int {
    kPrompt = 0,
    kDetached = 1,
    kUnresolved = 2
};

struct Vec2 {
    float x = 0.0f;
    float y = 0.0f;
};

struct Params {
    bool merge_detached_activities = true;
    bool save_all_component_graph_edges = false;

    unsigned min_component_hits = 6u;
    unsigned boundary_candidates_per_side = 3u;
    unsigned max_geom_neighbors = 16u;
    unsigned max_candidate_neighbors = 12u;

    float cell_size_scale = 1.0f;
    float candidate_radius_scale = 2.5f;
    float geom_radius_scale = 3.5f;
    float capsule_radius_scale = 1.0f;
    float tube_radius_scale = 1.0f;
    float start_radius_scale = 4.0f;
    float min_gap_scale = 2.0f;
    float max_gap_scale = 30.0f;
    float parent_search_radius_scale = 30.0f;

    float alpha_dist = 1.0f;
    float beta_angle = 1.0f;
    float gamma_support = 1.0f;
    float delta_charge = 0.25f;
    float tau_ch = 2.0f;

    float cone_angle_deg = 35.0f;
    float start_kappa = 1.5f;
    float min_start_asym = 0.70f;

    float min_pointing = 0.70f;
    float min_emptiness = 0.85f;
    float min_visibility = 0.80f;
    float min_path_break = 0.70f;

    float capsule_interior_frac = 0.10f;
    float corridor_end_exclude_frac = 0.15f;

    float max_secondary_vertex_residual = 2.0f;
    float max_secondary_vertex_back_projection = 2.0f;
};

struct PlaneGeometry {
    float wire_pitch = 1.0f;
    float drift_cm_per_tick = 1.0f;
    float wire_min = 0.0f;
    float wire_max = 0.0f;
    float drift_min = 0.0f;
    float drift_max = 0.0f;
    std::vector<std::pair<float, float>> dead_wire_intervals;
};

struct DerivedScales {
    float s = 1.0f;
    float cell_size = 1.0f;
    float candidate_radius = 2.0f;
    float geom_radius = 3.0f;
    float capsule_radius = 1.0f;
    float tube_radius = 1.0f;
    float start_radius = 4.0f;
    float min_gap = 2.0f;
    float max_gap = 30.0f;
    float parent_search_radius = 30.0f;
    float wire_min = 0.0f;
    float wire_max = 0.0f;
    float drift_min = 0.0f;
    float drift_max = 0.0f;
    std::vector<std::pair<float, float>> dead_wire_intervals;
};

struct HitInput {
    std::size_t hit_key = 0u;
    int slice_index = -1;
    int cryo = 0;
    int tpc = 0;
    int plane = 0;
    int wire = 0;
    float tick = 0.0f;
    float wire_cm = 0.0f;
    float drift_cm = 0.0f;
    float charge = 0.0f;
    float peakamp = 0.0f;
    float rms = 0.0f;
};

struct PlaneInput {
    int slice_index = -1;
    int cryo = 0;
    int tpc = 0;
    int plane = 0;
    PlaneGeometry geometry;
    std::vector<HitInput> hits;
};

struct HitResult : public HitInput {
    Vec2 dir {1.0f, 0.0f};
    float anis = 0.0f;
    float local_scale = 1.0f;
    int component_id = -1;
    int activity_id = -1;
    bool is_boundary = false;
    bool is_start = false;
    int label = kUnresolved;
};

struct ComponentResult {
    int id = -1;
    int activity_id = -1;
    int label = kPrompt;
    int cryo = 0;
    int tpc = 0;
    int plane = 0;
    std::vector<int> hit_indices;
    int start_hit = -1;
    Vec2 start_pos;
    Vec2 start_dir {1.0f, 0.0f};
    float start_asym = 0.0f;
    float qsum = 0.0f;
    float bbox_min_wire = 0.0f;
    float bbox_max_wire = 0.0f;
    float bbox_min_drift = 0.0f;
    float bbox_max_drift = 0.0f;
    int parent_id = -1;
    Vec2 parent_anchor_pos;
    float gap = qnanf();
    float emptiness = qnanf();
    float visibility = qnanf();
    float pointing = qnanf();
    float path_break = qnanf();
    float score = qnanf();
    Vec2 activity_vertex;
    float activity_residual = qnanf();
};

struct EdgeResult {
    int plane = 0;
    int src_component = -1;
    int dst_component = -1;
    int kind = 0;
    Vec2 src_pos;
    Vec2 dst_pos;
    float gap = qnanf();
    float emptiness = qnanf();
    float visibility = qnanf();
    float pointing = qnanf();
    float path_break = qnanf();
    float score = qnanf();
};

struct Result {
    PlaneInput input;
    std::vector<HitResult> hits;
    std::vector<ComponentResult> components;
    std::vector<EdgeResult> edges;
};

class GraphEngine {
public:
    explicit GraphEngine(Params params = Params {});

    void set_params(const Params& params);
    const Params& params() const;

    Result run(const PlaneInput& input) const;

private:
    struct LocalGeom {
        Vec2 dir {1.0f, 0.0f};
        float anis = 0.0f;
        float scale = 1.0f;
    };

    struct HitNode {
        std::size_t hit_key = 0u;
        int slice_index = -1;
        int cryo = 0;
        int tpc = 0;
        int plane = 0;
        int wire = 0;
        float tick = 0.0f;
        float wire_cm = 0.0f;
        float drift_cm = 0.0f;
        float charge = 0.0f;
        float peakamp = 0.0f;
        float rms = 0.0f;
        Vec2 pos;
        LocalGeom geom;
        int component_local = -1;
        int activity_id = -1;
        bool is_boundary = false;
        bool is_start = false;
    };

    struct CandidateEdge {
        int i = -1;
        int j = -1;
        float dist = 0.0f;
        float support = 0.0f;
        float cost = 0.0f;
    };

    struct CompGraphEdge {
        int a = -1;
        int b = -1;
        float break_weight = 0.0f;
        float min_cost = 0.0f;
    };

    struct ParentCandidate {
        bool valid = false;
        bool visible = false;
        int parent_local = -1;
        int parent_anchor_hit = -1;
        Vec2 parent_anchor_pos;
        float gap = qnanf();
        float emptiness = qnanf();
        float visibility = qnanf();
        float pointing = qnanf();
        float path_break = qnanf();
        float score = qnanf();
    };

    struct ComponentSummary {
        int local_id = -1;
        int activity_id = -1;
        int label = kPrompt;
        int cryo = 0;
        int tpc = 0;
        int plane = 0;
        std::vector<int> hit_indices;
        std::vector<int> boundary_hits;
        int start_hit = -1;
        Vec2 start_pos;
        Vec2 start_dir {1.0f, 0.0f};
        float start_asym = 0.0f;
        float qsum = 0.0f;
        float bbox_min_wire = 0.0f;
        float bbox_max_wire = 0.0f;
        float bbox_min_drift = 0.0f;
        float bbox_max_drift = 0.0f;
        int parent_local = -1;
        Vec2 parent_anchor_pos;
        float gap = qnanf();
        float emptiness = qnanf();
        float visibility = qnanf();
        float pointing = qnanf();
        float path_break = qnanf();
        float score = qnanf();
        Vec2 activity_vertex;
        float activity_residual = qnanf();
    };

    class UnionFind {
    public:
        explicit UnionFind(std::size_t n = 0u)
            : fParent(n)
            , fRank(n, 0u) {
            std::iota(fParent.begin(), fParent.end(), 0u);
        }

        std::size_t find(std::size_t x) {
            if (fParent[x] != x) fParent[x] = find(fParent[x]);
            return fParent[x];
        }

        void unite(std::size_t a, std::size_t b) {
            a = find(a);
            b = find(b);
            if (a == b) return;
            if (fRank[a] < fRank[b]) std::swap(a, b);
            fParent[b] = a;
            if (fRank[a] == fRank[b]) ++fRank[a];
        }

    private:
        std::vector<std::size_t> fParent;
        std::vector<std::size_t> fRank;
    };

    class SpatialGrid {
    public:
        SpatialGrid() = default;
        SpatialGrid(const std::vector<HitNode>* hits, float cell_size);

        void reset(const std::vector<HitNode>* hits, float cell_size);
        std::vector<int> queryRadius(const Vec2& p, float radius) const;
        std::vector<int> queryAABB(float xmin, float xmax, float ymin, float ymax) const;
        bool hasOccupancy(const Vec2& p,
                          float radius,
                          int exclude_hit_a,
                          int exclude_hit_b,
                          int exclude_comp_a,
                          int exclude_comp_b) const;

    private:
        const std::vector<HitNode>* fHits = nullptr;
        float fCellSize = 1.0f;

        struct CellKey {
            int ix = 0;
            int iy = 0;

            bool operator==(const CellKey& other) const {
                return ix == other.ix && iy == other.iy;
            }
        };

        struct CellKeyHash {
            std::size_t operator()(const CellKey& k) const {
                const std::uint64_t a = static_cast<std::uint32_t>(k.ix);
                const std::uint64_t b = static_cast<std::uint32_t>(k.iy);
                return static_cast<std::size_t>((a << 32) ^ b);
            }
        };

        std::unordered_map<CellKey, std::vector<int>, CellKeyHash> fBins;

        CellKey cell(const Vec2& p) const;
        void build();

        static Vec2 subtract(const Vec2& a, const Vec2& b) { return Vec2 {a.x - b.x, a.y - b.y}; }
        static float norm2(const Vec2& v) { return v.x * v.x + v.y * v.y; }
    };

    Params fParams;

    static Vec2 add(const Vec2& a, const Vec2& b) { return Vec2 {a.x + b.x, a.y + b.y}; }
    static Vec2 subtract(const Vec2& a, const Vec2& b) { return Vec2 {a.x - b.x, a.y - b.y}; }
    static Vec2 scale(const Vec2& a, float s) { return Vec2 {a.x * s, a.y * s}; }
    static float dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
    static float norm2(const Vec2& a) { return dot(a, a); }
    static float norm(const Vec2& a) { return std::sqrt(norm2(a)); }

    static Vec2 normalize(const Vec2& a) {
        const float n = norm(a);
        if (n < kEps) return Vec2 {1.0f, 0.0f};
        return scale(a, 1.0f / n);
    }

    static float acute_axis_angle(const Vec2& a, const Vec2& b) {
        const float na = norm(a);
        const float nb = norm(b);
        if (na < kEps || nb < kEps) return 0.5f * kPi;
        const float c = std::fabs(dot(a, b) / (na * nb));
        return std::acos(clamp_value(c, -1.0f, 1.0f));
    }

    static float point_to_bbox_distance(const Vec2& p,
                                        float xmin,
                                        float xmax,
                                        float ymin,
                                        float ymax) {
        const float dx = (p.x < xmin) ? (xmin - p.x) : ((p.x > xmax) ? (p.x - xmax) : 0.0f);
        const float dy = (p.y < ymin) ? (ymin - p.y) : ((p.y > ymax) ? (p.y - ymax) : 0.0f);
        return std::sqrt(dx * dx + dy * dy);
    }

    static float line_projection_t(const Vec2& p, const Vec2& a, const Vec2& b) {
        const Vec2 ab = subtract(b, a);
        const float ab2 = norm2(ab);
        if (ab2 < kEps) return 0.0f;
        return dot(subtract(p, a), ab) / ab2;
    }

    static float point_line_distance(const Vec2& p,
                                     const Vec2& a,
                                     const Vec2& b,
                                     float& t) {
        t = line_projection_t(p, a, b);
        const Vec2 proj = add(a, scale(subtract(b, a), t));
        return norm(subtract(p, proj));
    }

    static float merged_interval_length(std::vector<std::pair<float, float>> intervals,
                                        float lo,
                                        float hi);

    static bool clip_segment_to_box(const Vec2& a,
                                    const Vec2& b,
                                    float xmin,
                                    float xmax,
                                    float ymin,
                                    float ymax,
                                    float& t0,
                                    float& t1);

    static bool solve2x2(float a00,
                         float a01,
                         float a10,
                         float a11,
                         float b0,
                         float b1,
                         Vec2& x);

    DerivedScales make_scales(const PlaneGeometry& geometry) const;

    LocalGeom compute_local_geom_for_hit(int i,
                                         const std::vector<HitNode>& hits,
                                         const SpatialGrid& grid,
                                         const DerivedScales& scales) const;

    void compute_local_geometry(std::vector<HitNode>& hits,
                                const SpatialGrid& grid,
                                const DerivedScales& scales) const;

    float sample_support_fraction(const Vec2& a,
                                  const Vec2& b,
                                  const SpatialGrid& grid,
                                  const std::vector<HitNode>& hits,
                                  float radius,
                                  int exclude_hit_a,
                                  int exclude_hit_b) const;

    std::vector<CandidateEdge> build_candidate_edges(const std::vector<HitNode>& hits,
                                                     const SpatialGrid& grid,
                                                     const DerivedScales& scales) const;

    float edge_cost(const CandidateEdge& edge,
                    const std::vector<HitNode>& hits,
                    const DerivedScales& scales) const;

    std::vector<ComponentSummary> build_components(std::vector<HitNode>& hits,
                                                   const std::vector<CandidateEdge>& edges) const;

    Vec2 component_axis(const ComponentSummary& comp,
                        const std::vector<HitNode>& hits) const;

    Vec2 local_axis_at_hit(int hit_idx,
                           const ComponentSummary& comp,
                           const std::vector<HitNode>& hits,
                           const DerivedScales& scales) const;

    float start_asymmetry(int center_hit_idx,
                          const Vec2& axis,
                          const ComponentSummary& comp,
                          const std::vector<HitNode>& hits,
                          const DerivedScales& scales) const;

    void choose_component_start(ComponentSummary& comp,
                                std::vector<HitNode>& hits,
                                const DerivedScales& scales) const;

    std::vector<CompGraphEdge> build_component_graph(
        const std::vector<HitNode>& hits,
        const std::vector<CandidateEdge>& edges,
        const std::vector<ComponentSummary>& comps,
        std::vector<std::vector<std::pair<int, float>>>& adjacency) const;

    std::vector<float> dijkstra_component_graph(
        const std::vector<std::vector<std::pair<int, float>>>& adjacency,
        int source) const;

    bool inside_plane_bounds(const Vec2& p, const DerivedScales& scales) const;

    void corridor_metrics(const Vec2& a,
                          const Vec2& b,
                          const SpatialGrid& grid,
                          const std::vector<HitNode>& hits,
                          const DerivedScales& scales,
                          int exclude_comp_a,
                          int exclude_comp_b,
                          float& emptiness,
                          float& visibility) const;

    ParentCandidate best_parent_candidate(const ComponentSummary& child,
                                          const std::vector<ComponentSummary>& comps,
                                          const std::vector<HitNode>& hits,
                                          const SpatialGrid& grid,
                                          const DerivedScales& scales,
                                          const std::vector<float>& comp_break_dist) const;

    void classify_components(std::vector<ComponentSummary>& comps,
                             const std::vector<HitNode>& hits,
                             const SpatialGrid& grid,
                             const DerivedScales& scales,
                             const std::vector<std::vector<std::pair<int, float>>>& adjacency) const;

    bool fit_common_secondary_vertex(const std::vector<ComponentSummary*>& group,
                                     Vec2& vertex,
                                     float& mean_residual) const;

    void merge_detached_activities(std::vector<ComponentSummary>& comps) const;
};

} // namespace cg
} // namespace analysis

namespace analysis {
namespace cg {

inline GraphEngine::GraphEngine(Params params)
    : fParams(std::move(params)) {}

inline void GraphEngine::set_params(const Params& params) {
    fParams = params;
}

inline const Params& GraphEngine::params() const {
    return fParams;
}

inline GraphEngine::SpatialGrid::SpatialGrid(const std::vector<HitNode>* hits,
                                                 float cell_size)
    : fHits(hits)
    , fCellSize(std::max(cell_size, 0.05f)) {
    build();
}

inline void GraphEngine::SpatialGrid::reset(const std::vector<HitNode>* hits,
                                                float cell_size) {
    fHits = hits;
    fCellSize = std::max(cell_size, 0.05f);
    build();
}

inline std::vector<int> GraphEngine::SpatialGrid::queryRadius(const Vec2& p,
                                                                  float radius) const {
    std::vector<int> out;
    if (!fHits || fBins.empty()) return out;

    const int reach = static_cast<int>(std::ceil(radius / fCellSize));
    const CellKey c = cell(p);

    for (int dx = -reach; dx <= reach; ++dx) {
        for (int dy = -reach; dy <= reach; ++dy) {
            const CellKey key {c.ix + dx, c.iy + dy};
            const auto iter = fBins.find(key);
            if (iter == fBins.end()) continue;

            for (const int idx : iter->second) {
                const Vec2 d = subtract((*fHits)[idx].pos, p);
                if (norm2(d) <= radius * radius) out.push_back(idx);
            }
        }
    }
    return out;
}

inline std::vector<int> GraphEngine::SpatialGrid::queryAABB(float xmin,
                                                                float xmax,
                                                                float ymin,
                                                                float ymax) const {
    std::vector<int> out;
    if (!fHits || fBins.empty()) return out;

    const CellKey cmin {
        static_cast<int>(std::floor(std::min(xmin, xmax) / fCellSize)),
        static_cast<int>(std::floor(std::min(ymin, ymax) / fCellSize))
    };
    const CellKey cmax {
        static_cast<int>(std::floor(std::max(xmin, xmax) / fCellSize)),
        static_cast<int>(std::floor(std::max(ymin, ymax) / fCellSize))
    };

    for (int ix = cmin.ix; ix <= cmax.ix; ++ix) {
        for (int iy = cmin.iy; iy <= cmax.iy; ++iy) {
            const CellKey key {ix, iy};
            const auto iter = fBins.find(key);
            if (iter == fBins.end()) continue;
            out.insert(out.end(), iter->second.begin(), iter->second.end());
        }
    }
    return out;
}

inline bool GraphEngine::SpatialGrid::hasOccupancy(const Vec2& p,
                                                       float radius,
                                                       int exclude_hit_a,
                                                       int exclude_hit_b,
                                                       int exclude_comp_a,
                                                       int exclude_comp_b) const {
    if (!fHits || fBins.empty()) return false;

    const int reach = static_cast<int>(std::ceil(radius / fCellSize));
    const CellKey c = cell(p);

    for (int dx = -reach; dx <= reach; ++dx) {
        for (int dy = -reach; dy <= reach; ++dy) {
            const CellKey key {c.ix + dx, c.iy + dy};
            const auto iter = fBins.find(key);
            if (iter == fBins.end()) continue;

            for (const int idx : iter->second) {
                if (idx == exclude_hit_a || idx == exclude_hit_b) continue;
                const auto& hit = (*fHits)[idx];
                if (hit.component_local == exclude_comp_a || hit.component_local == exclude_comp_b) continue;

                const Vec2 d = subtract(hit.pos, p);
                if (norm2(d) <= radius * radius) return true;
            }
        }
    }
    return false;
}

inline GraphEngine::SpatialGrid::CellKey
GraphEngine::SpatialGrid::cell(const Vec2& p) const {
    return CellKey {
        static_cast<int>(std::floor(p.x / fCellSize)),
        static_cast<int>(std::floor(p.y / fCellSize))
    };
}

inline void GraphEngine::SpatialGrid::build() {
    fBins.clear();
    if (!fHits) return;
    for (std::size_t i = 0u; i < fHits->size(); ++i) {
        fBins[cell((*fHits)[i].pos)].push_back(static_cast<int>(i));
    }
}

inline float GraphEngine::merged_interval_length(std::vector<std::pair<float, float>> intervals,
                                                     float lo,
                                                     float hi) {
    if (!(hi > lo) || intervals.empty()) return 0.0f;

    for (auto& interval : intervals) {
        interval.first = clamp_value(interval.first, lo, hi);
        interval.second = clamp_value(interval.second, lo, hi);
        if (interval.second < interval.first) std::swap(interval.first, interval.second);
    }

    std::sort(intervals.begin(), intervals.end(),
              [](const std::pair<float, float>& a, const std::pair<float, float>& b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });

    float covered = 0.0f;
    float cur_lo = intervals.front().first;
    float cur_hi = intervals.front().second;

    for (std::size_t i = 1u; i < intervals.size(); ++i) {
        if (intervals[i].second <= cur_lo) continue;
        if (intervals[i].first > cur_hi) {
            covered += std::max(0.0f, cur_hi - cur_lo);
            cur_lo = intervals[i].first;
            cur_hi = intervals[i].second;
        } else {
            cur_hi = std::max(cur_hi, intervals[i].second);
        }
    }

    covered += std::max(0.0f, cur_hi - cur_lo);
    return clamp_value(covered, 0.0f, std::max(0.0f, hi - lo));
}

inline bool GraphEngine::clip_segment_to_box(const Vec2& a,
                                                 const Vec2& b,
                                                 float xmin,
                                                 float xmax,
                                                 float ymin,
                                                 float ymax,
                                                 float& t0,
                                                 float& t1) {
    t0 = 0.0f;
    t1 = 1.0f;

    const float dx = b.x - a.x;
    const float dy = b.y - a.y;

    auto clip = [&](float p, float q) -> bool {
        if (std::fabs(p) < 1.0e-8f) return (q >= 0.0f);

        const float r = q / p;
        if (p < 0.0f) {
            if (r > t1) return false;
            if (r > t0) t0 = r;
        } else {
            if (r < t0) return false;
            if (r < t1) t1 = r;
        }
        return true;
    };

    if (!clip(-dx, a.x - xmin)) return false;
    if (!clip(dx, xmax - a.x)) return false;
    if (!clip(-dy, a.y - ymin)) return false;
    if (!clip(dy, ymax - a.y)) return false;

    return (t1 >= t0);
}

inline bool GraphEngine::solve2x2(float a00,
                                      float a01,
                                      float a10,
                                      float a11,
                                      float b0,
                                      float b1,
                                      Vec2& x) {
    const float det = a00 * a11 - a01 * a10;
    if (std::fabs(det) < 1.0e-8f) return false;
    x.x = (b0 * a11 - a01 * b1) / det;
    x.y = (-b0 * a10 + a00 * b1) / det;
    return true;
}

inline DerivedScales GraphEngine::make_scales(const PlaneGeometry& geometry) const {
    DerivedScales out;
    out.s = std::max(geometry.wire_pitch, geometry.drift_cm_per_tick);

    out.cell_size = std::max(fParams.cell_size_scale * out.s, 0.05f);
    out.candidate_radius = std::max(fParams.candidate_radius_scale * out.s, 0.05f);
    out.geom_radius = std::max(fParams.geom_radius_scale * out.s, 0.05f);
    out.capsule_radius = std::max(fParams.capsule_radius_scale * out.s, 0.05f);
    out.tube_radius = std::max(fParams.tube_radius_scale * out.s, 0.05f);
    out.start_radius = std::max(fParams.start_radius_scale * out.s, 0.10f);
    out.min_gap = std::max(fParams.min_gap_scale * out.s, 0.10f);
    out.max_gap = std::max(fParams.max_gap_scale * out.s, out.min_gap + out.s);
    out.parent_search_radius = std::max(fParams.parent_search_radius_scale * out.s, out.max_gap);

    out.wire_min = geometry.wire_min;
    out.wire_max = geometry.wire_max;
    out.drift_min = geometry.drift_min;
    out.drift_max = geometry.drift_max;
    out.dead_wire_intervals = geometry.dead_wire_intervals;
    return out;
}

inline GraphEngine::LocalGeom GraphEngine::compute_local_geom_for_hit(
    int i,
    const std::vector<HitNode>& hits,
    const SpatialGrid& grid,
    const DerivedScales& scales) const {

    LocalGeom out;
    out.dir = Vec2 {1.0f, 0.0f};
    out.anis = 0.0f;
    out.scale = scales.s;

    auto neighbors = grid.queryRadius(hits[i].pos, scales.geom_radius);
    neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), i), neighbors.end());

    std::sort(neighbors.begin(), neighbors.end(),
              [&](const int a, const int b) {
                  return norm2(subtract(hits[a].pos, hits[i].pos)) <
                         norm2(subtract(hits[b].pos, hits[i].pos));
              });

    if (neighbors.size() > fParams.max_geom_neighbors) neighbors.resize(fParams.max_geom_neighbors);

    if (neighbors.empty()) return out;

    if (neighbors.size() == 1u) {
        out.dir = normalize(subtract(hits[neighbors.front()].pos, hits[i].pos));
        out.anis = 1.0f;
        out.scale = std::max(norm(subtract(hits[neighbors.front()].pos, hits[i].pos)), 0.5f * scales.s);
        return out;
    }

    float wsum = 0.0f;
    Vec2 mean {0.0f, 0.0f};
    for (const int j : neighbors) {
        const float d = std::max(norm(subtract(hits[j].pos, hits[i].pos)), 0.25f * scales.s);
        const float w = (std::max(hits[j].charge, 0.0f) + 1.0f) / d;
        mean = add(mean, scale(hits[j].pos, w));
        wsum += w;
    }
    if (wsum < kEps) return out;
    mean = scale(mean, 1.0f / wsum);

    float cxx = 0.0f;
    float cxy = 0.0f;
    float cyy = 0.0f;
    for (const int j : neighbors) {
        const float d = std::max(norm(subtract(hits[j].pos, hits[i].pos)), 0.25f * scales.s);
        const float w = (std::max(hits[j].charge, 0.0f) + 1.0f) / d;
        const Vec2 x = subtract(hits[j].pos, mean);
        cxx += w * x.x * x.x;
        cxy += w * x.x * x.y;
        cyy += w * x.y * x.y;
    }
    cxx /= wsum;
    cxy /= wsum;
    cyy /= wsum;

    const float trace = cxx + cyy;
    const float disc = std::sqrt(std::max((cxx - cyy) * (cxx - cyy) + 4.0f * cxy * cxy, 0.0f));
    const float l1 = 0.5f * (trace + disc);
    const float l2 = 0.5f * (trace - disc);

    Vec2 dir;
    if (std::fabs(cxy) > 1.0e-8f) {
        dir = Vec2 {l1 - cyy, cxy};
    } else {
        dir = (cxx >= cyy) ? Vec2 {1.0f, 0.0f} : Vec2 {0.0f, 1.0f};
    }

    out.dir = normalize(dir);
    out.anis = (l1 + l2 > kEps) ? clamp_value((l1 - l2) / (l1 + l2 + kEps), 0.0f, 1.0f) : 0.0f;
    out.scale = std::max(std::sqrt(std::max(l1, 0.0f)), 0.5f * scales.s);

    return out;
}

inline void GraphEngine::compute_local_geometry(std::vector<HitNode>& hits,
                                                    const SpatialGrid& grid,
                                                    const DerivedScales& scales) const {
    for (std::size_t i = 0u; i < hits.size(); ++i) {
        hits[i].geom = this->compute_local_geom_for_hit(static_cast<int>(i), hits, grid, scales);
    }
}

inline float GraphEngine::sample_support_fraction(const Vec2& a,
                                                      const Vec2& b,
                                                      const SpatialGrid& grid,
                                                      const std::vector<HitNode>& hits,
                                                      float radius,
                                                      int exclude_hit_a,
                                                      int exclude_hit_b) const {
    const Vec2 ab = subtract(b, a);
    const float L = norm(ab);
    if (L < kEps) return 1.0f;

    const float inner_lo = clamp_value(fParams.capsule_interior_frac, 0.0f, 0.49f);
    const float inner_hi = 1.0f - inner_lo;
    if (!(inner_hi > inner_lo)) return 1.0f;

    const float xmin = std::min(a.x, b.x) - radius;
    const float xmax = std::max(a.x, b.x) + radius;
    const float ymin = std::min(a.y, b.y) - radius;
    const float ymax = std::max(a.y, b.y) + radius;

    const auto candidates = grid.queryAABB(xmin, xmax, ymin, ymax);
    std::vector<std::pair<float, float>> intervals;
    intervals.reserve(candidates.size());

    for (const int idx : candidates) {
        if (idx == exclude_hit_a || idx == exclude_hit_b) continue;

        float t_center = 0.0f;
        const float d_perp = point_line_distance(hits[idx].pos, a, b, t_center);
        if (d_perp > radius) continue;

        const float dt = std::sqrt(std::max(radius * radius - d_perp * d_perp, 0.0f)) / L;
        const float t0 = t_center - dt;
        const float t1 = t_center + dt;
        if (t1 <= inner_lo || t0 >= inner_hi) continue;

        intervals.emplace_back(std::max(t0, inner_lo), std::min(t1, inner_hi));
    }

    const float covered = merged_interval_length(intervals, inner_lo, inner_hi);
    return clamp_value(covered / std::max(inner_hi - inner_lo, kEps), 0.0f, 1.0f);
}

inline std::vector<GraphEngine::CandidateEdge>
GraphEngine::build_candidate_edges(const std::vector<HitNode>& hits,
                                              const SpatialGrid& grid,
                                              const DerivedScales& scales) const {
    std::vector<CandidateEdge> edges;
    edges.reserve(hits.size() * 6u);

    for (std::size_t i = 0u; i < hits.size(); ++i) {
        auto neighbors = grid.queryRadius(hits[i].pos, scales.candidate_radius);

        std::vector<std::pair<float, int>> ordered;
        ordered.reserve(neighbors.size());

        for (const int j : neighbors) {
            if (static_cast<int>(i) >= j) continue;
            const float d = norm(subtract(hits[j].pos, hits[i].pos));
            ordered.emplace_back(d, j);
        }

        std::sort(ordered.begin(), ordered.end(),
                  [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
                      return a.first < b.first;
                  });

        if (ordered.size() > fParams.max_candidate_neighbors) {
            ordered.resize(fParams.max_candidate_neighbors);
        }

        for (const auto& dj : ordered) {
            CandidateEdge edge;
            edge.i = static_cast<int>(i);
            edge.j = dj.second;
            edge.dist = dj.first;
            edge.support = this->sample_support_fraction(hits[edge.i].pos,
                                                         hits[edge.j].pos,
                                                         grid,
                                                         hits,
                                                         scales.capsule_radius,
                                                         edge.i,
                                                         edge.j);
            edge.cost = this->edge_cost(edge, hits, scales);
            edges.push_back(edge);
        }
    }

    return edges;
}

inline float GraphEngine::edge_cost(const CandidateEdge& edge,
                                        const std::vector<HitNode>& hits,
                                        const DerivedScales& scales) const {
    const auto& hi = hits[edge.i];
    const auto& hj = hits[edge.j];

    const Vec2 r = subtract(hj.pos, hi.pos);
    const float d = std::max(norm(r), 0.01f);
    const Vec2 rhat = scale(r, 1.0f / d);

    const float theta = 0.5f * (
        acute_axis_angle(hi.geom.dir, rhat) +
        acute_axis_angle(hj.geom.dir, rhat)
    );

    const float scale_sum = std::max(hi.geom.scale + hj.geom.scale, 0.5f * scales.s);
    const float dist_term = d / scale_sum;
    const float angle_term = hi.geom.anis * hj.geom.anis * theta;
    const float dq_term = std::fabs(std::log((std::max(hi.charge, 0.0f) + 1.0e-3f) /
                                             (std::max(hj.charge, 0.0f) + 1.0e-3f)));

    return fParams.alpha_dist * dist_term
         + fParams.beta_angle * angle_term
         + fParams.gamma_support * (1.0f - edge.support)
         + fParams.delta_charge * dq_term;
}

inline std::vector<GraphEngine::ComponentSummary>
GraphEngine::build_components(std::vector<HitNode>& hits,
                                         const std::vector<CandidateEdge>& edges) const {
    UnionFind uf(hits.size());
    for (const auto& edge : edges) {
        if (edge.cost < fParams.tau_ch) {
            uf.unite(static_cast<std::size_t>(edge.i), static_cast<std::size_t>(edge.j));
        }
    }

    std::map<std::size_t, int> root_to_local;
    std::vector<ComponentSummary> comps;

    for (std::size_t i = 0u; i < hits.size(); ++i) {
        const std::size_t root = uf.find(i);
        auto iter = root_to_local.find(root);
        if (iter == root_to_local.end()) {
            const int new_id = static_cast<int>(comps.size());
            root_to_local[root] = new_id;
            ComponentSummary comp;
            comp.local_id = new_id;
            comp.activity_id = new_id;
            comp.cryo = hits[i].cryo;
            comp.tpc = hits[i].tpc;
            comp.plane = hits[i].plane;
            comp.bbox_min_wire = hits[i].wire_cm;
            comp.bbox_max_wire = hits[i].wire_cm;
            comp.bbox_min_drift = hits[i].drift_cm;
            comp.bbox_max_drift = hits[i].drift_cm;
            comps.push_back(comp);
            iter = root_to_local.find(root);
        }

        ComponentSummary& comp = comps[iter->second];
        comp.hit_indices.push_back(static_cast<int>(i));
        comp.qsum += std::max(hits[i].charge, 0.0f);
        comp.bbox_min_wire = std::min(comp.bbox_min_wire, hits[i].wire_cm);
        comp.bbox_max_wire = std::max(comp.bbox_max_wire, hits[i].wire_cm);
        comp.bbox_min_drift = std::min(comp.bbox_min_drift, hits[i].drift_cm);
        comp.bbox_max_drift = std::max(comp.bbox_max_drift, hits[i].drift_cm);

        hits[i].component_local = comp.local_id;
    }

    return comps;
}

inline Vec2 GraphEngine::component_axis(const ComponentSummary& comp,
                                            const std::vector<HitNode>& hits) const {
    if (comp.hit_indices.size() < 2u) return Vec2 {1.0f, 0.0f};

    Vec2 mean {0.0f, 0.0f};
    float qsum = 0.0f;
    for (const int idx : comp.hit_indices) {
        const float w = std::max(hits[idx].charge, 0.0f) + 1.0f;
        mean = add(mean, scale(hits[idx].pos, w));
        qsum += w;
    }
    mean = scale(mean, 1.0f / std::max(qsum, kEps));

    float cxx = 0.0f;
    float cxy = 0.0f;
    float cyy = 0.0f;
    for (const int idx : comp.hit_indices) {
        const float w = std::max(hits[idx].charge, 0.0f) + 1.0f;
        const Vec2 x = subtract(hits[idx].pos, mean);
        cxx += w * x.x * x.x;
        cxy += w * x.x * x.y;
        cyy += w * x.y * x.y;
    }
    cxx /= std::max(qsum, kEps);
    cxy /= std::max(qsum, kEps);
    cyy /= std::max(qsum, kEps);

    const float trace = cxx + cyy;
    const float disc = std::sqrt(std::max((cxx - cyy) * (cxx - cyy) + 4.0f * cxy * cxy, 0.0f));
    const float l1 = 0.5f * (trace + disc);

    Vec2 dir;
    if (std::fabs(cxy) > 1.0e-8f) dir = Vec2 {l1 - cyy, cxy};
    else dir = (cxx >= cyy) ? Vec2 {1.0f, 0.0f} : Vec2 {0.0f, 1.0f};

    return normalize(dir);
}

inline Vec2 GraphEngine::local_axis_at_hit(int hit_idx,
                                               const ComponentSummary& comp,
                                               const std::vector<HitNode>& hits,
                                               const DerivedScales& scales) const {
    std::vector<int> near_hits;
    near_hits.reserve(comp.hit_indices.size());

    for (const int idx : comp.hit_indices) {
        if (idx == hit_idx) continue;
        if (norm(subtract(hits[idx].pos, hits[hit_idx].pos)) <= scales.start_radius) {
            near_hits.push_back(idx);
        }
    }

    std::sort(near_hits.begin(), near_hits.end(),
              [&](const int a, const int b) {
                  return norm2(subtract(hits[a].pos, hits[hit_idx].pos)) <
                         norm2(subtract(hits[b].pos, hits[hit_idx].pos));
              });

    if (near_hits.size() > fParams.max_geom_neighbors) near_hits.resize(fParams.max_geom_neighbors);

    if (near_hits.size() < 2u) return this->component_axis(comp, hits);

    Vec2 mean {0.0f, 0.0f};
    float wsum = 0.0f;
    for (const int idx : near_hits) {
        const float d = std::max(norm(subtract(hits[idx].pos, hits[hit_idx].pos)), 0.25f * scales.s);
        const float w = (std::max(hits[idx].charge, 0.0f) + 1.0f) / d;
        mean = add(mean, scale(hits[idx].pos, w));
        wsum += w;
    }
    mean = scale(mean, 1.0f / std::max(wsum, kEps));

    float cxx = 0.0f;
    float cxy = 0.0f;
    float cyy = 0.0f;
    for (const int idx : near_hits) {
        const float d = std::max(norm(subtract(hits[idx].pos, hits[hit_idx].pos)), 0.25f * scales.s);
        const float w = (std::max(hits[idx].charge, 0.0f) + 1.0f) / d;
        const Vec2 x = subtract(hits[idx].pos, mean);
        cxx += w * x.x * x.x;
        cxy += w * x.x * x.y;
        cyy += w * x.y * x.y;
    }
    cxx /= std::max(wsum, kEps);
    cxy /= std::max(wsum, kEps);
    cyy /= std::max(wsum, kEps);

    const float trace = cxx + cyy;
    const float disc = std::sqrt(std::max((cxx - cyy) * (cxx - cyy) + 4.0f * cxy * cxy, 0.0f));
    const float l1 = 0.5f * (trace + disc);

    Vec2 dir;
    if (std::fabs(cxy) > 1.0e-8f) dir = Vec2 {l1 - cyy, cxy};
    else dir = (cxx >= cyy) ? Vec2 {1.0f, 0.0f} : Vec2 {0.0f, 1.0f};

    return normalize(dir);
}

inline float GraphEngine::start_asymmetry(int center_hit_idx,
                                              const Vec2& axis,
                                              const ComponentSummary& comp,
                                              const std::vector<HitNode>& hits,
                                              const DerivedScales& scales) const {
    const float cut = std::cos(fParams.cone_angle_deg * kPi / 180.0f);
    float forward = 0.0f;
    float backward = 0.0f;

    for (const int idx : comp.hit_indices) {
        if (idx == center_hit_idx) continue;

        const Vec2 rel = subtract(hits[idx].pos, hits[center_hit_idx].pos);
        const float r = norm(rel);
        if (r < kEps || r > scales.start_radius) continue;

        const float c = dot(rel, axis) / r;
        if (c >= cut) forward += std::max(hits[idx].charge, 0.0f);
        if (c <= -cut) backward += std::max(hits[idx].charge, 0.0f);
    }

    return forward / (forward + fParams.start_kappa * backward + kEps);
}

inline void GraphEngine::choose_component_start(ComponentSummary& comp,
                                                    std::vector<HitNode>& hits,
                                                    const DerivedScales& scales) const {
    const Vec2 axis = this->component_axis(comp, hits);

    std::vector<std::pair<float, int>> proj;
    proj.reserve(comp.hit_indices.size());
    for (const int idx : comp.hit_indices) {
        proj.emplace_back(dot(hits[idx].pos, axis), idx);
    }
    std::sort(proj.begin(), proj.end(),
              [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
                  return a.first < b.first;
              });

    std::unordered_set<int> candidate_set;
    const unsigned nside = std::max(1u, fParams.boundary_candidates_per_side);
    for (unsigned i = 0u; i < std::min<std::size_t>(nside, proj.size()); ++i) {
        candidate_set.insert(proj[i].second);
        candidate_set.insert(proj[proj.size() - 1u - i].second);
    }

    comp.boundary_hits.assign(candidate_set.begin(), candidate_set.end());

    float best_asym = -1.0f;
    int best_idx = -1;
    Vec2 best_dir {1.0f, 0.0f};

    for (const int idx : comp.boundary_hits) {
        hits[idx].is_boundary = true;

        Vec2 loc_axis = this->local_axis_at_hit(idx, comp, hits, scales);
        float a1 = this->start_asymmetry(idx, loc_axis, comp, hits, scales);
        float a2 = this->start_asymmetry(idx, scale(loc_axis, -1.0f), comp, hits, scales);

        if (a2 > a1) {
            a1 = a2;
            loc_axis = scale(loc_axis, -1.0f);
        }

        if (a1 > best_asym) {
            best_asym = a1;
            best_idx = idx;
            best_dir = normalize(loc_axis);
        }
    }

    if (best_idx < 0 && !comp.hit_indices.empty()) {
        best_idx = comp.hit_indices.front();
        best_dir = axis;
        best_asym = 0.0f;
    }

    comp.start_hit = best_idx;
    comp.start_pos = (best_idx >= 0) ? hits[best_idx].pos : Vec2 {};
    comp.start_dir = best_dir;
    comp.start_asym = std::max(best_asym, 0.0f);
    comp.activity_vertex = comp.start_pos;
    comp.activity_residual = qnanf();

    if (best_idx >= 0) hits[best_idx].is_start = true;
}

inline std::vector<GraphEngine::CompGraphEdge>
GraphEngine::build_component_graph(
    const std::vector<HitNode>& hits,
    const std::vector<CandidateEdge>& edges,
    const std::vector<ComponentSummary>& comps,
    std::vector<std::vector<std::pair<int, float>>>& adjacency) const {

    struct PairStats {
        float min_break = std::numeric_limits<float>::max();
        float min_cost = std::numeric_limits<float>::max();
    };

    std::unordered_map<std::uint64_t, PairStats> pair_map;

    auto make_key = [](int a, int b) {
        const std::uint32_t x = static_cast<std::uint32_t>(std::min(a, b));
        const std::uint32_t y = static_cast<std::uint32_t>(std::max(a, b));
        return (static_cast<std::uint64_t>(x) << 32) ^ static_cast<std::uint64_t>(y);
    };

    for (const auto& edge : edges) {
        const int ca = hits[edge.i].component_local;
        const int cb = hits[edge.j].component_local;
        if (ca < 0 || cb < 0 || ca == cb) continue;

        PairStats& stats = pair_map[make_key(ca, cb)];
        stats.min_break = std::min(stats.min_break, std::max(edge.cost - fParams.tau_ch, 0.0f));
        stats.min_cost = std::min(stats.min_cost, edge.cost);
    }

    adjacency.assign(comps.size(), {});
    std::vector<CompGraphEdge> out;
    out.reserve(pair_map.size());

    for (const auto& item : pair_map) {
        const int a = static_cast<int>(item.first >> 32);
        const int b = static_cast<int>(item.first & 0xffffffffu);
        const float w = item.second.min_break;

        adjacency[a].emplace_back(b, w);
        adjacency[b].emplace_back(a, w);

        CompGraphEdge rec;
        rec.a = a;
        rec.b = b;
        rec.break_weight = w;
        rec.min_cost = item.second.min_cost;
        out.push_back(rec);
    }

    return out;
}

inline std::vector<float> GraphEngine::dijkstra_component_graph(
    const std::vector<std::vector<std::pair<int, float>>>& adjacency,
    int source) const {

    const float inf = std::numeric_limits<float>::max();
    std::vector<float> dist(adjacency.size(), inf);

    using Node = std::pair<float, int>;
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

    dist[source] = 0.0f;
    pq.emplace(0.0f, source);

    while (!pq.empty()) {
        const auto current = pq.top();
        pq.pop();

        const float d = current.first;
        const int u = current.second;
        if (d > dist[u]) continue;

        for (const auto& edge : adjacency[u]) {
            const int v = edge.first;
            const float nd = d + edge.second;
            if (nd < dist[v]) {
                dist[v] = nd;
                pq.emplace(nd, v);
            }
        }
    }

    return dist;
}

inline bool GraphEngine::inside_plane_bounds(const Vec2& p,
                                                 const DerivedScales& scales) const {
    return (p.x >= scales.wire_min && p.x <= scales.wire_max &&
            p.y >= scales.drift_min && p.y <= scales.drift_max);
}

inline void GraphEngine::corridor_metrics(const Vec2& a,
                                              const Vec2& b,
                                              const SpatialGrid& grid,
                                              const std::vector<HitNode>& hits,
                                              const DerivedScales& scales,
                                              int exclude_comp_a,
                                              int exclude_comp_b,
                                              float& emptiness,
                                              float& visibility) const {
    emptiness = 0.0f;
    visibility = 0.0f;

    const Vec2 ab = subtract(b, a);
    const float L = norm(ab);
    if (L < kEps) {
        emptiness = 0.0f;
        visibility = 1.0f;
        return;
    }

    const float inner_lo = clamp_value(fParams.corridor_end_exclude_frac, 0.0f, 0.49f);
    const float inner_hi = 1.0f - inner_lo;
    if (!(inner_hi > inner_lo)) return;

    float box_lo = 0.0f;
    float box_hi = 1.0f;
    if (!clip_segment_to_box(a,
                             b,
                             scales.wire_min,
                             scales.wire_max,
                             scales.drift_min,
                             scales.drift_max,
                             box_lo,
                             box_hi)) {
        return;
    }

    const float vis_lo = std::max(inner_lo, box_lo);
    const float vis_hi = std::min(inner_hi, box_hi);
    const float interior_len = inner_hi - inner_lo;
    if (!(vis_hi > vis_lo) || !(interior_len > 0.0f)) {
        visibility = 0.0f;
        emptiness = 0.0f;
        return;
    }

    std::vector<std::pair<float, float>> dead_intervals;
    dead_intervals.reserve(scales.dead_wire_intervals.size());

    const float dx = b.x - a.x;
    if (std::fabs(dx) < kEps) {
        for (const auto& dead : scales.dead_wire_intervals) {
            if (a.x >= dead.first && a.x <= dead.second) {
                dead_intervals.emplace_back(vis_lo, vis_hi);
                break;
            }
        }
    } else {
        for (const auto& dead : scales.dead_wire_intervals) {
            float t0 = (dead.first - a.x) / dx;
            float t1 = (dead.second - a.x) / dx;
            if (t1 < t0) std::swap(t0, t1);
            if (t1 <= vis_lo || t0 >= vis_hi) continue;
            dead_intervals.emplace_back(std::max(t0, vis_lo), std::min(t1, vis_hi));
        }
    }

    std::vector<std::pair<float, float>> live_intervals;
    live_intervals.reserve(dead_intervals.size() + 1u);

    if (dead_intervals.empty()) {
        live_intervals.emplace_back(vis_lo, vis_hi);
    } else {
        for (auto& interval : dead_intervals) {
            interval.first = clamp_value(interval.first, vis_lo, vis_hi);
            interval.second = clamp_value(interval.second, vis_lo, vis_hi);
            if (interval.second < interval.first) std::swap(interval.first, interval.second);
        }

        std::sort(dead_intervals.begin(), dead_intervals.end(),
                  [](const std::pair<float, float>& lhs, const std::pair<float, float>& rhs) {
                      if (lhs.first != rhs.first) return lhs.first < rhs.first;
                      return lhs.second < rhs.second;
                  });

        float cursor = vis_lo;
        for (const auto& interval : dead_intervals) {
            if (interval.second <= cursor) continue;
            if (interval.first > cursor) {
                live_intervals.emplace_back(cursor, std::min(interval.first, vis_hi));
            }
            cursor = std::max(cursor, interval.second);
            if (cursor >= vis_hi) break;
        }
        if (cursor < vis_hi) live_intervals.emplace_back(cursor, vis_hi);
    }

    float live_visible_len = 0.0f;
    for (const auto& interval : live_intervals) {
        live_visible_len += std::max(0.0f, interval.second - interval.first);
    }

    visibility = clamp_value(live_visible_len / interior_len, 0.0f, 1.0f);

    const float xmin = std::min(a.x, b.x) - scales.tube_radius;
    const float xmax = std::max(a.x, b.x) + scales.tube_radius;
    const float ymin = std::min(a.y, b.y) - scales.tube_radius;
    const float ymax = std::max(a.y, b.y) + scales.tube_radius;

    const auto candidates = grid.queryAABB(xmin, xmax, ymin, ymax);
    std::vector<std::pair<float, float>> intervals;

    if (!(live_visible_len > 0.0f)) {
        // Entire in-bounds corridor is masked by dead channels:
        // unobservable, but not positively obstructed.
        emptiness = 1.0f;
        return;
    }

    intervals.reserve(candidates.size());

    for (const int idx : candidates) {
        const auto& hit = hits[idx];
        if (hit.component_local == exclude_comp_a || hit.component_local == exclude_comp_b) continue;

        float t_center = 0.0f;
        const float d_perp = point_line_distance(hit.pos, a, b, t_center);
        if (d_perp > scales.tube_radius) continue;

        const float dt = std::sqrt(std::max(scales.tube_radius * scales.tube_radius - d_perp * d_perp, 0.0f)) / L;
        const float t0 = t_center - dt;
        const float t1 = t_center + dt;
        for (const auto& live : live_intervals) {
            if (t1 <= live.first || t0 >= live.second) continue;
            intervals.emplace_back(std::max(t0, live.first), std::min(t1, live.second));
        }
    }

    const float covered = merged_interval_length(intervals, vis_lo, vis_hi);
    emptiness = clamp_value(1.0f - covered / std::max(live_visible_len, kEps), 0.0f, 1.0f);
}

inline GraphEngine::ParentCandidate GraphEngine::best_parent_candidate(
    const ComponentSummary& child,
    const std::vector<ComponentSummary>& comps,
    const std::vector<HitNode>& hits,
    const SpatialGrid& grid,
    const DerivedScales& scales,
    const std::vector<float>& comp_break_dist) const {

    ParentCandidate best_visible;
    ParentCandidate best_unresolved;
    float best_visible_score = -1.0e9f;
    float best_unresolved_score = -1.0e9f;

    for (const auto& parent : comps) {
        if (parent.local_id == child.local_id) continue;

        const float coarse = point_to_bbox_distance(child.start_pos,
                                                    parent.bbox_min_wire,
                                                    parent.bbox_max_wire,
                                                    parent.bbox_min_drift,
                                                    parent.bbox_max_drift);
        if (coarse > scales.parent_search_radius) continue;

        const float raw_path_break = comp_break_dist[parent.local_id];
        const bool have_path = (raw_path_break < std::numeric_limits<float>::max() * 0.5f);
        const float path_break = have_path ? raw_path_break : (5.0f * fParams.min_path_break);

        for (const int anchor_idx : parent.boundary_hits) {
            const Vec2 anchor = hits[anchor_idx].pos;
            const Vec2 bridge = subtract(child.start_pos, anchor);
            const float gap = norm(bridge);

            if (gap < scales.min_gap || gap > scales.max_gap) continue;
            if (path_break < fParams.min_path_break) continue;

            const Vec2 bridge_dir = normalize(bridge);
            const float pointing = dot(child.start_dir, bridge_dir);
            if (pointing < fParams.min_pointing) continue;

            float emptiness = 0.0f;
            float visibility = 0.0f;
            this->corridor_metrics(anchor,
                                   child.start_pos,
                                   grid,
                                   hits,
                                   scales,
                                   child.local_id,
                                   parent.local_id,
                                   emptiness,
                                   visibility);

            const float score = emptiness +
                                pointing +
                                std::min(path_break / std::max(fParams.min_path_break, kEps), 3.0f);

            ParentCandidate cand;
            cand.valid = true;
            cand.visible = (visibility >= fParams.min_visibility);
            cand.parent_local = parent.local_id;
            cand.parent_anchor_hit = anchor_idx;
            cand.parent_anchor_pos = anchor;
            cand.gap = gap;
            cand.emptiness = emptiness;
            cand.visibility = visibility;
            cand.pointing = pointing;
            cand.path_break = path_break;
            cand.score = score;

            if (emptiness >= fParams.min_emptiness && visibility >= fParams.min_visibility) {
                if (score > best_visible_score) {
                    best_visible_score = score;
                    best_visible = cand;
                }
            } else if (emptiness >= fParams.min_emptiness && visibility < fParams.min_visibility) {
                if (score > best_unresolved_score) {
                    best_unresolved_score = score;
                    best_unresolved = cand;
                }
            }
        }
    }

    if (best_visible.valid) return best_visible;
    if (best_unresolved.valid) return best_unresolved;
    return ParentCandidate {};
}

inline void GraphEngine::classify_components(
    std::vector<ComponentSummary>& comps,
    const std::vector<HitNode>& hits,
    const SpatialGrid& grid,
    const DerivedScales& scales,
    const std::vector<std::vector<std::pair<int, float>>>& adjacency) const {

    for (auto& comp : comps) {
        if (comp.hit_indices.size() < fParams.min_component_hits) {
            comp.label = kUnresolved;
            continue;
        }

        const std::vector<float> break_dist = this->dijkstra_component_graph(adjacency, comp.local_id);
        const ParentCandidate parent = this->best_parent_candidate(comp, comps, hits, grid, scales, break_dist);

        if (parent.valid && parent.visible && comp.start_asym >= fParams.min_start_asym) {
            comp.label = kDetached;
        } else if (parent.valid && !parent.visible) {
            comp.label = kUnresolved;
        } else {
            comp.label = kPrompt;
        }

        if (parent.valid) {
            comp.parent_local = parent.parent_local;
            comp.parent_anchor_pos = parent.parent_anchor_pos;
            comp.gap = parent.gap;
            comp.emptiness = parent.emptiness;
            comp.visibility = parent.visibility;
            comp.pointing = parent.pointing;
            comp.path_break = parent.path_break;
            comp.score = parent.score;

            if (comp.start_asym < fParams.min_start_asym && comp.label == kDetached) comp.label = kUnresolved;
            if (comp.start_asym < fParams.min_start_asym && comp.label == kPrompt && parent.valid) {
                comp.label = kUnresolved;
            }
        }
    }
}

inline bool GraphEngine::fit_common_secondary_vertex(const std::vector<ComponentSummary*>& group,
                                                         Vec2& vertex,
                                                         float& mean_residual) const {
    vertex = Vec2 {};
    mean_residual = qnanf();

    if (group.size() < 2u) return false;

    float a00 = 0.0f;
    float a01 = 0.0f;
    float a11 = 0.0f;
    float b0 = 0.0f;
    float b1 = 0.0f;

    for (const auto* comp : group) {
        const Vec2 u = normalize(comp->start_dir);
        const float p00 = 1.0f - u.x * u.x;
        const float p01 = -u.x * u.y;
        const float p11 = 1.0f - u.y * u.y;

        a00 += p00;
        a01 += p01;
        a11 += p11;

        b0 += p00 * comp->start_pos.x + p01 * comp->start_pos.y;
        b1 += p01 * comp->start_pos.x + p11 * comp->start_pos.y;
    }

    if (!solve2x2(a00, a01, a01, a11, b0, b1, vertex)) return false;

    float sum_res2 = 0.0f;
    for (const auto* comp : group) {
        const Vec2 u = normalize(comp->start_dir);
        const Vec2 d = subtract(vertex, comp->start_pos);
        const float proj = dot(d, u);
        if (proj > fParams.max_secondary_vertex_back_projection) return false;

        const Vec2 perp = subtract(d, scale(u, proj));
        const float r = norm(perp);
        if (r > fParams.max_secondary_vertex_residual) return false;
        sum_res2 += r * r;
    }

    mean_residual = std::sqrt(sum_res2 / static_cast<float>(group.size()));
    return true;
}

inline void GraphEngine::merge_detached_activities(std::vector<ComponentSummary>& comps) const {
    for (auto& comp : comps) {
        comp.activity_id = comp.local_id;
        comp.activity_vertex = comp.start_pos;
        comp.activity_residual = qnanf();
    }

    if (!fParams.merge_detached_activities) return;

    std::map<int, std::vector<int>> by_parent;
    for (std::size_t i = 0u; i < comps.size(); ++i) {
        if (comps[i].label == kDetached && comps[i].parent_local >= 0) {
            by_parent[comps[i].parent_local].push_back(static_cast<int>(i));
        }
    }

    for (const auto& item : by_parent) {
        const auto& members = item.second;
        if (members.size() < 2u) continue;

        UnionFind uf(members.size());

        for (std::size_t i = 0u; i < members.size(); ++i) {
            for (std::size_t j = i + 1u; j < members.size(); ++j) {
                std::vector<ComponentSummary*> pair {
                    &comps[members[i]],
                    &comps[members[j]]
                };
                Vec2 vtx;
                float resid = 0.0f;
                if (this->fit_common_secondary_vertex(pair, vtx, resid)) {
                    uf.unite(i, j);
                }
            }
        }

        std::map<std::size_t, std::vector<int>> groups;
        for (std::size_t i = 0u; i < members.size(); ++i) {
            groups[uf.find(i)].push_back(members[i]);
        }

        for (const auto& g : groups) {
            std::vector<ComponentSummary*> comps_ptr;
            comps_ptr.reserve(g.second.size());
            for (const int idx : g.second) comps_ptr.push_back(&comps[idx]);

            Vec2 vtx;
            float resid = 0.0f;
            if (!this->fit_common_secondary_vertex(comps_ptr, vtx, resid)) continue;

            const int activity_id = comps[g.second.front()].local_id;
            for (const int idx : g.second) {
                comps[idx].activity_id = activity_id;
                comps[idx].activity_vertex = vtx;
                comps[idx].activity_residual = resid;
            }
        }
    }
}

inline Result GraphEngine::run(const PlaneInput& input) const {
    Result out;
    out.input = input;

    if (input.hits.empty()) return out;

    const DerivedScales scales = this->make_scales(input.geometry);

    std::vector<HitNode> hits;
    hits.reserve(input.hits.size());
    for (const auto& src : input.hits) {
        HitNode node;
        node.hit_key = src.hit_key;
        node.slice_index = src.slice_index;
        node.cryo = src.cryo;
        node.tpc = src.tpc;
        node.plane = src.plane;
        node.wire = src.wire;
        node.tick = src.tick;
        node.wire_cm = src.wire_cm;
        node.drift_cm = src.drift_cm;
        node.charge = src.charge;
        node.peakamp = src.peakamp;
        node.rms = src.rms;
        node.pos = Vec2 {src.wire_cm, src.drift_cm};
        hits.push_back(node);
    }

    SpatialGrid grid(&hits, scales.cell_size);
    this->compute_local_geometry(hits, grid, scales);
    const auto edges = this->build_candidate_edges(hits, grid, scales);
    auto comps = this->build_components(hits, edges);

    for (auto& comp : comps) {
        this->choose_component_start(comp, hits, scales);
    }

    std::vector<std::vector<std::pair<int, float>>> adjacency;
    const auto comp_graph_edges = this->build_component_graph(hits, edges, comps, adjacency);

    this->classify_components(comps, hits, grid, scales, adjacency);
    this->merge_detached_activities(comps);

    for (auto& hit : hits) {
        if (hit.component_local >= 0 && hit.component_local < static_cast<int>(comps.size())) {
            hit.activity_id = comps[hit.component_local].activity_id;
        }
    }

    out.hits.reserve(hits.size());
    for (const auto& hit : hits) {
        HitResult rec;
        rec.hit_key = hit.hit_key;
        rec.slice_index = hit.slice_index;
        rec.cryo = hit.cryo;
        rec.tpc = hit.tpc;
        rec.plane = hit.plane;
        rec.wire = hit.wire;
        rec.tick = hit.tick;
        rec.wire_cm = hit.wire_cm;
        rec.drift_cm = hit.drift_cm;
        rec.charge = hit.charge;
        rec.peakamp = hit.peakamp;
        rec.rms = hit.rms;
        rec.dir = hit.geom.dir;
        rec.anis = hit.geom.anis;
        rec.local_scale = hit.geom.scale;
        rec.component_id = hit.component_local;
        rec.activity_id = hit.activity_id;
        rec.is_boundary = hit.is_boundary;
        rec.is_start = hit.is_start;
        rec.label = (hit.component_local >= 0 && hit.component_local < static_cast<int>(comps.size()))
                        ? comps[hit.component_local].label
                        : kUnresolved;
        out.hits.push_back(rec);
    }

    out.components.reserve(comps.size());
    for (const auto& comp : comps) {
        ComponentResult rec;
        rec.id = comp.local_id;
        rec.activity_id = comp.activity_id;
        rec.label = comp.label;
        rec.cryo = comp.cryo;
        rec.tpc = comp.tpc;
        rec.plane = comp.plane;
        rec.hit_indices = comp.hit_indices;
        rec.start_hit = comp.start_hit;
        rec.start_pos = comp.start_pos;
        rec.start_dir = comp.start_dir;
        rec.start_asym = comp.start_asym;
        rec.qsum = comp.qsum;
        rec.bbox_min_wire = comp.bbox_min_wire;
        rec.bbox_max_wire = comp.bbox_max_wire;
        rec.bbox_min_drift = comp.bbox_min_drift;
        rec.bbox_max_drift = comp.bbox_max_drift;
        rec.parent_id = comp.parent_local;
        rec.parent_anchor_pos = comp.parent_anchor_pos;
        rec.gap = comp.gap;
        rec.emptiness = comp.emptiness;
        rec.visibility = comp.visibility;
        rec.pointing = comp.pointing;
        rec.path_break = comp.path_break;
        rec.score = comp.score;
        rec.activity_vertex = comp.activity_vertex;
        rec.activity_residual = comp.activity_residual;
        out.components.push_back(rec);
    }

    out.edges.reserve(comp_graph_edges.size() + comps.size());
    for (const auto& comp : comps) {
        if (comp.parent_local < 0) continue;

        EdgeResult edge;
        edge.plane = comp.plane;
        edge.src_component = comp.parent_local;
        edge.dst_component = comp.local_id;
        edge.kind = (comp.label == kDetached) ? 1 : 3;
        edge.src_pos = comp.parent_anchor_pos;
        edge.dst_pos = comp.start_pos;
        edge.gap = comp.gap;
        edge.emptiness = comp.emptiness;
        edge.visibility = comp.visibility;
        edge.pointing = comp.pointing;
        edge.path_break = comp.path_break;
        edge.score = comp.score;
        out.edges.push_back(edge);
    }

    if (fParams.save_all_component_graph_edges) {
        for (const auto& edge : comp_graph_edges) {
            if (edge.a < 0 || edge.b < 0 ||
                edge.a >= static_cast<int>(comps.size()) || edge.b >= static_cast<int>(comps.size())) {
                continue;
            }

            EdgeResult rec;
            rec.plane = comps[edge.a].plane;
            rec.src_component = edge.a;
            rec.dst_component = edge.b;
            rec.kind = 0;
            rec.src_pos = comps[edge.a].start_pos;
            rec.dst_pos = comps[edge.b].start_pos;
            rec.gap = norm(subtract(comps[edge.a].start_pos, comps[edge.b].start_pos));
            rec.emptiness = qnanf();
            rec.visibility = qnanf();
            rec.pointing = qnanf();
            rec.path_break = edge.break_weight;
            rec.score = edge.min_cost;
            out.edges.push_back(rec);
        }
    }

    return out;
}

} // namespace cg
} // namespace analysis

#endif
