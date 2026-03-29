#ifndef ANALYSIS_GRAPHENGINE_H
#define ANALYSIS_GRAPHENGINE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
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

#endif
