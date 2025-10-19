#ifndef DATAPRODUCTS_VERTEX_DISPLACEMENT_H
#define DATAPRODUCTS_VERTEX_DISPLACEMENT_H

#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct VDPeak {
  std::uint16_t col{0};
  std::uint16_t row{0};
  float         score{0.f};
};

struct VDPlaneSummary {
  int           view{0};                  // geo::kU/V/W as int
  std::uint32_t width{0}, height{0};

  // Beam fed to the model for this plane (pixel units before the script’s internal normalization)
  float         beam_col_pix{0.f};
  float         beam_row_pix{0.f};

  // PV + scale
  std::uint16_t pv_col{0}, pv_row{0};
  float         R_pix{0.f};

  // Plane-level response
  float         top_score{0.f};

  // Up to 10 DV peaks (sorted by score desc; same cap as the script)
  std::vector<VDPeak> peaks;
};

struct VDPerf {
  // Outer (C++) timings
  float t_write_req_ms{0.f};
  float t_exec_total_ms{0.f};
  float t_read_resp_ms{0.f};
  // Child (Python) timings
  float t_child_total_ms{0.f};
  float t_child_setup_ms{0.f};
  float t_child_infer_ms{0.f};
  float t_child_post_ms{0.f};
  // Optional memory stats if the script writes them
  float child_max_rss_mb{0.f};
  float child_cuda_mem_mb{0.f};
};

struct VertexDisplacement {
  std::string   algorithm;                // e.g. "pvdv2d_me"
  float         response{0.f};            // global max over planes

  VDPlaneSummary U, V, W;                 // per-plane summaries
  VDPerf         perf;                    // timings

  // For provenance (fixed in this minimal setup)
  float          theta_xz_deg{8.0f};      // ≈ NuMI off-axis horizontally
  float          theta_yz_deg{-3.0f};     // ≈ NuMI downward tilt
};

} // namespace image

#endif
