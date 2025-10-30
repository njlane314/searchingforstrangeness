#ifndef ANALYSIS_INFERENCE_PERF_H
#define ANALYSIS_INFERENCE_PERF_H

#include <string>
#include <vector>

namespace image {

struct ModelPerf {
  std::string model;

  float t_write_req_ms{0.f};
  float t_exec_total_ms{0.f};
  float t_read_resp_ms{0.f};

  float t_child_total_ms{0.f};
  float t_child_setup_ms{0.f};
  float t_child_infer_ms{0.f};
  float t_child_post_ms{0.f};

  float child_max_rss_mb{0.f};
};

struct InferencePerfProduct {
  std::vector<ModelPerf> per_model;
};

} // namespace image

#endif // ANALYSIS_INFERENCE_PERF_H
