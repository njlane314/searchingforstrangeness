#ifndef ANALYSIS_INFERENCE_PERF_H
#define ANALYSIS_INFERENCE_PERF_H

#include <string>
#include <vector>

namespace image {

struct ModelPerf {
    std::string model; ///< Model identifier or name.

    float t_write_req_ms{0.f}; ///< Time to write the inference request (ms).
    float t_exec_total_ms{0.f}; ///< Total inference execution time (ms).
    float t_read_resp_ms{0.f}; ///< Time to read the inference response (ms).

    float t_child_total_ms{0.f}; ///< Total child process time (ms).
    float t_child_setup_ms{0.f}; ///< Child process setup time (ms).
    float t_child_infer_ms{0.f}; ///< Child process inference time (ms).
    float t_child_post_ms{0.f}; ///< Child process post-processing time (ms).

    float child_max_rss_mb{0.f}; ///< Peak child memory usage (MB RSS).
};

struct InferencePerfProduct {
    std::vector<ModelPerf> per_model; ///< Performance metrics per model.
};

} // namespace image

#endif // ANALYSIS_INFERENCE_PERF_H
