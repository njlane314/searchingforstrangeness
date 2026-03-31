#ifndef ANALYSIS_INFERENCE_METRICS_H
#define ANALYSIS_INFERENCE_METRICS_H

#include <string>
#include <vector>

namespace image {

struct ModelMetrics {
    std::string model; ///< Model identifier or name.

    float t_exec_total_ms{0.f}; ///< Total inference execution time (ms).
    float t_child_infer_ms{0.f}; ///< Child process inference time (ms).
    float child_max_rss_mb{0.f}; ///< Peak child memory usage (MB RSS).
};

struct InferenceMetrics {
    std::vector<ModelMetrics> per_model; ///< Performance metrics per model.
};

} // namespace image

#endif // ANALYSIS_INFERENCE_METRICS_H
