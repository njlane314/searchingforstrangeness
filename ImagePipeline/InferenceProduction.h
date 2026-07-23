#ifndef IMAGEPIPELINE_INFERENCEPRODUCTION_H
#define IMAGEPIPELINE_INFERENCEPRODUCTION_H

#include "ImagePipeline/InferencePlaneTriplet.h"

#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace image {

class InferenceProduction {
  public:
    struct Perf {
        double t_write_req_ms{0.0};
        double t_exec_total_ms{0.0};
        double t_read_resp_ms{0.0};
        double t_child_total_ms{0.0};
        double t_child_setup_ms{0.0};
        double t_child_infer_ms{0.0};
        double t_child_post_ms{0.0};
        double child_max_rss_mb{0.0};
    };

    struct Result {
        std::vector<float> cls;
        Perf perf;
    };

    struct Config {
        std::string scratch_dir;
        std::string context;
        std::string arch;
        std::string weights_file;
        std::string inference_wrapper;
        std::string runtime_base_dir;
        std::string container_image;
        std::chrono::milliseconds timeout{0};
        std::size_t diagnostic_bytes{64U * 1024U};
        std::size_t max_logits{1024U * 1024U};
    };

    static Result runInference(const ValidatedPlaneTriplet &detector_images,
                               const Config &config);
};

} // namespace image

#endif // IMAGEPIPELINE_INFERENCEPRODUCTION_H
