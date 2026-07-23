#ifndef IMAGEPIPELINE_INFERENCECONFIGURATION_H
#define IMAGEPIPELINE_INFERENCECONFIGURATION_H

#include "ImagePipeline/InferenceProduction.h"

#include "fhiclcpp/ParameterSet.h"

#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct InferenceModelConfig {
    std::string label;
    std::string arch;
    std::string weights;
    std::string weights_dir;
    std::string checkpoint;
    std::string wrapper;
    std::string runtime;
    std::string container;
    double timeout_seconds{0.0};
    uint32_t diagnostic_bytes{64U * 1024U};
    uint32_t max_logits{1024U * 1024U};

    std::string productName() const;
};

/// Parses and validates the inference FHiCL configuration once per module.
class InferenceConfiguration {
  public:
    explicit InferenceConfiguration(
        const fhicl::ParameterSet &parameters);

    const std::vector<InferenceModelConfig> &models() const noexcept {
        return models_;
    }
    const std::string &scratchDirectory() const noexcept {
        return scratch_dir_;
    }

    InferenceProduction::Config
    makeRunConfig(const InferenceModelConfig &model,
                  const std::string &event_context) const;

  private:
    static InferenceModelConfig
    mergeModelConfig(const fhicl::ParameterSet &parameters,
                     InferenceModelConfig defaults);

    std::string runtime_base_dir_;
    std::string default_wrapper_;
    std::string scratch_dir_;
    std::vector<InferenceModelConfig> models_;
};

} // namespace image

#endif // IMAGEPIPELINE_INFERENCECONFIGURATION_H
