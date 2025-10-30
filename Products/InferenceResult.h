#ifndef ANALYSIS_INFERENCE_RESULT_H
#define ANALYSIS_INFERENCE_RESULT_H

#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct ModelResult {
  std::string model;

  std::vector<float> cls_scores;

  std::uint32_t segmentation_width{0};
  std::uint32_t segmentation_height{0};

  std::vector<std::uint8_t> segmentation_u;
  std::vector<std::uint8_t> segmentation_v;
  std::vector<std::uint8_t> segmentation_w;

  std::vector<float> confidence_u;
  std::vector<float> confidence_v;
  std::vector<float> confidence_w;

  std::vector<float> features;
  std::uint32_t feature_seed{0};
};

struct InferenceResultProduct {
  std::vector<ModelResult> per_model;
};

} // namespace image

#endif // ANALYSIS_INFERENCE_RESULT_H
