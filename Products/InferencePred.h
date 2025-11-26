#ifndef ANALYSIS_INFERENCE_PRED_H
#define ANALYSIS_INFERENCE_PRED_H

#include <string>
#include <vector>

namespace image {

struct ModelPred {
  std::string model;
  std::vector<float> scores;
};

struct InferencePredProduct {
  std::vector<ModelPred> per_model;
};

} // namespace image

#endif // ANALYSIS_INFERENCE_PRED_H
