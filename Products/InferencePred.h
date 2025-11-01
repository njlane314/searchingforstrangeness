#ifndef ANALYSIS_INFERENCE_PRED_H
#define ANALYSIS_INFERENCE_PRED_H

#include <string>
#include <vector>

namespace image {

struct ModelPred {
  std::string model;          // model label or arch
  std::vector<float> scores;  // K float32 scores read from results.bin
};

struct InferencePredProduct {
  std::vector<ModelPred> per_model;
};

} // namespace image

#endif // ANALYSIS_INFERENCE_PRED_H
