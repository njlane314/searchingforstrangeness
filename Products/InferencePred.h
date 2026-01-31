#ifndef ANALYSIS_INFERENCE_PRED_H
#define ANALYSIS_INFERENCE_PRED_H

#include <string>
#include <vector>

namespace image {

struct ModelPred {
    std::string model; ///< Model identifier or name.
    std::vector<float> scores; ///< Prediction scores for the model output.
};

struct InferencePredProduct {
    std::vector<ModelPred> per_model; ///< Predictions grouped by model.
};

} // namespace image

#endif // ANALYSIS_INFERENCE_PRED_H
