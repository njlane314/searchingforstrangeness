#ifndef ANALYSIS_INFERENCE_PREDICTIONS_H
#define ANALYSIS_INFERENCE_PREDICTIONS_H

#include <string>
#include <vector>

namespace image {

struct ModelPrediction {
    std::string model; ///< Model identifier or name.
    std::vector<float> scores; ///< Raw model outputs for the model (currently logits).
};

struct InferencePredictions {
    std::vector<ModelPrediction> per_model; ///< Predictions grouped by model.
};

} // namespace image

#endif // ANALYSIS_INFERENCE_PREDICTIONS_H
