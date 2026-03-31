#ifndef ANALYSIS_INFERENCE_PRED_COMPAT_H
#define ANALYSIS_INFERENCE_PRED_COMPAT_H

#include "Products/InferencePredictions.h"

namespace image {

// Backward-compatible aliases for stale dictionary dependencies and includes.
using ModelPred = ModelPrediction;
using InferencePredProduct = InferencePredictions;

} // namespace image

#endif // ANALYSIS_INFERENCE_PRED_COMPAT_H
