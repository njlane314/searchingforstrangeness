#ifndef ANALYSIS_INFERENCE_PERF_COMPAT_H
#define ANALYSIS_INFERENCE_PERF_COMPAT_H

#include "Products/InferenceMetrics.h"

namespace image {

// Backward-compatible aliases for stale dictionary dependencies and includes.
using ModelPerf = ModelMetrics;
using InferencePerfProduct = InferenceMetrics;

} // namespace image

#endif // ANALYSIS_INFERENCE_PERF_COMPAT_H
