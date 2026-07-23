#ifndef IMAGEPIPELINE_INFERENCEPROTOCOL_H
#define IMAGEPIPELINE_INFERENCEPROTOCOL_H

#include "ImagePipeline/InferencePlaneTriplet.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct InferenceChildMetrics {
    double total_ms{0.0};
    double setup_ms{0.0};
    double infer_ms{0.0};
    double post_ms{0.0};
    double max_rss_mb{0.0};
};

/// Encoder/decoder for the IASP v2 request and IAOK v1 result protocols.
class InferenceProtocol {
  public:
    static constexpr uint32_t requestVersion() noexcept { return 2U; }
    static constexpr uint32_t resultVersion() noexcept { return 1U; }
    static constexpr std::size_t defaultMaxLogits() noexcept {
        return 1024U * 1024U;
    }

    static void writeRequest(const std::string &path,
                             const ValidatedPlaneTriplet &planes);

    static std::vector<float>
    readResult(const std::string &path,
               std::size_t max_logits = defaultMaxLogits());

    /// Parse the complete metrics file emitted by infer_bin.py.
    ///
    /// All required performance keys must occur exactly once and every line
    /// must be a recognised, strictly parsed key/value pair.
    static InferenceChildMetrics readMetrics(const std::string &path);
};

} // namespace image

#endif // IMAGEPIPELINE_INFERENCEPROTOCOL_H
