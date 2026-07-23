#ifndef ANALYSIS_IMAGE_FEATURES_H
#define ANALYSIS_IMAGE_FEATURES_H

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

inline constexpr uint32_t kImageFeatureDimension = 2U;
inline constexpr uint32_t kImageFeatureDimensionWithDeadMask = 3U;

/// One plane in a three-element U/V/W image collection.
///
/// A window with no active pixels is represented by three metadata-bearing
/// planes with empty coords/features/semantic arrays, never by an absent or
/// zero-length collection.
class ImageFeatures {
  public:
    int view{static_cast<int>(geo::kUnknown)}; ///< Readout plane/view identifier.
    uint32_t width{0}; ///< Image width in pixels.
    uint32_t height{0}; ///< Image height in pixels.
    uint32_t hit_count{0}; ///< Unique contributing hits with retained ADC.
    float origin_x{0.f}; ///< World x-coordinate of the image origin.
    float origin_y{0.f}; ///< World y-coordinate of the image origin.
    float pixel_w{0.f}; ///< Pixel width in world units.
    float pixel_h{0.f}; ///< Pixel height in world units.
    int32_t vertex_row{-1}; ///< Reco-vertex row, or -1 when unavailable.
    int32_t vertex_col{-1}; ///< Reco-vertex column, or -1 when unavailable.
    uint32_t feature_dim{kImageFeatureDimension}; ///< Features stored per pixel.
    std::vector<int32_t> coords; ///< Interleaved (row, column) pairs.
    std::vector<float> features; ///< Feature vectors aligned with coords.
    std::vector<uint8_t> semantic; ///< Optional labels aligned with coords.
};

struct ImageValidationResult {
    std::string error;

    explicit operator bool() const noexcept { return error.empty(); }
};

/// Validate the persisted sparse-plane contract without depending on art.
ImageValidationResult
validateImageFeatures(const ImageFeatures &plane);

struct ValidatedImagePlaneTriplet {
    std::array<const ImageFeatures *, 3> planes{{nullptr, nullptr, nullptr}};
    std::string error;

    explicit operator bool() const noexcept { return error.empty(); }
};

/// Validate a collection and return planes in stable U/V/W order.
ValidatedImagePlaneTriplet
validateImagePlaneTriplet(
    const std::vector<ImageFeatures> &images);

} // namespace image

#endif // ANALYSIS_IMAGE_FEATURES_H
