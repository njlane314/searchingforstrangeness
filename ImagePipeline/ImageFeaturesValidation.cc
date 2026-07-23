#include "Products/ImageFeatures.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

namespace image {

ImageValidationResult
validateImageFeatures(const ImageFeatures &plane) {
    if (plane.view != static_cast<int>(geo::kU) &&
        plane.view != static_cast<int>(geo::kV) &&
        plane.view != static_cast<int>(geo::kW)) {
        return {"unsupported image view " +
                std::to_string(plane.view)};
    }
    if (plane.width == 0U || plane.height == 0U)
        return {"image dimensions must both be non-zero"};

    if (!std::isfinite(plane.origin_x) ||
        !std::isfinite(plane.origin_y) ||
        !std::isfinite(plane.pixel_w) ||
        !std::isfinite(plane.pixel_h) ||
        !(plane.pixel_w > 0.f) ||
        !(plane.pixel_h > 0.f)) {
        return {"image geometry contains non-finite or "
                "non-positive values"};
    }
    if (plane.feature_dim != kImageFeatureDimension &&
        plane.feature_dim != kImageFeatureDimensionWithDeadMask) {
        return {"image feature_dim must be " +
                std::to_string(kImageFeatureDimension) +
                " for [adc, nu-slice] or " +
                std::to_string(kImageFeatureDimensionWithDeadMask) +
                " for [adc, nu-slice, dead], got " +
                std::to_string(plane.feature_dim)};
    }
    if ((plane.coords.size() % 2U) != 0U)
        return {"coords must contain complete (row, col) pairs"};

    bool const vertex_missing =
        plane.vertex_row == -1 && plane.vertex_col == -1;
    bool const vertex_in_bounds =
        plane.vertex_row >= 0 && plane.vertex_col >= 0 &&
        static_cast<uint32_t>(plane.vertex_row) < plane.height &&
        static_cast<uint32_t>(plane.vertex_col) < plane.width;
    if (!vertex_missing && !vertex_in_bounds) {
        return {"vertex pixel must be (-1, -1) or lie inside the "
                "declared image dimensions"};
    }

    const std::size_t active_pixels =
        plane.coords.size() / 2U;
    const uint64_t image_capacity =
        static_cast<uint64_t>(plane.width) *
        static_cast<uint64_t>(plane.height);
    if (static_cast<uint64_t>(active_pixels) > image_capacity) {
        return {"active sparse-pixel count exceeds image capacity"};
    }
    if (active_pixels >
            std::numeric_limits<std::size_t>::max() /
                plane.feature_dim ||
        plane.features.size() !=
            active_pixels * plane.feature_dim) {
        return {"features are not aligned with coords and "
                "feature_dim"};
    }
    if (!plane.semantic.empty() &&
        plane.semantic.size() != active_pixels) {
        return {"semantic labels are not aligned with active pixels"};
    }

    for (std::size_t i = 0U; i < active_pixels; ++i) {
        int32_t const row = plane.coords[2U * i];
        int32_t const column = plane.coords[2U * i + 1U];
        if (row < 0 || column < 0 ||
            static_cast<uint32_t>(row) >= plane.height ||
            static_cast<uint32_t>(column) >= plane.width) {
            return {"sparse coordinate lies outside the declared "
                    "image dimensions"};
        }

        for (std::size_t feature = 0U;
             feature < plane.feature_dim; ++feature) {
            if (!std::isfinite(
                    plane.features[
                        i * plane.feature_dim + feature])) {
                return {"sparse feature contains a non-finite value"};
            }
        }
        float const adc =
            plane.features[i * plane.feature_dim];
        float const slice =
            plane.features[i * plane.feature_dim + 1U];
        if (!(adc > 0.f))
            return {"active sparse pixels must have positive ADC"};
        if (slice != 0.f && slice != 1.f) {
            return {"nu-slice feature must be exactly zero or one"};
        }
    }
    return {};
}

ValidatedImagePlaneTriplet
validateImagePlaneTriplet(
    const std::vector<ImageFeatures> &images) {
    ValidatedImagePlaneTriplet result;

    for (auto const &plane : images) {
        auto const validation = validateImageFeatures(plane);
        if (!validation) {
            result.error = validation.error;
            return result;
        }

        std::size_t index = 0U;
        if (plane.view == static_cast<int>(geo::kV))
            index = 1U;
        else if (plane.view == static_cast<int>(geo::kW))
            index = 2U;

        if (result.planes[index] != nullptr) {
            result.error = "duplicate image view " +
                           std::to_string(plane.view);
            return result;
        }
        result.planes[index] = &plane;
    }

    if (images.size() != result.planes.size() ||
        result.planes[0] == nullptr ||
        result.planes[1] == nullptr ||
        result.planes[2] == nullptr) {
        result.error =
            "image collection must contain exactly one U, V, and W plane";
    }
    return result;
}

} // namespace image
