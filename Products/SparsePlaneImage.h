#ifndef ANALYSIS_SPARSE_PLANE_IMAGE_H
#define ANALYSIS_SPARSE_PLANE_IMAGE_H

#include <cstdint>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

class SparsePlaneImage {
    public:
        int view{static_cast<int>(geo::kUnknown)}; ///< Readout plane/view identifier.
        uint32_t width{0}; ///< Image width in pixels.
        uint32_t height{0}; ///< Image height in pixels.
        uint32_t hit_count{0}; ///< Number of hits contributing to this plane image.
        float origin_x{0.f}; ///< World x-coordinate of the image origin.
        float origin_y{0.f}; ///< World y-coordinate of the image origin.
        float pixel_w{0.f}; ///< Pixel width in world units.
        float pixel_h{0.f}; ///< Pixel height in world units.
        int32_t vertex_row{-1}; ///< Row index of the reco neutrino vertex in this plane (-1 if unavailable).
        int32_t vertex_col{-1}; ///< Column index of the reco neutrino vertex in this plane (-1 if unavailable).
        uint32_t feature_dim{3}; ///< Number of features stored per active pixel (adc, nu-slice mask, dead mask).
        std::vector<int32_t> coords; ///< Interleaved (row, col) coordinates for active pixels, including dead-only pixels if encoded.
        std::vector<float> features; ///< Interleaved feature vectors aligned with coords.
        std::vector<uint8_t> semantic; ///< Optional semantic labels aligned with the active-pixel ordering.
};

} // namespace image

#endif // ANALYSIS_SPARSE_PLANE_IMAGE_H
