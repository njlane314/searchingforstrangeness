#ifndef DATAPRODUCTS_IMAGE_PRODUCTS_H
#define DATAPRODUCTS_IMAGE_PRODUCTS_H

#include <cstdint>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

class ImageProduct {
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
        // Sparse payload: index[k] maps to adc[k] and semantic[k] for one active pixel.
        std::vector<uint32_t> index; ///< Flattened pixel index in row-major order for active pixels.
        std::vector<float> adc; ///< ADC values for active pixels.
        std::vector<uint8_t> semantic; ///< Semantic labels aligned with adc/index.
};

} // namespace image

#endif // DATAPRODUCTS_IMAGE_PRODUCTS_H
