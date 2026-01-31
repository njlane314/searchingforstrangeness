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
        float origin_x{0.f}; ///< World x-coordinate of the image origin.
        float origin_y{0.f}; ///< World y-coordinate of the image origin.
        float pixel_w{0.f}; ///< Pixel width in world units.
        float pixel_h{0.f}; ///< Pixel height in world units.
        std::vector<float> adc; ///< Flattened ADC values in row-major order.
        std::vector<int32_t> semantic; ///< Flattened semantic labels aligned with adc.
};

} // namespace image

#endif // DATAPRODUCTS_IMAGE_PRODUCTS_H
