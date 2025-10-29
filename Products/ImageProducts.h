#ifndef DATAPRODUCTS_IMAGE_PRODUCTS_H
#define DATAPRODUCTS_IMAGE_PRODUCTS_H

#include <cstdint>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

struct PlaneImage {
  int view{static_cast<int>(geo::kUnknown)};
  uint32_t width{0};
  uint32_t height{0};
  float origin_x{0.f};
  float origin_y{0.f};
  float pixel_w{0.f};
  float pixel_h{0.f};
  std::vector<float> adc;
  std::vector<int32_t> semantic;
};

} // namespace image

#endif // DATAPRODUCTS_IMAGE_PRODUCTS_H
