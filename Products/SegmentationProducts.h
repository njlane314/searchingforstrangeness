#ifndef DATAPRODUCTS_SEGMENTATION_PRODUCTS_H
#define DATAPRODUCTS_SEGMENTATION_PRODUCTS_H

#include <cstdint>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

struct PlaneSegmentation {
  int view{static_cast<int>(geo::kUnknown)};
  uint32_t width{0};
  uint32_t height{0};
  uint32_t nclasses{0};
  std::vector<uint8_t> labels;
  std::vector<float> confidence;
};

} // namespace image

#endif // DATAPRODUCTS_SEGMENTATION_PRODUCTS_H
