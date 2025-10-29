#ifndef DATAPRODUCTS_SEGMENTATION_PRODUCTS_H
#define DATAPRODUCTS_SEGMENTATION_PRODUCTS_H

#include <cstdint>
#include <string>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

struct PlaneSegmentation {
  std::string model;                             // Inference model that produced this segmentation
  int         view{static_cast<int>(geo::kUnknown)};
  uint32_t    width{0};
  uint32_t    height{0};
  uint32_t    nclasses{0};
  std::vector<uint8_t> labels;                  // Flattened label image (row-major)
  std::vector<float>   confidence;              // Optional per-pixel confidence (same layout as labels)
};

} // namespace image

#endif // DATAPRODUCTS_SEGMENTATION_PRODUCTS_H
