#ifndef ANALYSIS_SEGMENTATION_PRODUCTS_H
#define ANALYSIS_SEGMENTATION_PRODUCTS_H

#include <cstdint>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace analysis {

struct PlaneSegmentationProduct {
  int view{static_cast<int>(geo::kUnknown)};
  uint32_t width{0};
  uint32_t height{0};
  uint32_t nclasses{0};
  std::vector<uint8_t> labels;
  std::vector<float> confidence;
};

} // namespace analysis

#endif // ANALYSIS_SEGMENTATION_PRODUCTS_H
