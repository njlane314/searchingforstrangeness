#ifndef DATAPRODUCTS_IMAGE_PRODUCTS_H
#define DATAPRODUCTS_IMAGE_PRODUCTS_H

#include <cstdint>
#include <string>
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

struct InferenceScores {
  std::string               model;   // Name of the inference model that generated the scores
  std::vector<std::string>  names;   // Optional human-readable class labels
  std::vector<float>        scores;  // Raw classification scores (aligned with `names` when provided)
};

} // namespace image

#endif // DATAPRODUCTS_IMAGE_PRODUCTS_H
