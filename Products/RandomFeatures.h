#ifndef DATAPRODUCTS_RANDOM_FEATURES_H
#define DATAPRODUCTS_RANDOM_FEATURES_H

#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct RandomFeatures {
  std::string model;               // Inference model that produced these features
  std::uint32_t seed{0};           // RNG seed used by the wrapper (if available)
  std::vector<float> values;       // Flat feature vector returned by the wrapper
};

} // namespace image

#endif // DATAPRODUCTS_RANDOM_FEATURES_H
