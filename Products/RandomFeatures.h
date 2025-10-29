#ifndef ANALYSIS_RANDOM_FEATURES_H
#define ANALYSIS_RANDOM_FEATURES_H

#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct RandomFeatures {
  std::string model;
  std::vector<float> values;
  std::uint32_t seed{0};
};

} // namespace image

#endif // ANALYSIS_RANDOM_FEATURES_H
