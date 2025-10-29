#ifndef IMAGE_RANDOMFEATURES_H
#define IMAGE_RANDOMFEATURES_H

#include <cstdint>
#include <string>
#include <vector>

namespace image {

struct RandomFeatures {
  std::string model;       // model key, e.g. "rand_sparse_me"
  std::uint32_t dim{0};    // length of 'values'
  std::uint32_t seed{0};   // RNG seed used by the wrapper (optional meta)
  std::vector<float> values;  // feature vector
};

using RandomFeaturesCollection = std::vector<RandomFeatures>;

} // namespace image

#endif // IMAGE_RANDOMFEATURES_H
