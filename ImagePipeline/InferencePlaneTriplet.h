#ifndef IMAGEPIPELINE_INFERENCEPLANE_TRIPLET_H
#define IMAGEPIPELINE_INFERENCEPLANE_TRIPLET_H

#include "Products/ImageFeatures.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace image {

/// A validated, non-owning U/V/W view over an ImageFeatures collection.
///
/// Construction is the single inference-side boundary for persisted image
/// products. The underlying collection must outlive this object.
class ValidatedPlaneTriplet {
  public:
    explicit ValidatedPlaneTriplet(
        const std::vector<ImageFeatures> &detector_images);

    const ImageFeatures &u() const noexcept { return *planes_[0]; }
    const ImageFeatures &v() const noexcept { return *planes_[1]; }
    const ImageFeatures &w() const noexcept { return *planes_[2]; }
    const ImageFeatures &operator[](std::size_t index) const;

    const std::array<const ImageFeatures *, 3> &planes() const noexcept {
        return planes_;
    }

    uint32_t width() const noexcept { return u().width; }
    uint32_t height() const noexcept { return u().height; }
    uint32_t featureDimension() const noexcept { return u().feature_dim; }
    bool empty() const noexcept;

  private:
    std::array<const ImageFeatures *, 3> planes_{{nullptr, nullptr, nullptr}};
};

} // namespace image

#endif // IMAGEPIPELINE_INFERENCEPLANE_TRIPLET_H
