#include "ImagePipeline/InferencePlaneTriplet.h"

#include <sstream>
#include <stdexcept>

namespace image {

ValidatedPlaneTriplet::ValidatedPlaneTriplet(
    const std::vector<ImageFeatures> &detector_images) {
    auto validated = validateImagePlaneTriplet(detector_images);
    if (!validated) {
        throw std::invalid_argument(
            "Invalid inference image collection: " + validated.error);
    }

    planes_ = validated.planes;
    if (u().width != v().width || u().width != w().width ||
        u().height != v().height || u().height != w().height) {
        std::ostringstream message;
        message << "Inference U/V/W dimensions must match: "
                << "U=" << u().width << "x" << u().height << ", "
                << "V=" << v().width << "x" << v().height << ", "
                << "W=" << w().width << "x" << w().height;
        throw std::invalid_argument(message.str());
    }

    if (u().feature_dim != v().feature_dim ||
        u().feature_dim != w().feature_dim) {
        std::ostringstream message;
        message << "Inference U/V/W feature dimensions must match: "
                << "U=" << u().feature_dim << ", "
                << "V=" << v().feature_dim << ", "
                << "W=" << w().feature_dim;
        throw std::invalid_argument(message.str());
    }
}

const ImageFeatures &
ValidatedPlaneTriplet::operator[](const std::size_t index) const {
    if (index >= planes_.size()) {
        throw std::out_of_range("Inference plane index must be in [0, 3)");
    }
    return *planes_[index];
}

bool ValidatedPlaneTriplet::empty() const noexcept {
    return u().coords.empty() && v().coords.empty() && w().coords.empty();
}

} // namespace image
