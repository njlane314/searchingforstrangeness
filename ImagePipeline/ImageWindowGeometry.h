#ifndef IMAGEPIPELINE_IMAGEWINDOWGEOMETRY_H
#define IMAGEPIPELINE_IMAGEWINDOWGEOMETRY_H

#include "ImagePipeline/Image.h"

#include "larcorealg/Geometry/GeometryCore.h"

#include <TVector3.h>

#include <cstddef>
#include <vector>

namespace image {

geo::View_t canonicalImageView(geo::View_t view);
std::size_t imageViewIndex(geo::View_t view);
TVector3 projectToImageView(
    float x, float y, float z, geo::View_t view);

struct ImageWindowGeometryConfig {
    std::size_t width{512U};
    std::size_t height{512U};
    double drift_step{0.0};
    double pitch_u{0.0};
    double pitch_v{0.0};
    double pitch_w{0.0};
};

class ImageWindowGeometry {
  public:
    ImageWindowGeometry(
        geo::GeometryCore const &geometry,
        ImageWindowGeometryConfig config);

    std::vector<ImageProperties>
    croppedWindowProperties(const TVector3 &centre) const;
    const std::vector<ImageProperties> &
    fullWindowProperties() const noexcept;
    double trimmingRadius() const noexcept;

  private:
    std::vector<ImageProperties>
    computeFullWindowProperties() const;

    geo::GeometryCore const *geometry_{nullptr};
    ImageWindowGeometryConfig config_;
    std::vector<ImageProperties> full_window_properties_;
};

} // namespace image

#endif
