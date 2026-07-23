#include "ImagePipeline/ImageWindowGeometry.h"

#include <cetlib_except/exception.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace {

image::ImageProperties makePropertiesFromBounds(
    double min_wire_coordinate,
    double max_wire_coordinate,
    double min_drift,
    double max_drift,
    double pixel_width,
    double pixel_height,
    geo::View_t view) {
    if (!std::isfinite(min_wire_coordinate) ||
        !std::isfinite(max_wire_coordinate) ||
        !std::isfinite(min_drift) ||
        !std::isfinite(max_drift)) {
        throw cet::exception("ImageWindowGeometry")
            << "Cannot build image properties from non-finite bounds.";
    }
    if (!(pixel_width > 0.0) ||
        !(pixel_height > 0.0)) {
        throw cet::exception("ImageWindowGeometry")
            << "Cannot build image properties with non-positive "
               "pixel size.";
    }
    if (max_wire_coordinate < min_wire_coordinate ||
        max_drift < min_drift) {
        throw cet::exception("ImageWindowGeometry")
            << "Cannot build image properties from inverted bounds.";
    }

    auto const width = std::max<std::size_t>(
        1U,
        static_cast<std::size_t>(std::ceil(
            (max_wire_coordinate - min_wire_coordinate) /
            pixel_width)) +
            1U);
    auto const height = std::max<std::size_t>(
        1U,
        static_cast<std::size_t>(std::ceil(
            (max_drift - min_drift) / pixel_height)) +
            1U);
    double const center_wire =
        min_wire_coordinate +
        0.5 * static_cast<double>(width) * pixel_width;
    double const center_drift =
        min_drift +
        0.5 * static_cast<double>(height) * pixel_height;

    return image::ImageProperties(
        center_wire, center_drift, width, height,
        pixel_height, pixel_width, view);
}

} // namespace

namespace image {

geo::View_t canonicalImageView(geo::View_t view) {
    return view == geo::kY ? geo::kW : view;
}

std::size_t imageViewIndex(geo::View_t view) {
    switch (canonicalImageView(view)) {
    case geo::kU:
        return 0U;
    case geo::kV:
        return 1U;
    case geo::kW:
        return 2U;
    default:
        throw cet::exception("ImageWindowGeometry")
            << "Unsupported geo::View_t: "
            << static_cast<int>(view);
    }
}

TVector3 projectToImageView(
    float x, float y, float z, geo::View_t view) {
    float constexpr u_angle = 1.04719758034f;
    float constexpr v_angle = -1.04719758034f;

    float wire_coordinate = z;
    switch (canonicalImageView(view)) {
    case geo::kU:
        wire_coordinate =
            z * std::cos(u_angle) - y * std::sin(u_angle);
        break;
    case geo::kV:
        wire_coordinate =
            z * std::cos(v_angle) - y * std::sin(v_angle);
        break;
    case geo::kW:
        break;
    default:
        throw cet::exception("ImageWindowGeometry")
            << "Unsupported geo::View_t: "
            << static_cast<int>(view);
    }
    return TVector3(x, 0.f, wire_coordinate);
}

ImageWindowGeometry::ImageWindowGeometry(
    geo::GeometryCore const &geometry,
    ImageWindowGeometryConfig config)
    : geometry_{&geometry}, config_{std::move(config)}
{
    if (config_.width == 0U || config_.height == 0U ||
        !std::isfinite(config_.drift_step) ||
        !std::isfinite(config_.pitch_u) ||
        !std::isfinite(config_.pitch_v) ||
        !std::isfinite(config_.pitch_w) ||
        !(config_.drift_step > 0.0) ||
        !(config_.pitch_u > 0.0) ||
        !(config_.pitch_v > 0.0) ||
        !(config_.pitch_w > 0.0)) {
        throw cet::exception("ImageWindowGeometry")
            << "Image dimensions and pixel sizes must be finite "
               "and positive.";
    }
    full_window_properties_ = computeFullWindowProperties();
}

std::vector<ImageProperties>
ImageWindowGeometry::croppedWindowProperties(
    const TVector3 &center) const {
    if (!std::isfinite(center.X()) ||
        !std::isfinite(center.Y()) ||
        !std::isfinite(center.Z())) {
        throw cet::exception("ImageWindowGeometry")
            << "Cannot center a cropped image on a non-finite point.";
    }

    auto const projected_u = projectToImageView(
        static_cast<float>(center.X()),
        static_cast<float>(center.Y()),
        static_cast<float>(center.Z()), geo::kU);
    auto const projected_v = projectToImageView(
        static_cast<float>(center.X()),
        static_cast<float>(center.Y()),
        static_cast<float>(center.Z()), geo::kV);
    auto const projected_w = projectToImageView(
        static_cast<float>(center.X()),
        static_cast<float>(center.Y()),
        static_cast<float>(center.Z()), geo::kW);

    std::vector<ImageProperties> properties;
    properties.reserve(3U);
    properties.emplace_back(
        projected_u.Z(), projected_u.X(),
        config_.width, config_.height, config_.drift_step,
        config_.pitch_u, geo::kU);
    properties.emplace_back(
        projected_v.Z(), projected_v.X(),
        config_.width, config_.height, config_.drift_step,
        config_.pitch_v, geo::kV);
    properties.emplace_back(
        projected_w.Z(), projected_w.X(),
        config_.width, config_.height, config_.drift_step,
        config_.pitch_w, geo::kW);
    return properties;
}

const std::vector<ImageProperties> &
ImageWindowGeometry::fullWindowProperties() const noexcept {
    return full_window_properties_;
}

double ImageWindowGeometry::trimmingRadius() const noexcept {
    return 0.5 * std::min(
        static_cast<double>(config_.height) *
            config_.drift_step,
        static_cast<double>(config_.width) *
            config_.pitch_w);
}

std::vector<ImageProperties>
ImageWindowGeometry::computeFullWindowProperties() const {
    auto const &tpc = geometry_->TPC();
    auto const bounds = tpc.ActiveBoundingBox();
    double const min_drift = bounds.MinX();
    double const max_drift = bounds.MaxX();

    auto make_view = [&](unsigned int plane_index,
                         geo::View_t view, double pitch) {
        auto const &plane = tpc.Plane(plane_index);
        double min_wire =
            std::numeric_limits<double>::max();
        double max_wire =
            std::numeric_limits<double>::lowest();

        for (std::size_t wire_index = 0U;
             wire_index < plane.Nwires(); ++wire_index) {
            TVector3 const center =
                plane.Wire(
                    static_cast<unsigned int>(wire_index))
                    .GetCenter();
            TVector3 const projected = projectToImageView(
                static_cast<float>(center.X()),
                static_cast<float>(center.Y()),
                static_cast<float>(center.Z()), view);
            min_wire = std::min(
                min_wire,
                static_cast<double>(projected.Z()));
            max_wire = std::max(
                max_wire,
                static_cast<double>(projected.Z()));
        }

        return makePropertiesFromBounds(
            min_wire, max_wire, min_drift, max_drift,
            pitch, config_.drift_step, view);
    };

    std::vector<ImageProperties> properties;
    properties.reserve(3U);
    properties.emplace_back(
        make_view(0U, geo::kU, config_.pitch_u));
    properties.emplace_back(
        make_view(1U, geo::kV, config_.pitch_v));
    properties.emplace_back(
        make_view(2U, geo::kW, config_.pitch_w));
    return properties;
}

} // namespace image
