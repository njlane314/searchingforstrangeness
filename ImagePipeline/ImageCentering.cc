#include "ImagePipeline/ImageCentering.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

namespace image {

std::optional<TVector3>
weightedCentroid(
    const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
    const std::vector<double> &weights) {
    TVector3 sum(0., 0., 0.);
    double total_weight = 0.0;

    const std::size_t count =
        std::min(spacepoints.size(), weights.size());
    for (std::size_t i = 0U; i < count; ++i) {
        auto const &spacepoint = spacepoints[i];
        double const weight = weights[i];
        if (!spacepoint || !std::isfinite(weight) ||
            !(weight > 0.0)) {
            continue;
        }

        auto const *coordinates = spacepoint->XYZ();
        TVector3 const point(
            coordinates[0], coordinates[1], coordinates[2]);
        if (!std::isfinite(point.X()) ||
            !std::isfinite(point.Y()) ||
            !std::isfinite(point.Z())) {
            continue;
        }

        sum += weight * point;
        total_weight += weight;
    }

    if (!std::isfinite(total_weight) ||
        !(total_weight > 0.0)) {
        return std::nullopt;
    }

    TVector3 const center = sum * (1.0 / total_weight);
    if (!std::isfinite(center.X()) ||
        !std::isfinite(center.Y()) ||
        !std::isfinite(center.Z())) {
        return std::nullopt;
    }
    return center;
}

TVector3 trimmedCentroid(
    const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
    const std::vector<double> &weights,
    const TVector3 &seed,
    double radius) {
    if (!std::isfinite(radius) || !(radius > 0.0) ||
        !std::isfinite(seed.X()) ||
        !std::isfinite(seed.Y()) ||
        !std::isfinite(seed.Z())) {
        return seed;
    }

    double const radius_squared = radius * radius;
    TVector3 sum(0., 0., 0.);
    double total_weight = 0.0;

    const std::size_t count =
        std::min(spacepoints.size(), weights.size());
    for (std::size_t i = 0U; i < count; ++i) {
        auto const &spacepoint = spacepoints[i];
        double const weight = weights[i];
        if (!spacepoint || !std::isfinite(weight) ||
            !(weight > 0.0)) {
            continue;
        }

        auto const *coordinates = spacepoint->XYZ();
        TVector3 const point(
            coordinates[0], coordinates[1], coordinates[2]);
        if (!std::isfinite(point.X()) ||
            !std::isfinite(point.Y()) ||
            !std::isfinite(point.Z()) ||
            (point - seed).Mag2() > radius_squared) {
            continue;
        }

        sum += weight * point;
        total_weight += weight;
    }

    if (!(total_weight > 0.0))
        return seed;
    return sum * (1.0 / total_weight);
}

ImageCentering::ImageCentering(
    art::InputTag hit_producer,
    art::InputTag spacepoint_producer)
    : hit_producer_{std::move(hit_producer)}
    , spacepoint_producer_{std::move(spacepoint_producer)}
{}

ImageCenter ImageCentering::compute(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
    const std::optional<TVector3> &vertex,
    double trimming_radius) const {
    if (neutrino_hits.empty()) {
        if (vertex) {
            return ImageCenter{
                *vertex, ImageCenterSeed::Vertex};
        }
        return ImageCenter{
            TVector3(0., 0., 0.),
            ImageCenterSeed::OriginFallback};
    }

    std::map<art::Ptr<recob::SpacePoint>, double>
        spacepoint_charge;
    auto const hit_handle =
        event.getValidHandle<std::vector<recob::Hit>>(
            hit_producer_);
    art::FindManyP<recob::SpacePoint> hit_to_spacepoints(
        hit_handle, event, spacepoint_producer_);

    for (auto const &hit : neutrino_hits) {
        if (!hit)
            continue;

        double const charge = hit->Integral();
        if (!std::isfinite(charge) || !(charge > 0.0))
            continue;

        auto const &spacepoints =
            hit_to_spacepoints.at(hit.key());
        if (spacepoints.empty())
            continue;

        double const charge_per_spacepoint =
            charge / static_cast<double>(spacepoints.size());
        for (auto const &spacepoint : spacepoints) {
            if (spacepoint) {
                spacepoint_charge[spacepoint] +=
                    charge_per_spacepoint;
            }
        }
    }

    std::vector<art::Ptr<recob::SpacePoint>> spacepoints;
    std::vector<double> weights;
    spacepoints.reserve(spacepoint_charge.size());
    weights.reserve(spacepoint_charge.size());
    for (auto const &entry : spacepoint_charge) {
        spacepoints.push_back(entry.first);
        weights.push_back(entry.second);
    }

    ImageCenter result;
    TVector3 reference(0., 0., 0.);
    if (vertex) {
        reference = *vertex;
        result.seed = ImageCenterSeed::Vertex;
    } else if (auto const centroid =
                   weightedCentroid(spacepoints, weights)) {
        reference = *centroid;
        result.seed = ImageCenterSeed::WeightedCentroid;
    }
    result.position = reference;

    if (result.seed != ImageCenterSeed::OriginFallback &&
        std::isfinite(trimming_radius) &&
        trimming_radius > 0.0 && !spacepoints.empty()) {
        result.position = trimmedCentroid(
            spacepoints, weights, reference, trimming_radius);
    }

    if (!std::isfinite(result.position.X()) ||
        !std::isfinite(result.position.Y()) ||
        !std::isfinite(result.position.Z())) {
        result.position.SetXYZ(0., 0., 0.);
        result.seed = ImageCenterSeed::OriginFallback;
    }
    return result;
}

} // namespace image
