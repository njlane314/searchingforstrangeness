#include "ImagePipeline/ImageProduction.h"

#include <cetlib_except/exception.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace image {

ImageEventContext::ImageEventContext(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<art::Ptr<recob::Hit>> &slice_hits,
    PixelImageOptions const &options)
    : wires_{event.getValidHandle<std::vector<recob::Wire>>(
          options.producers.wire)}
    , all_hits_{event.getValidHandle<std::vector<recob::Hit>>(
          options.producers.hit)}
    , wire_hit_assoc_{wires_, event, options.producers.hit}
{
    if (!wire_hit_assoc_.isValid()) {
        throw cet::exception("ImageProduction")
            << "Wire-to-hit associations are unavailable for tag '"
            << options.producers.hit.encode() << "'.";
    }

    for (auto const &hit : hits) {
        if (hit)
            accepted_hits_.insert(hit);
    }
    for (auto const &hit : slice_hits) {
        if (hit)
            slice_hit_keys_.insert(
                static_cast<std::size_t>(hit.key()));
    }

    if (options.semantic == nullptr)
        return;

    const bool found_mcps =
        event.getByLabel(options.producers.mcp, mcp_vector_);
    if (!found_mcps || !mcp_vector_.isValid()) {
        throw cet::exception("ImageProduction")
            << "Semantic labeling is enabled, but the MCParticle "
               "collection '"
            << options.producers.mcp.encode()
            << "' is unavailable. Set ImageProducer IsData=true for "
               "data input.";
    }

    mcp_bkth_assoc_ = std::make_unique<BacktrackerAssociations>(
        all_hits_, event, options.producers.bkt);
    if (!mcp_bkth_assoc_->isValid()) {
        throw cet::exception("ImageProduction")
            << "Hit-to-MCParticle associations are unavailable for tag '"
            << options.producers.bkt.encode() << "'.";
    }
    semantic_labels_ =
        options.semantic->classifyParticles(event);

    for (std::size_t i = 0U; i < mcp_vector_->size(); ++i)
        trackid_to_index_[mcp_vector_->at(i).TrackId()] = i;
}

ImageEventContext::~ImageEventContext() = default;

const std::vector<recob::Wire> &
ImageEventContext::wires() const noexcept {
    return *wires_;
}

const art::FindManyP<recob::Hit> &
ImageEventContext::wireHitAssociations() const noexcept {
    return wire_hit_assoc_;
}

bool ImageEventContext::acceptsHit(
    const art::Ptr<recob::Hit> &hit) const {
    return hit && accepted_hits_.count(hit) != 0U;
}

bool ImageEventContext::isSliceHit(
    const art::Ptr<recob::Hit> &hit) const {
    return hit &&
           slice_hit_keys_.count(
               static_cast<std::size_t>(hit.key())) != 0U;
}

bool ImageEventContext::hasSemantics() const noexcept {
    return mcp_vector_.isValid() && mcp_bkth_assoc_ != nullptr;
}

const ImageEventContext::BacktrackerAssociations *
ImageEventContext::backtrackerAssociations() const noexcept {
    return mcp_bkth_assoc_.get();
}

const std::map<int, std::size_t> &
ImageEventContext::trackIdToIndex() const noexcept {
    return trackid_to_index_;
}

const std::vector<image::SemanticClassifier::SemanticLabel> &
ImageEventContext::semanticLabels() const noexcept {
    return semantic_labels_;
}

RasterizedImageWindow::RasterizedImageWindow(
    std::vector<ImageProperties> properties)
    : properties_{std::move(properties)}
{
    detector_images_.reserve(properties_.size());
    semantic_images_.reserve(properties_.size());
    slice_images_.reserve(properties_.size());
    contributing_hit_keys_.resize(properties_.size());

    for (auto const &property : properties_) {
        detector_images_.emplace_back(property);
        semantic_images_.emplace_back(property);
        slice_images_.emplace_back(property);
        semantic_images_.back().clear(
            static_cast<int>(
                image::SemanticClassifier::SemanticLabel::Empty));
    }
}

const std::vector<ImageProperties> &
RasterizedImageWindow::properties() const noexcept {
    return properties_;
}

const std::vector<Image<float>> &
RasterizedImageWindow::detectorImages() const noexcept {
    return detector_images_;
}

const std::vector<Image<int>> &
RasterizedImageWindow::semanticImages() const noexcept {
    return semantic_images_;
}

const std::vector<Image<uint8_t>> &
RasterizedImageWindow::sliceImages() const noexcept {
    return slice_images_;
}

std::size_t RasterizedImageWindow::contributingHitCount(
    std::size_t plane_index) const {
    return contributing_hit_keys_.at(plane_index).size();
}

struct ImageRasterizer::WirePreparation {
    geo::PlaneID plane_id;
    geo::View_t view{geo::kUnknown};
    double wire_coordinate{0.0};
    std::vector<art::Ptr<recob::Hit>> hits;
};

struct ImageRasterizer::Destination {
    RasterizedImageWindow *window{nullptr};
    std::size_t plane_index{0U};
    std::size_t column{0U};
};

ImageRasterizer::ImageRasterizer(
    geo::GeometryCore const &geometry)
    : geometry_{&geometry}
{}

void ImageRasterizer::rasterize(
    const ImageEventContext &event_context,
    const std::vector<RasterizedImageWindow *> &windows,
    detinfo::DetectorProperties const *detector_properties) const {
    if (detector_properties == nullptr) {
        throw cet::exception("ImageProduction")
            << "ImageRasterizer called with null "
               "DetectorProperties pointer.";
    }
    for (auto const *window : windows) {
        if (window == nullptr) {
            throw cet::exception("ImageProduction")
                << "ImageRasterizer called with a null output window.";
        }
    }

    auto const &event_wires = event_context.wires();
    for (std::size_t wire_index = 0U;
         wire_index < event_wires.size(); ++wire_index) {
        fillWire(wire_index, event_context, windows,
                 *detector_properties);
    }
}

std::optional<ImageRasterizer::WirePreparation>
ImageRasterizer::prepareWire(
    std::size_t wire_index,
    const ImageEventContext &context) const {
    auto hits_for_wire =
        context.wireHitAssociations().at(wire_index);

    std::vector<art::Ptr<recob::Hit>> filtered_hits;
    filtered_hits.reserve(hits_for_wire.size());
    for (auto const &hit : hits_for_wire) {
        if (context.acceptsHit(hit))
            filtered_hits.push_back(hit);
    }
    if (filtered_hits.empty())
        return std::nullopt;

    geo::WireID const &wire_id =
        filtered_hits.front()->WireID();
    geo::PlaneID const plane_id = wire_id.planeID();
    geo::View_t const view =
        canonicalImageView(geometry_->View(plane_id));
    (void)imageViewIndex(view);

    geo::WireGeo const &wire_geometry =
        geometry_->WireIDToWireGeo(wire_id);
    TVector3 const wire_center = wire_geometry.GetCenter();
    TVector3 const projected = projectToImageView(
        static_cast<float>(wire_center.X()),
        static_cast<float>(wire_center.Y()),
        static_cast<float>(wire_center.Z()),
        view);

    return WirePreparation{
        plane_id, view, projected.Z(), std::move(filtered_hits)};
}

std::vector<ImageRasterizer::Destination>
ImageRasterizer::destinationsForWire(
    geo::View_t view, double wire_coordinate,
    const std::vector<RasterizedImageWindow *> &windows) {
    std::vector<Destination> destinations;
    destinations.reserve(windows.size());

    for (auto *window : windows) {
        for (std::size_t plane_index = 0U;
             plane_index < window->properties_.size();
             ++plane_index) {
            auto const &property =
                window->properties_[plane_index];
            if (canonicalImageView(property.view()) != view)
                continue;

            auto const column = property.col(wire_coordinate);
            if (column) {
                destinations.push_back(
                    Destination{window, plane_index, *column});
            }
            break;
        }
    }
    return destinations;
}

image::SemanticClassifier::SemanticLabel
ImageRasterizer::semanticLabel(
    const art::Ptr<recob::Hit> &hit,
    const ImageEventContext &context) {
    using SemanticLabel =
        image::SemanticClassifier::SemanticLabel;
    SemanticLabel const fallback = SemanticLabel::Cosmic;
    if (!context.hasSemantics() || !hit)
        return fallback;

    std::vector<art::Ptr<simb::MCParticle>> particles;
    std::vector<anab::BackTrackerHitMatchingData const *>
        matching_data;
    context.backtrackerAssociations()->get(
        static_cast<std::size_t>(hit.key()), particles,
        matching_data);

    float best_fraction = -1.0f;
    int best_track_id = -1;
    const std::size_t associations =
        std::min(particles.size(), matching_data.size());
    for (std::size_t i = 0U; i < associations; ++i) {
        auto const *data = matching_data[i];
        auto const &particle = particles[i];
        if (data == nullptr || !particle ||
            !std::isfinite(data->ideFraction) ||
            !(data->ideFraction > best_fraction))
            continue;
        best_fraction = data->ideFraction;
        best_track_id = particle->TrackId();
    }
    if (best_track_id == -1)
        return fallback;

    auto const track =
        context.trackIdToIndex().find(best_track_id);
    if (track == context.trackIdToIndex().end() ||
        track->second >= context.semanticLabels().size())
        return fallback;

    if (best_fraction <= 0.5f)
        return SemanticLabel::Ambiguous;
    return context.semanticLabels()[track->second];
}

void ImageRasterizer::mergeSemanticPixel(
    Image<int> &image, std::size_t row, std::size_t column,
    image::SemanticClassifier::SemanticLabel label) {
    int const empty = static_cast<int>(
        image::SemanticClassifier::SemanticLabel::Empty);
    int const ambiguous = static_cast<int>(
        image::SemanticClassifier::SemanticLabel::Ambiguous);
    int const next = static_cast<int>(label);
    int const current = image.get(row, column);

    if (current == empty) {
        image.set(row, column, next, false);
    } else if (current != next && current != ambiguous) {
        image.set(row, column, ambiguous, false);
    }
}

void ImageRasterizer::fillWire(
    std::size_t wire_index,
    const ImageEventContext &context,
    const std::vector<RasterizedImageWindow *> &windows,
    detinfo::DetectorProperties const &detector_properties) const {
    auto prepared = prepareWire(wire_index, context);
    if (!prepared)
        return;
    auto const &wire = context.wires()[wire_index];
    auto const &wire_data = *prepared;

    auto destinations = destinationsForWire(
        wire_data.view, wire_data.wire_coordinate, windows);
    if (destinations.empty())
        return;

    int min_tick = std::numeric_limits<int>::max();
    int max_tick = std::numeric_limits<int>::min();
    for (auto const &hit : wire_data.hits) {
        if (!hit)
            continue;
        min_tick = std::min(min_tick, hit->StartTick());
        max_tick = std::max(max_tick, hit->EndTick());
    }
    if (min_tick >= max_tick)
        return;

    const std::size_t tick_count =
        static_cast<std::size_t>(max_tick - min_tick);
    std::vector<unsigned char> active(tick_count, 0U);
    std::vector<unsigned char> in_slice(tick_count, 0U);
    std::vector<std::vector<std::size_t>> hit_keys(tick_count);

    int const empty_label = static_cast<int>(
        image::SemanticClassifier::SemanticLabel::Empty);
    int const ambiguous_label = static_cast<int>(
        image::SemanticClassifier::SemanticLabel::Ambiguous);
    std::vector<int> semantic;
    if (context.hasSemantics())
        semantic.assign(tick_count, empty_label);

    for (auto const &hit : wire_data.hits) {
        if (!hit)
            continue;

        bool const slice_hit = context.isSliceHit(hit);
        auto const label = semanticLabel(hit, context);
        int const tick_start =
            std::max(hit->StartTick(), min_tick);
        int const tick_end =
            std::min(hit->EndTick(), max_tick);
        if (tick_start >= tick_end)
            continue;

        for (int tick = tick_start; tick < tick_end; ++tick) {
            std::size_t const index =
                static_cast<std::size_t>(tick - min_tick);
            active[index] = 1U;
            if (slice_hit)
                in_slice[index] = 1U;
            hit_keys[index].push_back(
                static_cast<std::size_t>(hit.key()));

            if (context.hasSemantics()) {
                int &slot = semantic[index];
                int const next = static_cast<int>(label);
                if (slot == empty_label)
                    slot = next;
                else if (slot != next)
                    slot = ambiguous_label;
            }
        }
    }

    for (auto const &range : wire.SignalROI().get_ranges()) {
        int const range_begin = range.begin_index();
        auto const &adcs = range.data();
        int const range_end =
            range_begin + static_cast<int>(adcs.size());
        int const sample_begin = std::max(range_begin, min_tick);
        int const sample_end = std::min(range_end, max_tick);
        if (sample_begin >= sample_end)
            continue;

        for (int tick = sample_begin; tick < sample_end; ++tick) {
            std::size_t const tick_index =
                static_cast<std::size_t>(tick - min_tick);
            if (active[tick_index] == 0U)
                continue;

            float const adc = adcs[tick - range_begin];
            if (!std::isfinite(adc) ||
                !(adc > kAdcThreshold))
                continue;

            double const drift =
                detector_properties.ConvertTicksToX(
                    static_cast<double>(tick),
                    wire_data.plane_id);
            if (!std::isfinite(drift))
                continue;

            for (auto const &destination : destinations) {
                auto const row =
                    destination.window
                        ->properties_[destination.plane_index]
                        .row(drift);
                if (!row)
                    continue;

                destination.window
                    ->detector_images_[destination.plane_index]
                    .set(*row, destination.column, adc, true);
                if (in_slice[tick_index] != 0U) {
                    destination.window
                        ->slice_images_[destination.plane_index]
                        .set(*row, destination.column,
                             static_cast<uint8_t>(1U), false);
                }
                if (context.hasSemantics() &&
                    semantic[tick_index] != empty_label) {
                    mergeSemanticPixel(
                        destination.window
                            ->semantic_images_[
                                destination.plane_index],
                        *row, destination.column,
                        static_cast<
                            image::SemanticClassifier::SemanticLabel>(
                            semantic[tick_index]));
                }

                auto &contributors =
                    destination.window
                        ->contributing_hit_keys_[
                            destination.plane_index];
                contributors.insert(
                    hit_keys[tick_index].begin(),
                    hit_keys[tick_index].end());
            }
        }
    }
}

std::vector<ImageFeatures>
SparseImagePacker::pack(
    const RasterizedImageWindow &window,
    const std::optional<TVector3> &vertex,
    bool include_semantics) {
    auto const &properties = window.properties();
    auto const &detector = window.detectorImages();
    auto const &semantic = window.semanticImages();
    auto const &slice = window.sliceImages();

    if (properties.size() != detector.size() ||
        properties.size() != semantic.size() ||
        properties.size() != slice.size()) {
        throw cet::exception("ImageProduction")
            << "Rasterized image components are not aligned.";
    }

    std::vector<ImageFeatures> output;
    output.reserve(properties.size());
    for (std::size_t i = 0U; i < properties.size(); ++i) {
        output.push_back(packPlane(
            detector[i], semantic[i], slice[i], properties[i],
            vertex, include_semantics,
            window.contributingHitCount(i)));
    }
    return output;
}

std::pair<int32_t, int32_t>
SparseImagePacker::vertexPixel(
    const ImageProperties &property,
    const std::optional<TVector3> &vertex) {
    if (!vertex)
        return {-1, -1};

    auto const projected = projectToImageView(
        static_cast<float>(vertex->X()),
        static_cast<float>(vertex->Y()),
        static_cast<float>(vertex->Z()), property.view());
    auto const column = property.col(projected.Z());
    auto const row = property.row(projected.X());
    if (!column || !row)
        return {-1, -1};
    return {static_cast<int32_t>(*row),
            static_cast<int32_t>(*column)};
}

ImageFeatures SparseImagePacker::packPlane(
    const Image<float> &detector,
    const Image<int> &semantic,
    const Image<uint8_t> &slice,
    const ImageProperties &property,
    const std::optional<TVector3> &vertex,
    bool include_semantics,
    std::size_t contributing_hits) {
    if (property.width() >
            std::numeric_limits<uint32_t>::max() ||
        property.height() >
            std::numeric_limits<uint32_t>::max()) {
        throw cet::exception("ImageProduction")
            << "Image dimensions exceed uint32_t.";
    }
    if (contributing_hits >
        std::numeric_limits<uint32_t>::max()) {
        throw cet::exception("ImageProduction")
            << "Contributing hit count exceeds uint32_t.";
    }

    ImageFeatures output;
    output.view = static_cast<int>(
        canonicalImageView(property.view()));
    output.width = static_cast<uint32_t>(property.width());
    output.height = static_cast<uint32_t>(property.height());
    output.hit_count =
        static_cast<uint32_t>(contributing_hits);
    output.origin_x =
        static_cast<float>(property.origin_x());
    output.origin_y =
        static_cast<float>(property.origin_y());
    output.pixel_w =
        static_cast<float>(property.pixel_w());
    output.pixel_h =
        static_cast<float>(property.pixel_h());
    auto const vertex_pixel = vertexPixel(property, vertex);
    output.vertex_row = vertex_pixel.first;
    output.vertex_col = vertex_pixel.second;
    output.feature_dim = kImageFeatureDimension;

    auto const &adc_pixels = detector.data();
    auto const &semantic_pixels = semantic.data();
    auto const &slice_pixels = slice.data();
    std::size_t active_pixels = 0U;
    for (float adc : adc_pixels) {
        if (std::isfinite(adc) && adc > 0.f)
            ++active_pixels;
    }

    output.coords.reserve(active_pixels * 2U);
    output.features.reserve(
        active_pixels * output.feature_dim);
    if (include_semantics)
        output.semantic.reserve(active_pixels);

    for (std::size_t index = 0U;
         index < adc_pixels.size(); ++index) {
        float const adc = adc_pixels[index];
        if (!std::isfinite(adc) || !(adc > 0.f))
            continue;

        output.coords.push_back(static_cast<int32_t>(
            index / property.width()));
        output.coords.push_back(static_cast<int32_t>(
            index % property.width()));
        output.features.push_back(adc);
        output.features.push_back(
            static_cast<float>(slice_pixels[index]));
        if (include_semantics) {
            output.semantic.push_back(static_cast<uint8_t>(
                semantic_pixels[index]));
        }
    }
    return output;
}

ImageProduction::ImageProduction(
    geo::GeometryCore const &geometry,
    PixelImageOptions const &options)
    : geometry_{&geometry}, options_{options}
{}

void ImageProduction::build(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<art::Ptr<recob::Hit>> &slice_hits,
    const std::vector<ImageProperties> &properties,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images,
    std::vector<Image<uint8_t>> &slice_images,
    detinfo::DetectorProperties const *detector_properties) const {
    ImageEventContext context(
        event, hits, slice_hits, options_);
    RasterizedImageWindow window(properties);
    ImageRasterizer(*geometry_).rasterize(
        context, {&window}, detector_properties);
    detector_images = window.detectorImages();
    semantic_images = window.semanticImages();
    slice_images = window.sliceImages();
}

} // namespace image
