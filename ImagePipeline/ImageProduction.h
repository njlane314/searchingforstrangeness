#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

#include "ImagePipeline/Image.h"
#include "ImagePipeline/SemanticClassifier.h"
#include "ImagePipeline/ImageWindowGeometry.h"
#include "Products/ImageFeatures.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <TVector3.h>

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace image {

struct PixelImageOptions {
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
    } producers;
    image::SemanticClassifier *semantic{nullptr};
};

/// Event-scoped handles, associations, and truth lookup tables.
class ImageEventContext {
  public:
    using BacktrackerAssociations =
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>;

    ImageEventContext(
        const art::Event &event,
        const std::vector<art::Ptr<recob::Hit>> &hits,
        const std::vector<art::Ptr<recob::Hit>> &slice_hits,
        PixelImageOptions const &options);
    ~ImageEventContext();

    const std::vector<recob::Wire> &wires() const noexcept;
    const art::FindManyP<recob::Hit> &
    wireHitAssociations() const noexcept;
    bool acceptsHit(const art::Ptr<recob::Hit> &hit) const;
    bool isSliceHit(const art::Ptr<recob::Hit> &hit) const;
    bool hasSemantics() const noexcept;
    const BacktrackerAssociations *
    backtrackerAssociations() const noexcept;
    const std::map<int, std::size_t> &
    trackIdToIndex() const noexcept;
    const std::vector<image::SemanticClassifier::SemanticLabel> &
    semanticLabels() const noexcept;

  private:
    art::ValidHandle<std::vector<recob::Wire>> wires_;
    art::ValidHandle<std::vector<recob::Hit>> all_hits_;
    art::FindManyP<recob::Hit> wire_hit_assoc_;
    std::set<art::Ptr<recob::Hit>> accepted_hits_;
    std::set<std::size_t> slice_hit_keys_;
    art::Handle<std::vector<simb::MCParticle>> mcp_vector_;
    std::unique_ptr<BacktrackerAssociations> mcp_bkth_assoc_;
    std::map<int, std::size_t> trackid_to_index_;
    std::vector<image::SemanticClassifier::SemanticLabel>
        semantic_labels_;
};

/// Dense images and contributing-hit bookkeeping for one image window.
class RasterizedImageWindow {
  public:
    explicit RasterizedImageWindow(
        std::vector<ImageProperties> properties);

    const std::vector<ImageProperties> &
    properties() const noexcept;
    const std::vector<Image<float>> &
    detectorImages() const noexcept;
    const std::vector<Image<int>> &
    semanticImages() const noexcept;
    const std::vector<Image<uint8_t>> &
    sliceImages() const noexcept;
    std::size_t contributingHitCount(
        std::size_t plane_index) const;

  private:
    friend class ImageRasterizer;

    std::vector<ImageProperties> properties_;
    std::vector<Image<float>> detector_images_;
    std::vector<Image<int>> semantic_images_;
    std::vector<Image<uint8_t>> slice_images_;
    std::vector<std::set<std::size_t>> contributing_hit_keys_;
};

/// Rasterize any number of windows in one wire/ROI/tick traversal.
class ImageRasterizer {
  public:
    explicit ImageRasterizer(
        geo::GeometryCore const &geometry);

    void rasterize(
        const ImageEventContext &event_context,
        const std::vector<RasterizedImageWindow *> &windows,
        detinfo::DetectorProperties const *detector_properties) const;

  private:
    static constexpr float kAdcThreshold = 0.0f;

    struct WirePreparation;
    struct Destination;

    std::optional<WirePreparation>
    prepareWire(std::size_t wire_index,
                const ImageEventContext &context) const;
    static std::vector<Destination>
    destinationsForWire(
        geo::View_t view, double wire_coordinate,
        const std::vector<RasterizedImageWindow *> &windows);
    static image::SemanticClassifier::SemanticLabel
    semanticLabel(const art::Ptr<recob::Hit> &hit,
                  const ImageEventContext &context);
    static void mergeSemanticPixel(
        Image<int> &image, std::size_t row, std::size_t column,
        image::SemanticClassifier::SemanticLabel label);
    void fillWire(
        std::size_t wire_index,
        const ImageEventContext &context,
        const std::vector<RasterizedImageWindow *> &windows,
        detinfo::DetectorProperties const &detector_properties) const;

    geo::GeometryCore const *geometry_{nullptr};
};

/// Convert a rasterized window into the persisted sparse product.
class SparseImagePacker {
  public:
    static std::vector<ImageFeatures>
    pack(const RasterizedImageWindow &window,
         const std::optional<TVector3> &vertex,
         bool include_semantics);

  private:
    static std::pair<int32_t, int32_t>
    vertexPixel(const ImageProperties &property,
                const std::optional<TVector3> &vertex);
    static ImageFeatures
    packPlane(const Image<float> &detector,
              const Image<int> &semantic,
              const Image<uint8_t> &slice,
              const ImageProperties &property,
              const std::optional<TVector3> &vertex,
              bool include_semantics,
              std::size_t contributing_hits);
};

/// Backward-compatible single-window facade.
class ImageProduction {
  public:
    ImageProduction(geo::GeometryCore const &geometry,
                    PixelImageOptions const &options);

    void build(
        const art::Event &event,
        const std::vector<art::Ptr<recob::Hit>> &hits,
        const std::vector<art::Ptr<recob::Hit>> &slice_hits,
        const std::vector<ImageProperties> &properties,
        std::vector<Image<float>> &detector_images,
        std::vector<Image<int>> &semantic_images,
        std::vector<Image<uint8_t>> &slice_images,
        detinfo::DetectorProperties const *detector_properties) const;

  private:
    geo::GeometryCore const *geometry_{nullptr};
    PixelImageOptions options_;
};

} // namespace image

#endif
