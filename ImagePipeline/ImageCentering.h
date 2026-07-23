#ifndef IMAGEPIPELINE_IMAGECENTERING_H
#define IMAGEPIPELINE_IMAGECENTERING_H

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include <TVector3.h>

#include <optional>
#include <vector>

namespace image {

enum class ImageCenterSeed {
    Vertex,
    WeightedCentroid,
    OriginFallback
};

struct ImageCenter {
    TVector3 position;
    ImageCenterSeed seed{ImageCenterSeed::OriginFallback};
};

std::optional<TVector3>
weightedCentroid(
    const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
    const std::vector<double> &weights);

TVector3 trimmedCentroid(
    const std::vector<art::Ptr<recob::SpacePoint>> &spacepoints,
    const std::vector<double> &weights,
    const TVector3 &seed,
    double radius);

class ImageCentering {
  public:
    ImageCentering(art::InputTag hit_producer,
                   art::InputTag spacepoint_producer);

    ImageCenter compute(
        const art::Event &event,
        const std::vector<art::Ptr<recob::Hit>> &neutrino_hits,
        const std::optional<TVector3> &vertex,
        double trimming_radius) const;

  private:
    art::InputTag hit_producer_;
    art::InputTag spacepoint_producer_;
};

} // namespace image

#endif
