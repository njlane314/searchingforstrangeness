#ifndef IMAGEPIPELINE_IMAGECENTRING_H
#define IMAGEPIPELINE_IMAGECENTRING_H

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include <TVector3.h>

#include <optional>
#include <vector>

namespace image {

enum class ImageCentreSeed {
    Vertex,
    WeightedCentroid,
    OriginFallback
};

struct ImageCentre {
    TVector3 position;
    ImageCentreSeed seed{ImageCentreSeed::OriginFallback};
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

class ImageCentring {
  public:
    ImageCentring(art::InputTag hit_producer,
                  art::InputTag spacepoint_producer);

    ImageCentre compute(
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
