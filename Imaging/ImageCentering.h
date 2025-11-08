#ifndef IMAGING_IMAGECENTERING_H
#define IMAGING_IMAGECENTERING_H
#include <cmath>
#include <set>
#include <utility>
#include <vector>

#include "Common/PandoraUtilities.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"

#include <TVector3.h>

namespace image {
class ImageCentering {
  public:
    static std::pair<double, double>
    centroidWithinRadius(const art::Event &event, common::PandoraView view,
                         const std::vector<art::Ptr<recob::Hit>> &hits,
                         double radius,
                         const std::set<unsigned int> &bad_channels,
                         double vtx_z, double vtx_x) {
        double W = 0.0;
        double Zs = 0.0;
        double Xs = 0.0;
        for (auto const &h : hits) {
            if (bad_channels.count(h->Channel()))
                continue;
            if (common::GetPandoraView(h) != view)
                continue;
            double q = std::max(0.f, h->Integral());
            if (q <= 0)
                continue;
            TVector3 p = common::GetPandoraHitPosition(event, h, view);
            double d = std::hypot(p.Z() - vtx_z, p.X() - vtx_x);
            if (d <= radius) {
                W += q;
                Zs += q * p.Z();
                Xs += q * p.X();
            }
        }
        if (W == 0.0)
            return {vtx_z, vtx_x};
        return {Zs / W, Xs / W};
    }
};
}

#endif // IMAGING_IMAGECENTERING_H
