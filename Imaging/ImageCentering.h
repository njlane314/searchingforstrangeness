#ifndef IMAGING_IMAGECENTERING_H
#define IMAGING_IMAGECENTERING_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include <TVector3.h>

namespace image {

inline TVector3
trimmedCentroid3D(const std::vector<art::Ptr<recob::SpacePoint>> &sps,
                  const std::vector<double> &weights,
                  const TVector3 &seed,
                  double radius)
{
    const double R2 = radius * radius;

    TVector3 sum(0., 0., 0.);
    double   W   = 0.0;

    const std::size_t n = std::min(sps.size(), weights.size());

    for (std::size_t i = 0; i < n; ++i) {
        auto const &sp = sps[i];
        if (!sp) continue;

        double w = weights[i];
        if (!(w > 0.0)) continue;

        auto const *xyz = sp->XYZ();
        TVector3 p(xyz[0], xyz[1], xyz[2]);

        if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
            continue;

        if ((p - seed).Mag2() > R2)
            continue;

        sum += w * p;
        W   += w;
    }

    if (W == 0.0) return seed;
    return sum * (1.0 / W);
}

}

#endif
