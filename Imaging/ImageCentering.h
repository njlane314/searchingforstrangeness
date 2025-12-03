#ifndef IMAGING_IMAGECENTERING_H
#define IMAGING_IMAGECENTERING_H

#include <cmath>
#include <vector>

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include <TVector3.h>

namespace image {

inline TVector3
trimmedMeanCenterFromSpacePoints(
    const std::vector<art::Ptr<recob::SpacePoint>> &sps,
    const TVector3 &seed,
    double radius)
{
    const double R2 = radius * radius;

    TVector3 sum(0., 0., 0.);
    double   W = 0.0;

    for (auto const &sp : sps) {
        if (!sp) continue;

        auto const *xyz = sp->XYZ();
        TVector3 p(xyz[0], xyz[1], xyz[2]);

        if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
            continue;

        if ((p - seed).Mag2() > R2)
            continue;

        double w = 1.0;

        sum += w * p;
        W   += w;
    }

    if (W == 0.0) return seed;
    return sum * (1.0 / W);
}

}

#endif
