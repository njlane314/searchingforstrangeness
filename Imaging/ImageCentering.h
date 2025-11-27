#ifndef IMAGING_IMAGECENTERING_H
#define IMAGING_IMAGECENTERING_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include "Common/PandoraUtilities.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"

#include <TVector3.h>

namespace image {

inline std::pair<double, double>
centroidWithinRadius(const art::Event &event, common::PandoraView view,
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     double radius,
                     double vtx_z, double vtx_x) {
    double W = 0.0;
    double Zs = 0.0;
    double Xs = 0.0;
    for (auto const &h : hits) {
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

inline TVector3
robustGlobalCentroid(const art::Event &event,
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     double min_charge = 0.0,
                     double trim_factor = 2.0)
{
    struct HitPoint {
        TVector3 pos;
        double   q;
    };

    std::vector<HitPoint> hp;
    hp.reserve(hits.size());

    for (auto const &h : hits) {
        if (!h) continue;

        double q = std::max(0.f, h->Integral());
        if (q <= min_charge) continue;

        auto view = common::GetPandoraView(h);
        TVector3 p = common::GetPandoraHitPosition(event, h, view);

        if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
            continue;

        hp.push_back({p, q});
    }

    if (hp.empty()) {
        return TVector3(std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN());
    }

    double   sumw = 0.0;
    TVector3 c1(0., 0., 0.);
    for (auto const &hpt : hp) {
        sumw += hpt.q;
        c1   += hpt.q * hpt.pos;
    }
    if (sumw > 0.0) c1 *= (1.0 / sumw);

    std::vector<double> dists;
    dists.reserve(hp.size());
    for (auto const &hpt : hp) {
        dists.push_back((hpt.pos - c1).Mag());
    }

    auto mid = dists.begin() + dists.size() / 2;
    std::nth_element(dists.begin(), mid, dists.end());
    double med = *mid;
    if (med <= 0.0) {
        return c1;
    }

    const double R = trim_factor * med;

    sumw = 0.0;
    TVector3 c2(0., 0., 0.);
    for (auto const &hpt : hp) {
        if ((hpt.pos - c1).Mag() > R) continue;
        sumw += hpt.q;
        c2   += hpt.q * hpt.pos;
    }

    if (sumw <= 0.0)
        return c1;

    c2 *= (1.0 / sumw);
    return c2;
}

}

#endif
