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
                     double min_charge /* ADC integral threshold */ = 0.0,
                     double trim_factor = 2.0 /* unused, kept for API compat */)
{
    (void)trim_factor;

    struct Sample {
        TVector3 pos;
        double   w;
    };

    std::vector<Sample> pts;
    pts.reserve(hits.size());

    for (auto const &h : hits) {
        if (!h) continue;

        double q = std::max(0.f, h->Integral());
        if (q <= min_charge) continue;

        auto view = common::GetPandoraView(h);
        TVector3 p = common::GetPandoraHitPosition(event, h, view);

        if (!std::isfinite(p.X()) || !std::isfinite(p.Y()) || !std::isfinite(p.Z()))
            continue;

        pts.push_back(Sample{p, q});
    }

    if (pts.empty()) {
        return TVector3(std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN());
    }

    TVector3 m(0., 0., 0.);
    double   W = 0.0;
    for (auto const &s : pts) {
        m += s.w * s.pos;
        W += s.w;
    }
    if (W > 0.0) m *= (1.0 / W);

    constexpr int    kMaxIter = 25;
    constexpr double kTol     = 0.1;

    for (int iter = 0; iter < kMaxIter; ++iter) {
        TVector3 num(0., 0., 0.);
        double   denom = 0.0;
        bool     singular = false;
        TVector3 singular_pos;

        for (auto const &s : pts) {
            TVector3 diff = m - s.pos;
            double   d    = diff.Mag();

            if (d < 1e-3) {
                singular      = true;
                singular_pos  = s.pos;
                break;
            }

            double w_over_d = s.w / d;
            num   += w_over_d * s.pos;
            denom += w_over_d;
        }

        if (singular || denom == 0.0) {
            if (singular) m = singular_pos;
            break;
        }

        TVector3 m_new = num * (1.0 / denom);
        if ((m_new - m).Mag() < kTol) {
            m = m_new;
            break;
        }
        m = m_new;
    }

    return m;
}

}

#endif
