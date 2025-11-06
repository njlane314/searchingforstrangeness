#pragma once
#include <algorithm>
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

inline double median3(double a, double b, double c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return b;
}

inline double median2(double a, double b) { return 0.5 * (a + b); }

struct FusedAnchor {
    double x_star;                  
    double Y_star, Z_star;          
    double wU_star, wV_star, wW_star; 
};

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

inline FusedAnchor fuse_and_project(double xU, double wU,
                                    double xV, double wV,
                                    double xW, double wW)
{
    const double x_star = median3(xU, xV, xW);
    const double Z_star = median2(wW, (wU + wV));
    const double Y_star = (wV - wU) / std::sqrt(3.0);

    const double wU_star = 0.5 * Z_star - (std::sqrt(3.0)/2.0) * Y_star;
    const double wV_star = 0.5 * Z_star + (std::sqrt(3.0)/2.0) * Y_star;
    const double wW_star = Z_star;

    return {x_star, Y_star, Z_star, wU_star, wV_star, wW_star};
}

} 
