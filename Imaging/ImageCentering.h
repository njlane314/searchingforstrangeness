#pragma once
#include <algorithm>
#include <cmath>

namespace image {

// Median of three (robust to one outlier)
inline double median3(double a, double b, double c) {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
    return b;
}

// "Median" of two is their average in common practice
inline double median2(double a, double b) { return 0.5 * (a + b); }

struct FusedAnchor {
    double x_star;                  // X*
    double Y_star, Z_star;          // Y*, Z*
    double wU_star, wV_star, wW_star; // projected to U,V,W
};

// Implements the equations from your spec:
//   x*  = median(x_U, x_V, x_W)
//   Z*  = median(w_W, w_U + w_V)
//   Y*  = (w_V - w_U) / sqrt(3)
// and
//   w_U* = 0.5 Z* - (sqrt(3)/2) Y*
//   w_V* = 0.5 Z* + (sqrt(3)/2) Y*
//   w_W* = Z*
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

} // namespace image
