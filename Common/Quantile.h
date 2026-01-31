#ifndef COMMON_QUANTILE_H
#define COMMON_QUANTILE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace common {

template <typename T>
inline T clamp(T v, T lo, T hi) { return std::max(lo, std::min(v, hi)); }

inline float quantile_linear(std::vector<float> v, double p) {
    if (v.empty())
        return std::numeric_limits<float>::quiet_NaN();
    p = clamp(p, 0.0, 1.0);
    double pos = p * (v.size() - 1);
    size_t lo = static_cast<size_t>(std::floor(pos));
    size_t hi = static_cast<size_t>(std::ceil(pos));
    std::nth_element(v.begin(), v.begin() + lo, v.end());
    float qlo = v[lo];
    if (hi == lo)
        return qlo;
    std::nth_element(v.begin(), v.begin() + hi, v.end());
    float qhi = v[hi];
    double frac = pos - lo;
    return static_cast<float>((1.0 - frac) * qlo + frac * qhi);
}

} // namespace common

#endif // COMMON_QUANTILE_H
