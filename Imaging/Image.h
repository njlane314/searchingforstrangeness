#ifndef IMAGE_H
#define IMAGE_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "canvas/Utilities/InputTag.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TDirectoryFile.h"
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace image {
class ImageProperties {
    public:
    ImageProperties() = default;
    ImageProperties(double center_x, double center_y, size_t width, size_t height, double pixel_h, double pixel_w, geo::View_t view)
        : center_x_(center_x), center_y_(center_y), height_(height), width_(width), pixel_w_(pixel_w), pixel_h_(pixel_h), view_(view) {
        origin_x_ = center_x - (width * pixel_w) / 2.0;
        origin_y_ = center_y - (height * pixel_h) / 2.0;
    }
    std::optional<size_t> index(size_t row, size_t col) const {
        if (row >= height_ || col >= width_)
            return std::nullopt;
        return row * width_ + col;
    }
    std::optional<size_t> col(double x) const {
        if (x < origin_x_ || x >= origin_x_ + width_ * pixel_w_)
            return std::nullopt;
        return static_cast<size_t>((x - origin_x_) / pixel_w_);
    }
    std::optional<size_t> row(double y) const {
        if (y < origin_y_ || y >= origin_y_ + height_ * pixel_h_)
            return std::nullopt;
        return static_cast<size_t>((y - origin_y_) / pixel_h_);
    }
    size_t height() const { return height_; }
    size_t width() const { return width_; }
    double pixel_w() const { return pixel_w_; }
    double pixel_h() const { return pixel_h_; }
    geo::View_t view() const { return view_; }
    double origin_x() const { return origin_x_; }
    double origin_y() const { return origin_y_; }
    double max_x() const { return origin_x_ + width_ * pixel_w_; }
    double max_y() const { return origin_y_ + height_ * pixel_h_; }

    private:
    double center_x_, center_y_;
    double origin_x_, origin_y_;
    size_t height_, width_;
    double pixel_w_, pixel_h_;
    geo::View_t view_{geo::kUnknown};
};

template <typename T>
class Image {
    public:
    Image() = default;
    Image(const ImageProperties &prop)
        : prop_(prop), pixels_(prop.height() * prop.width(), T(0)) {}
    void set(size_t row, size_t col, T value, bool accumulate = true) {
        auto idx = prop_.index(row, col);
        if (idx) {
            if (accumulate)
                pixels_[*idx] += value;
            else
                pixels_[*idx] = value;
        }
    }
    T get(size_t row, size_t col) const {
        auto idx = prop_.index(row, col);
        return idx ? pixels_[*idx] : T(0);
    }
    void clear(T value = T(0)) {
        std::fill(pixels_.begin(), pixels_.end(), value);
    }
    std::vector<T> data() const {
        return pixels_;
    }
    geo::View_t view() const { return prop_.view(); }
    size_t height() const { return prop_.height(); }
    size_t width() const { return prop_.width(); }
    const ImageProperties& properties() const { return prop_; }

    Image<T>& blur(float sigma_px);

    private:
    ImageProperties prop_;
    std::vector<T> pixels_;
};

template <typename T>
inline Image<T>& Image<T>::blur(float sigma_px)
{
    if (!(sigma_px > 0.f)) return *this;

    const int R = std::max(1, static_cast<int>(std::ceil(3.0f * sigma_px)));
    std::vector<float> k(2 * R + 1);
    double sum = 0.0;
    for (int i = -R; i <= R; ++i) {
        const double v = std::exp(-(1.0 * i * i) / (2.0 * sigma_px * sigma_px));
        k[i + R] = static_cast<float>(v);
        sum += v;
    }
    const float inv = static_cast<float>(1.0 / sum);
    for (auto& v : k) v *= inv;

    const size_t H = height();
    const size_t W = width();

    Image<T> tmp(properties());
    tmp.clear(T{});
    for (size_t r = 0; r < H; ++r) {
        for (size_t c = 0; c < W; ++c) {
            double acc = 0.0;
            for (int dr = -R; dr <= R; ++dr) {
                const int rr = std::clamp<int>(static_cast<int>(r) + dr, 0, static_cast<int>(H) - 1);
                acc += static_cast<double>(k[dr + R]) *
                       static_cast<double>(get(static_cast<size_t>(rr), c));
            }
            tmp.set(r, c, static_cast<T>(acc), false);
        }
    }

    Image<T> out(properties());
    out.clear(T{});
    for (size_t r = 0; r < H; ++r) {
        for (size_t c = 0; c < W; ++c) {
            double acc = 0.0;
            for (int dc = -R; dc <= R; ++dc) {
                const int cc = std::clamp<int>(static_cast<int>(c) + dc, 0, static_cast<int>(W) - 1);
                acc += static_cast<double>(k[dc + R]) *
                       static_cast<double>(tmp.get(r, static_cast<size_t>(cc)));
            }
            out.set(r, c, static_cast<T>(acc), false);
        }
    }

    *this = std::move(out);
    return *this;
}
} 

#endif 
