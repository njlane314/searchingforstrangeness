#ifndef IMAGE_H
#define IMAGE_H

#include <algorithm>
#include <cstddef>
#include <optional>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace image {

class ImageProperties {
public:
    ImageProperties() = default;
    ImageProperties(double center_x, double center_y, size_t width, size_t height, double pixel_h, double pixel_w, geo::View_t view)
        : origin_x_(center_x - (width * pixel_w) / 2.0)
        , origin_y_(center_y - (height * pixel_h) / 2.0)
        , height_(height)
        , width_(width)
        , pixel_w_(pixel_w)
        , pixel_h_(pixel_h)
        , view_(view)
    {}

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
    double origin_x_{0.0};
    double origin_y_{0.0};
    size_t height_{0U};
    size_t width_{0U};
    double pixel_w_{0.0};
    double pixel_h_{0.0};
    geo::View_t view_{geo::kUnknown};
};

template <typename T>
class Image {
public:
    Image() = default;
    explicit Image(const ImageProperties &prop)
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
    const std::vector<T>& data() const noexcept { return pixels_; }

    geo::View_t view() const { return prop_.view(); }
    size_t height() const { return prop_.height(); }
    size_t width() const { return prop_.width(); }
    const ImageProperties& properties() const { return prop_; }

private:
    ImageProperties prop_;
    std::vector<T> pixels_;
};
}

#endif
