#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H

#include <vector>
#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <TFile.h>
#include <TTree.h>
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "canvas/Utilities/InputTag.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "TDirectoryFile.h"
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/WireReadout.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"

namespace image 
{
    class ImageProperties {
    public:
        ImageProperties() = default;
        ImageProperties(double center_x, double center_y, size_t width, size_t height, double pixel_h, double pixel_w, geo::View_t view)
            : center_x_(center_x), center_y_(center_y), height_(height), width_(width), pixel_w_(pixel_w), pixel_h_(pixel_h), view_(view) {
            origin_x_ = center_x - (width * pixel_w) / 2.0;
            origin_y_ = center_y - (height * pixel_h) / 2.0;
        }
        size_t index(size_t row, size_t col) const {
            if (row >= height_ || col >= width_) return static_cast<size_t>(-1);
            return col * height_ + row;
        }
        size_t col(double x) const {
            if (x < origin_x_ || x >= origin_x_ + width_ * pixel_w_) return static_cast<size_t>(-1);
            return static_cast<size_t>((x - origin_x_) / pixel_w_);
        }
        size_t row(double y) const {
            if (y < origin_y_ || y >= origin_y_ + height_ * pixel_h_) return static_cast<size_t>(-1);
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
        geo::View_t view_;
    };

    class Image {
    public:
        Image() = default;
        Image(const ImageProperties& prop)
            : prop_(prop), pixels_(prop.height() * prop.width(), 0.0f) {}
        void set(size_t row, size_t col, float value, bool accumulate = true) {
            size_t idx = prop_.index(row, col);
            if (idx != static_cast<size_t>(-1)) {
                if (accumulate)
                    pixels_[idx] += value;
                else
                    pixels_[idx] = value;
            }
        }
        float get(size_t row, size_t col) const {
            size_t idx = prop_.index(row, col);
            return (idx != static_cast<size_t>(-1)) ? pixels_[idx] : 0.0f;
        }
        void clear(float value = 0.0f) {
            std::fill(pixels_.begin(), pixels_.end(), value);
        }
        std::vector<float> data() const {
            return pixels_;
        }
        geo::View_t view() const { return prop_.view(); }
        size_t height() const { return prop_.height(); }
        size_t width() const { return prop_.width(); }
    private:
        ImageProperties prop_;
        std::vector<float> pixels_;
    };
} 

#endif