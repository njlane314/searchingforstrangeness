#ifndef IMAGE_H
#define IMAGE_H

#include <vector>
#include <string>
#include <map>
#include <array>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <limits>
#include <optional>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "CommonDefs/Pandora.h"

namespace analysis 
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
        geo::View_t view_ {geo::kUnknown};
    };

    template <typename T>
    class Image {
    public:
        Image() = default;
        Image(const ImageProperties& prop)
        : prop_(prop), pixels_(prop.height() * prop.width(), T(0)) {}
        void set(size_t row, size_t col, T value, bool accumulate = true) {
            size_t idx = prop_.index(row, col);
            if (idx != static_cast<size_t>(-1)) {
                if (accumulate)
                    pixels_[idx] += value;
                else
                    pixels_[idx] = value;
            }
        }
        T get(size_t row, size_t col) const {
            size_t idx = prop_.index(row, col);
            return (idx != static_cast<size_t>(-1)) ? pixels_[idx] : T(0);
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
    private:
        ImageProperties prop_;
        std::vector<T> pixels_;
    };

}

#endif