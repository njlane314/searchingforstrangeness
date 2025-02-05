#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H

#include <torch/torch.h>
#include <torch/script.h>

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

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TDirectoryFile.h"

namespace image {

class ImageProperties {
public:
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

    const ImageProperties& properties() const { return prop_; }

    torch::Tensor tensor() const {
        return torch::from_blob(const_cast<float*>(pixels_.data()),
                                {1, 1, (long)prop_.height(), (long)prop_.width()},
                                torch::kFloat32).clone();
    }

private:
    ImageProperties prop_;
    std::vector<float> pixels_;
};

std::vector<Image> ConvertWiresToImages(const std::vector<ImageProperties>& properties, const std::vector<art::Ptr<recob::Wire>>& wires, const geo::GeometryCore& geo) {
    std::vector<Image> images;
    for (const auto& prop : properties) images.emplace_back(prop);

    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    for (const auto& wire : wires) {
        auto ch_id = wire->Channel();
        auto wire_ids = geo.ChannelToWire(ch_id);
        if (wire_ids.empty()) continue;

        auto view = geo.View(wire_ids.front());
        for (const auto& range : wire->SignalROI().get_ranges()) {
            const auto& adcs = range.data();
            int start_idx = range.begin_index();

            for (size_t idx = 0; idx < adcs.size(); ++idx) {
                geo::PlaneID plane_id(wire_ids.front());
                int tdc = start_idx + idx;
                double x = detProp->ConvertTicksToX(tdc, plane_id);
                TVector3 wire_center = geo.Cryostat(wire_ids.front().Cryostat).TPC(wire_ids.front().TPC).Plane(wire_ids.front().Plane).Wire(wire_ids.front().Wire).GetCenter();

                float coord = (view == geo::kW) ? wire_center.Z() :
                              (wire_center.Z() * std::cos((view == geo::kU ? 1 : -1) * 1.04719758034)) - 
                              (wire_center.Y() * std::sin((view == geo::kU ? 1 : -1) * 1.04719758034));

                for (auto& img : images) {
                    if (img.properties().view() == view) {
                        size_t row = img.properties().row(x);
                        size_t col = img.properties().col(coord);
                        img.set(row, col, adcs[idx]);
                    }
                }
            }
        }
    }
    return images;
}

class ImageTrainingHandler {
public:
    ImageTrainingHandler(TTree* tree, const geo::GeometryCore& geo)
        : tree_(tree), geo_(geo) {
        tree_->Branch("run", &run_);
        tree_->Branch("subrun", &subrun_);
        tree_->Branch("event", &event_);
        tree_->Branch("is_signal", &is_signal_);
        tree_->Branch("planes", &planes_);
        tree_->Branch("image_data", &image_data_);
    }

    void reset() {
        run_ = subrun_ = event_ = 0;
        is_signal_ = false;
        planes_.clear();
        image_data_.clear();
    }

    void add(const art::Event& e, 
             const bool is_signal, 
             const std::vector<art::Ptr<recob::Wire>>& wires,
             const std::vector<ImageProperties>& properties) {
        reset();

        run_ = e.run();
        subrun_ = e.subRun();
        event_ = e.event();

        is_signal_ = is_signal;

        auto images = ConvertWiresToImages(properties, wires, geo_);
        for (const auto& img : images) {
            planes_.push_back(static_cast<int>(img.properties().view()));
            auto pixels = img.data();
            image_data_.insert(image_data_.end(), pixels.begin(), pixels.end());
        }

        tree_->Fill();
    }

private:
    const geo::GeometryCore& geo_;

    TTree* tree_;
    int run_, subrun_, event_;
    bool is_signal_;
    std::vector<int> planes_;
    std::vector<float> image_data_;
};

} 

#endif // IMAGEPROCESSOR_H