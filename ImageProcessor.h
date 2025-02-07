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
#include "canvas/Utilities/InputTag.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TDirectoryFile.h"
#include "EventClassifier.h"

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

std::vector<Image> constructInputImages(const std::vector<ImageProperties>& properties, const std::vector<art::Ptr<recob::Wire>>& wires, const geo::GeometryCore& geo) {
    std::vector<Image> input;
    for (const auto& prop : properties) input.emplace_back(prop);

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

                for (auto& img : input) {
                    if (img.properties().view() == view) {
                        size_t row = img.properties().row(x);
                        size_t col = img.properties().col(coord);
                        img.set(row, col, adcs[idx]);
                    }
                }
            }
        }
    }
    return input;
}

std::vector<Image> constructTruthImages(const std::vector<ImageProperties>& properties, 
                                          const std::vector<art::Ptr<recob::Hit>>& hits, 
                                          const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
                                          const size_t kernel_size, 
                                          const signature::Pattern& pattern, 
                                          const geo::GeometryCore& geo) {
    std::vector<Image> truth;
    truth.reserve(properties.size());
    for (const auto& prop : properties) {
        Image img(prop);
        img.clear(0.0f);
        truth.push_back(std::move(img));
    }

    auto const* detProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    std::unordered_map<int, simb::MCParticle> truth_particles;
    for (const auto& [type, signature] : pattern) {
        for (const auto& particle : signature) {
            truth_particles[particle->TrackId()] = *particle;
        }
    }

    constexpr float angle = 1.04719758034f; // ~pi/3
    auto computeCoordinate = [angle](geo::View_t view, const TVector3& wire_center) -> float {
        if (view == geo::kW) {
            return wire_center.Z();
        }
        const float sign = (view == geo::kU) ? 1.0f : -1.0f;
        return (wire_center.Z() * std::cos(sign * angle)) - (wire_center.Y() * std::sin(sign * angle));
    };

    auto applyKernel = [kernel_size](Image& img, size_t row, size_t col, float value) {
        const int half = static_cast<int>(kernel_size) / 2;
        for (int dr = -half; dr <= half; ++dr) {
            for (int dc = -half; dc <= half; ++dc) {
                img.set(row + dr, col + dc, value, /*accumulate=*/false);
            }
        }
    };

    for (const auto& hit : hits) {
        const auto& assocParticles = mcp_bkth_assoc.at(hit.key());
        const auto& assocData = mcp_bkth_assoc.data(hit.key());
        if (assocParticles.empty() || assocData.empty())
            continue;

        int track_id = -1;
        bool valid_hit = false;
        for (size_t i = 0; i < assocData.size(); ++i) {
            if (assocData[i]->isMaxIDE == 1) {
                track_id = assocParticles[i]->TrackId();
                valid_hit = true;
                break;
            }
        }
        if (!valid_hit || truth_particles.find(track_id) == truth_particles.end())
            continue;

        double x = detProp->ConvertTicksToX(hit->PeakTime(), hit->WireID());
        TVector3 wire_center = geo.Cryostat(hit->WireID().Cryostat)
                                    .TPC(hit->WireID().TPC)
                                    .Plane(hit->WireID().Plane)
                                    .Wire(hit->WireID().Wire)
                                    .GetCenter();
        geo::View_t view = hit->View();
        float transformed_coord = computeCoordinate(view, wire_center);

        for (auto& img : truth) {
            if (img.properties().view() != view)
                continue;

            size_t row = img.properties().row(x);
            size_t col = img.properties().col(transformed_coord);
            if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1))
                continue;
            applyKernel(img, row, col, 1.0f);
        }
    }

    return truth;
}

class ImageTrainingHandler {
public:
    ImageTrainingHandler(TTree* tree, const geo::GeometryCore& geo)
        : tree_(tree), geo_(geo) {
        tree_->Branch("run", &run_);
        tree_->Branch("subrun", &subrun_);
        tree_->Branch("event", &event_);
        tree_->Branch("event_type", &event_type_);
        tree_->Branch("planes", &planes_);
        tree_->Branch("width", &width_);
        tree_->Branch("height", &height_);
        tree_->Branch("input_data", &input_data_);
        tree_->Branch("truth_data", &truth_data_);
    }

    void reset() {
        run_ = subrun_ = event_ = 0;
        event_type_ = static_cast<int>(signature::EventType::kOther);
        planes_.clear();
        width_.clear();
        height_.clear();
        input_data_.clear();
        truth_data_.clear();
    }

    void add(const art::Event& e, 
             const signature::EventType event_type,
             const std::vector<art::Ptr<recob::Wire>>& wires,
             const std::vector<art::Ptr<recob::Hit>>& hits,
             const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
             const signature::Pattern& pattern,
             const size_t& kernel_size,
             const std::vector<ImageProperties>& properties) {
        this->reset();

        run_ = e.run();
        subrun_ = e.subRun();
        event_ = e.event();

        event_type_ = static_cast<int>(event_type);

        std::vector<Image> input = constructInputImages(properties, wires, geo_);
        std::vector<Image> truth = constructTruthImages(properties, hits, mcp_bkth_assoc, kernel_size, pattern, geo_);

        std::vector<geo::View_t> ordered_views = { geo::kU, geo::kV, geo::kW };

        for (auto view : ordered_views) {
            auto input_iter = std::find_if(input.begin(), input.end(),
                                        [view](const Image& img) {
                                            return img.properties().view() == view;
                                        });

            auto truth_iter = std::find_if(truth.begin(), truth.end(),
                                        [view](const Image& img) {
                                            return img.properties().view() == view;
                                        });

            if (input_iter != input.end() && truth_iter != truth.end()) {
                if (input_iter->properties().width() != truth_iter->properties().width() || input_iter->properties().height() != truth_iter->properties().height())
                    continue;

                planes_.push_back(static_cast<int>(view));
                width_.push_back(input_iter->properties().width());
                height_.push_back(input_iter->properties().height());
                input_data_.push_back(input_iter->data());
                truth_data_.push_back(truth_iter->data());
            }
        }
        
        tree_->Fill();
    }

private:
    const geo::GeometryCore& geo_;

    TTree* tree_;
    int run_, subrun_, event_;
    int event_type_;
    std::vector<int> planes_;
    std::vector<int> width_; 
    std::vector<int> height_;
    std::vector<std::vector<float>> input_data_;
    std::vector<std::vector<float>> truth_data_; 
};

} 

#endif // IMAGEPROCESSOR_H