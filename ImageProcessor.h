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
    ImageProperties(double center_x, double center_y, size_t width, size_t height, double pixel_h, double pixel_w, geo::View_t view, int kernel)
        : center_x_(center_x), center_y_(center_y), height_(height), width_(width), pixel_w_(pixel_w), pixel_h_(pixel_h), view_(view), kernel_(kernel) {
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
    int kernel() const { return kernel_; }
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
    int kernel_; 
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

std::vector<Image> constructInputImages(const std::vector<ImageProperties>& properties,
                                          const std::vector<art::Ptr<recob::Wire>>& wires,
                                          const geo::GeometryCore& geo) {
    std::vector<Image> input;
    for (const auto& prop : properties) input.emplace_back(prop);
    auto const* det_props = lar::providerFrom<detinfo::DetectorPropertiesService>();
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
                double x = det_props->ConvertTicksToX(tdc, plane_id);
                TVector3 wire_center = geo.Cryostat(wire_ids.front().Cryostat)
                                            .TPC(wire_ids.front().TPC)
                                            .Plane(wire_ids.front().Plane)
                                            .Wire(wire_ids.front().Wire)
                                            .GetCenter();
                float coord = (view == geo::kW) ? wire_center.Z() :
                              (wire_center.Z() * std::cos((view == geo::kU ? 1 : -1) * 1.04719758034))
                              - (wire_center.Y() * std::sin((view == geo::kU ? 1 : -1) * 1.04719758034));
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
                                          const signature::Pattern& pattern,
                                          const signature::EventClassifier& classifier,
                                          const geo::GeometryCore& geo) {
    std::vector<Image> truth;
    truth.reserve(properties.size());
    for (const auto& prop : properties) {
        Image img(prop);
        img.clear(signature::kEmptySignature);
        truth.push_back(std::move(img));
    }
    auto const* det_props = lar::providerFrom<detinfo::DetectorPropertiesService>();
    std::unordered_map<int, signature::SignatureType> truth_particles_U;
    std::unordered_map<int, signature::SignatureType> truth_particles_V;
    std::unordered_map<int, signature::SignatureType> truth_particles_W;
    for (const auto& [type, sig] : pattern) {
        if (classifier.passClarity(sig, common::TPC_VIEW_U))
            for (const auto& particle : sig)
                truth_particles_U[particle->TrackId()] = type;
        if (classifier.passClarity(sig, common::TPC_VIEW_V))
            for (const auto& particle : sig)
                truth_particles_V[particle->TrackId()] = type;
        if (classifier.passClarity(sig, common::TPC_VIEW_W))
            for (const auto& particle : sig)
                truth_particles_W[particle->TrackId()] = type;
    }
    constexpr float angle = 1.04719758034f; // ~pi/3
    auto computeCoord = [angle](geo::View_t view, const TVector3& wire_center) -> float {
        if (view == geo::kW) return wire_center.Z();
        const float sign = (view == geo::kU) ? 1.0f : -1.0f;
        return (wire_center.Z() * std::cos(sign * angle)) - (wire_center.Y() * std::sin(sign * angle));
    };
    auto applyKernel = [](Image& img, size_t row, size_t col, size_t kernel, float value) {
        int half = static_cast<int>(kernel) / 2;
        for (int dr = -half; dr <= half; ++dr) {
            for (int dc = -half; dc <= half; ++dc) {
                int new_row = static_cast<int>(row) + dr;
                int new_col = static_cast<int>(col) + dc;
                if (new_row < 0 || new_row >= static_cast<int>(img.properties().height()))
                    continue;
                if (new_col < 0 || new_col >= static_cast<int>(img.properties().width()))
                    continue;
                img.set(new_row, new_col, value, /*accumulate=*/false);
            }
        }
    };
    for (const auto& hit : hits) {
        geo::View_t view = hit->View();
        double x = det_props->ConvertTicksToX(hit->PeakTime(), hit->WireID());
        TVector3 wire_center = geo.Cryostat(hit->WireID().Cryostat)
                                    .TPC(hit->WireID().TPC)
                                    .Plane(hit->WireID().Plane)
                                    .Wire(hit->WireID().Wire)
                                    .GetCenter();
        float transformed_coord = computeCoord(view, wire_center);
        for (auto& img : truth) {
            if (img.properties().view() != view) continue;
            size_t row = img.properties().row(x);
            size_t col = img.properties().col(transformed_coord);
            if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1))
                continue;
            size_t ker = img.properties().kernel();
            applyKernel(img, row, col, ker, static_cast<float>(signature::kNoiseSignature));
        }
    }
    for (const auto& hit : hits) {
        const auto& assoc_part = mcp_bkth_assoc.at(hit.key());
        const auto& assoc_data = mcp_bkth_assoc.data(hit.key());
        if (assoc_part.empty() || assoc_data.empty())
            continue;
        int track_id = -1;
        bool valid_hit = false;
        for (size_t i = 0; i < assoc_data.size(); ++i) {
            if (assoc_data[i]->isMaxIDE == 1) {
                track_id = assoc_part[i]->TrackId();
                valid_hit = true;
                break;
            }
        }
        if (!valid_hit) continue;
        geo::View_t hit_view = hit->View();
        signature::SignatureType label = signature::kEmptySignature;
        bool found = false;
        if (hit_view == geo::kU) {
            auto it = truth_particles_U.find(track_id);
            if (it == truth_particles_U.end()) continue;
            label = it->second;
            found = true;
        } else if (hit_view == geo::kV) {
            auto it = truth_particles_V.find(track_id);
            if (it == truth_particles_V.end()) continue;
            label = it->second;
            found = true;
        } else if (hit_view == geo::kW) {
            auto it = truth_particles_W.find(track_id);
            if (it == truth_particles_W.end()) continue;
            label = it->second;
            found = true;
        }
        if (!found) continue;
        double x = det_props->ConvertTicksToX(hit->PeakTime(), hit->WireID());
        TVector3 wire_center = geo.Cryostat(hit->WireID().Cryostat)
                                    .TPC(hit->WireID().TPC)
                                    .Plane(hit->WireID().Plane)
                                    .Wire(hit->WireID().Wire)
                                    .GetCenter();
        geo::View_t view = hit->View();
        float transformed_coord = computeCoord(view, wire_center);
        for (auto& img : truth) {
            if (img.properties().view() != view) continue;
            size_t row = img.properties().row(x);
            size_t col = img.properties().col(transformed_coord);
            if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1))
                continue;
            size_t ker = img.properties().kernel();
            applyKernel(img, row, col, ker, static_cast<float>(label));
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
        tree_->Branch("type", &type_);
        tree_->Branch("signatures", &signatures_);
        tree_->Branch("clarity", &clarity_);
        tree_->Branch("input", &input_);
        tree_->Branch("truth", &truth_);
    }
    void reset() {
        run_ = subrun_ = event_ = 0;
        type_ = static_cast<int>(signature::EventType::kOther);
        signatures_.clear();
        clarity_.clear();
        input_.clear();
        truth_.clear();
    }
    void add(const art::Event& e, 
             const std::vector<art::Ptr<recob::Wire>>& wires,
             const std::vector<art::Ptr<recob::Hit>>& hits,
             const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
             const signature::EventClassifier& classifier,
             const std::vector<ImageProperties>& properties) {
        this->reset();
        run_ = e.run();
        subrun_ = e.subRun();
        event_ = e.event();
        signature::EventType event_type = classifier.classifyEvent(e);
        type_ = static_cast<int>(event_type);
        std::vector<Image> input_images = constructInputImages(properties, wires, geo_);
        signature::Pattern pattern = classifier.getPattern(e);
        std::vector<Image> truth_images = constructTruthImages(properties, hits, mcp_bkth_assoc, pattern, classifier, geo_);
        for (const auto& [type, sig] : pattern) {
            signatures_.push_back(static_cast<int>(type));
            clarity_.push_back({
                classifier.passClarity(sig, common::TPC_VIEW_U) ? 1 : 0,
                classifier.passClarity(sig, common::TPC_VIEW_V) ? 1 : 0,
                classifier.passClarity(sig, common::TPC_VIEW_W) ? 1 : 0
            });
        }
        std::vector<geo::View_t> ordered_views = { geo::kU, geo::kV, geo::kW };
        for (auto view : ordered_views) {
            auto input_iter = std::find_if(input_images.begin(), input_images.end(),
                                        [view](const Image& img) {
                                            return img.properties().view() == view;
                                        });
            auto truth_iter = std::find_if(truth_images.begin(), truth_images.end(),
                                        [view](const Image& img) {
                                            return img.properties().view() == view;
                                        });                       
            if (input_iter != input_images.end() && truth_iter != truth_images.end()) {
                if (input_iter->properties().width() != truth_iter->properties().width() ||
                    input_iter->properties().height() != truth_iter->properties().height())
                    continue;
                input_.push_back(input_iter->data());
                truth_.push_back(truth_iter->data());
            }
        }
        tree_->Fill();
    }
private:
    const geo::GeometryCore& geo_;
    TTree* tree_;
    int run_, subrun_, event_;
    int type_;
    std::vector<int> signatures_;
    std::vector<std::vector<int>> clarity_;
    std::vector<int> planes_;
    std::vector<int> width_; 
    std::vector<int> height_;
    std::vector<std::vector<float>> input_;
    std::vector<std::vector<float>> truth_; 
};

} 

#endif // IMAGEPROCESSOR_H