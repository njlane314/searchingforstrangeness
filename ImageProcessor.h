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
#include "EventClassifier.h"
#include "LabelClassifier.h"
#include "art/Framework/Principal/Event.h"
#include "larcore/Geometry/WireReadout.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"

namespace image 
{
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
        geo::View_t view() const { return prop_.view(); }
        size_t height() const { return prop_.height(); }
        size_t width() const { return prop_.width(); }
    private:
        ImageProperties prop_;
        std::vector<float> pixels_;
    };

    void constructAllWireImages(const art::Event& e,
                                const std::vector<ImageProperties>& properties,
                                const std::vector<art::Ptr<recob::Wire>>& wires,
                                const art::FindManyP<recob::Hit>& wire_hit_assoc,
                                const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
                                const signature::Pattern& pattern,
                                const std::vector<signature::Label>& particle_labels,
                                const geo::GeometryCore* geom,
                                const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
                                const std::vector<bool>& bad_channel_mask,
                                std::vector<Image>& input_images,
                                std::vector<Image>& truth_images,
                                std::vector<Image>& label_images,
                                float adc_threshold = 10.0f) {
        input_images.clear();
        truth_images.clear();
        label_images.clear();

        for (const auto& prop : properties) {
            input_images.emplace_back(prop); 
            Image truth_img(prop);
            truth_img.clear(signature::kEmptySignature); 
            truth_images.push_back(std::move(truth_img));
            Image label_img(prop);
            label_img.clear(static_cast<float>(signature::Label::undefined)); 
            label_images.push_back(std::move(label_img));
        }

        std::array<std::unordered_map<int, signature::SignatureType>, 3> truth_maps;
        for (const auto& particle : mcp_vec) {
            for (auto& map : truth_maps) {
                map[particle->TrackId()] = signature::kNeutrinoSignature; 
            }
        }
        for (const auto& [type, sig] : pattern) {
            for (const auto& particle : sig) {
                for (auto& map : truth_maps) {
                    map[particle->TrackId()] = type; 
                }
            }
        }

        std::unordered_map<int, signature::Label> trackid_to_label;
        for (size_t i = 0; i < mcp_vec.size() && i < particle_labels.size(); ++i) {
            trackid_to_label[mcp_vec[i]->TrackId()] = particle_labels[i];
        }

        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

        for (size_t wire_idx = 0; wire_idx < wires.size(); ++wire_idx) {
            const auto& wire = wires[wire_idx];
            auto ch_id = wire->Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue; 

            std::vector<geo::WireID> wire_ids = wireReadout.ChannelToWire(ch_id);
            if (wire_ids.empty()) continue; 
            geo::View_t view = wireReadout.Plane(wire_ids.front().planeID()).View();
            size_t view_idx = view - geo::kU; 

            geo::Point_t center = wireReadout.Wire(wire_ids.front()).GetCenter();
            TVector3 wire_center(center.X(), center.Y(), center.Z());
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                    (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            const auto& hits = wire_hit_assoc.at(wire_idx);

            for (const auto& range : wire->SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t idx = 0; idx < adcs.size(); ++idx) {
                    int tick = start_tick + idx;
                    double x = det_props.ConvertTicksToX(tick, wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(x);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue; 

                    signature::SignatureType truth_label = signature::kEmptySignature;
                    signature::Label particle_label = signature::Label::undefined;

                    for (const auto& hit : hits) {
                        if (tick >= hit->StartTick() && tick < hit->EndTick()) {
                            auto data = mcp_bkth_assoc.data(hit.key());
                            if (data.empty()) {
                                truth_label = signature::kCosmicSignature;
                                particle_label = signature::Label::cosmic;
                            }
                            else {
                                for (size_t i = 0; i < data.size(); ++i) {
                                    if (data[i]->isMaxIDE == 1) {
                                        int track_id = mcp_bkth_assoc.at(hit.key())[i]->TrackId();
                                        if (auto it = truth_maps[view_idx].find(track_id); it != truth_maps[view_idx].end()) {
                                            truth_label = it->second;
                                        }
                                        if (auto it = trackid_to_label.find(track_id); it != trackid_to_label.end()) {
                                            particle_label = it->second;
                                        }
                                        break;
                                    }
                                }
                            }
                            break;
                        }
                    }

                    if (adcs[idx] > adc_threshold) {
                        input_images[view_idx].set(row, col, adcs[idx]);
                        truth_images[view_idx].set(row, col, static_cast<float>(truth_label), false);
                        label_images[view_idx].set(row, col, static_cast<float>(particle_label), false);
                    }
                }
            }
        }
    }

    /*std::vector<Image> constructInputWireImages(const art::Event& e,
                                                const std::vector<ImageProperties>& properties,
                                                const std::vector<art::Ptr<recob::Wire>>& wires,
                                                const geo::GeometryCore* geom,
                                                const std::vector<bool>& bad_channel_mask) {
        std::vector<Image> input;
        for (const auto& prop : properties) input.emplace_back(prop);
        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();
        for (const auto& wire : wires) {
            auto ch_id = wire->Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue;
            std::vector<geo::WireID> wire_ids = wireReadout.ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = wireReadout.Plane(wire_ids.front().planeID()).View();
            TVector3 wire_center = wireReadout.Wire(wire_ids.front()).GetCenter();
            float coord = (view == geo::kW) ? wire_center.Z() :
                        (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                            (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));
            for (const auto& range : wire->SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_idx = range.begin_index();
                for (size_t idx = 0; idx < adcs.size(); ++idx) {
                    int tdc = start_idx + idx;
                    double x = det_props.ConvertTicksToX(tdc, wire_ids.front().planeID());
                    for (auto& img : input) {
                        if (img.view() == view) {
                            size_t row = properties[view - geo::kU].row(x);
                            size_t col = properties[view - geo::kU].col(coord);
                            if (row != static_cast<size_t>(-1) && col != static_cast<size_t>(-1)) {
                                img.set(row, col, adcs[idx]);
                            }
                        }
                    }
                }
            }
        }
        return input;
    }

    std::vector<Image> constructTruthWireImages(const art::Event& e,
                                                const std::vector<ImageProperties>& properties,
                                                const std::vector<art::Ptr<recob::Wire>>& wires,
                                                const art::FindManyP<recob::Hit>& wire_hit_assoc,
                                                const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
                                                const signature::Pattern& pattern,
                                                const signature::EventClassifier& classifier,
                                                const geo::GeometryCore* geom,
                                                const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
                                                const std::vector<bool>& bad_channel_mask) {
        std::vector<Image> truth;
        truth.reserve(properties.size());
        for (const auto& prop : properties) {
            Image img(prop);
            img.clear(signature::kEmptySignature);
            truth.push_back(std::move(img));
        }

        std::array<std::unordered_map<int, signature::SignatureType>, 3> truth_maps;
        for (const auto& particle : mcp_vec) {
            for (auto& map : truth_maps) map[particle->TrackId()] = signature::kNeutrinoSignature;
        }
        for (const auto& [type, sig] : pattern) {
            for (const auto& particle : sig) {
                for (auto& map : truth_maps) map[particle->TrackId()] = type;
            }
        }

        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

        for (size_t wire_idx = 0; wire_idx < wires.size(); ++wire_idx) {
            const auto& wire = wires[wire_idx];
            auto ch_id = wire->Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue;

            const auto& hits = wire_hit_assoc.at(wire_idx);
            if (hits.empty()) continue;

            std::vector<geo::WireID> wire_ids = wireReadout.ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = wireReadout.Plane(wire_ids.front().planeID()).View();
            size_t view_idx = view - geo::kU;

            signature::SignatureType label = signature::kCosmicSignature;
            for (const auto& hit : hits) {
                if (auto data = mcp_bkth_assoc.data(hit.key()); !data.empty()) {
                    auto& map = truth_maps[view_idx];
                    for (size_t i = 0; i < data.size(); ++i) {
                        if (data[i]->isMaxIDE == 1) {
                            int track_id = mcp_bkth_assoc.at(hit.key())[i]->TrackId();
                            if (map.find(track_id) != map.end()) {
                                label = map[track_id];
                                break;
                            }
                        }
                    }
                }
                if (label != signature::kCosmicSignature) break;
            }

            TVector3 wire_center = wireReadout.Wire(wire_ids.front()).GetCenter();
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                    (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            for (const auto& range : wire->SignalROI().get_ranges()) {
                int start_tick = range.begin_index();
                double start_x = det_props.ConvertTicksToX(start_tick, wire_ids.front().planeID());
                double end_x = det_props.ConvertTicksToX(start_tick + range.data().size() - 1, wire_ids.front().planeID());

                for (auto& img : truth) {
                    if (img.view() != view) continue;
                    auto min_x = img.origin_x();
                    auto max_x = img.max_x();
                    if (start_x >= max_x || end_x <= min_x) continue;

                    size_t r_start = properties[view_idx].row(std::max(start_x, min_x));
                    size_t r_end = properties[view_idx].row(std::min(end_x, max_x));
                    size_t col = properties[view_idx].col(wire_coord);

                    if (r_start == static_cast<size_t>(-1) || r_end == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    for (size_t r = r_start; r <= r_end && r < img.height() && col < img.width(); ++r) {
                        img.set(r, col, static_cast<float>(label), false);
                    }
                }
            }
        }

        return truth;
    }

    std::vector<Image> constructLabelWireImages(const art::Event& e,
                                                const std::vector<ImageProperties>& properties,
                                                const std::vector<art::Ptr<recob::Wire>>& wires, 
                                                const art::FindManyP<recob::Hit>& wire_hit_assoc,
                                                const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
                                                const std::vector<signature::Label>& particle_labels, 
                                                const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
                                                const geo::GeometryCore* geom, 
                                                const std::vector<bool>& bad_channel_mask) {
        std::vector<Image> label_images; 
        for (const auto& prop : properties) {
            Image img(prop); 
            img.clear(static_cast<float>(signature::Label::undefined));
            label_images.push_back(std::move(img));
        }

        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

        std::unordered_map<int, signature::Label> trackid_to_label;
        for (size_t i = 0; i < mcp_vec.size() && i < particle_labels.size(); ++i) {
            trackid_to_label[mcp_vec[i]->TrackId()] = particle_labels[i];  
        }

        for (size_t wire_idx = 0; wire_idx < wires.size(); ++wire_idx) {
            const auto& wire = wires[wire_idx];
            auto ch_id = wire->Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue; 

            const auto& hits = wire_hit_assoc.at(wire_idx);
            if (hits.empty()) continue;  

            std::vector<geo::WireID> wire_ids = wireReadout.ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = wireReadout.Plane(wire_ids.front().planeID()).View();

            signature::Label label = signature::Label::undefined;
            for (const auto& hit : hits) {
                try {
                    const auto& assoc_part = mcp_bkth_assoc.at(hit.key());
                    const auto& assoc_data = mcp_bkth_assoc.data(hit.key());
                    if (assoc_part.empty() || assoc_data.empty()) continue;

                    for (size_t i = 0; i < assoc_data.size(); ++i) {
                        if (assoc_data[i]->isMaxIDE == 1) { 
                            int track_id = assoc_part[i]->TrackId();
                            auto it = trackid_to_label.find(track_id);
                            if (it != trackid_to_label.end()) {
                                label = it->second;
                                break;
                            }
                        }
                    }
                } catch (const std::out_of_range&) {
                    continue; 
                }
                if (label != signature::Label::undefined) break;  
            }

            if (label == signature::Label::undefined) continue;  

            TVector3 wire_center = wireReadout.Wire(wire_ids.front()).GetCenter();
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                    (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            for (const auto& range : wire->SignalROI().get_ranges()) {
                int start_tick = range.begin_index();
                double start_x = det_props.ConvertTicksToX(start_tick, wire_ids.front().planeID());
                double end_x = det_props.ConvertTicksToX(start_tick + range.data().size() - 1, wire_ids.front().planeID());

                for (auto& img : label_images) {
                    if (img.view() != view) continue;  
                    auto min_x = img.origin_x();
                    auto max_x = img.max_x();
                    if (start_x >= max_x || end_x <= min_x) continue;  

                    size_t view_idx = view - geo::kU;
                    size_t r_start = properties[view_idx].row(std::max(start_x, min_x));
                    size_t r_end = properties[view_idx].row(std::min(end_x, max_x));
                    size_t col = properties[view_idx].col(wire_coord);

                    if (r_start == static_cast<size_t>(-1) || r_end == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    for (size_t r = r_start; r <= r_end && r < img.height() && col < img.width(); ++r) {
                        img.set(r, col, static_cast<float>(label), false);
                    }
                }
            }
        }

        return label_images;
    }*/

    /*std::vector<Image> constructTruthHitImages(const art::Event& e,
                                               const std::vector<ImageProperties>& properties,
                                               const std::vector<art::Ptr<recob::Hit>>& hits,
                                               const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
                                               const signature::Pattern& pattern,
                                               const signature::EventClassifier& classifier,
                                               const geo::GeometryCore* geom,
                                               const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
                                               const std::vector<bool>& bad_channel_mask) {
        std::vector<Image> truth;
        truth.reserve(properties.size());
        for (const auto& prop : properties) {
            Image img(prop);
            img.clear(signature::kEmptySignature);
            truth.push_back(std::move(img));
        }
        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        std::unordered_map<int, signature::SignatureType> truth_particles_U;
        std::unordered_map<int, signature::SignatureType> truth_particles_V;
        std::unordered_map<int, signature::SignatureType> truth_particles_W;
        for (const auto& particle : mcp_vec) {
            truth_particles_U[particle->TrackId()] = signature::kNeutrinoSignature;
            truth_particles_V[particle->TrackId()] = signature::kNeutrinoSignature;
            truth_particles_W[particle->TrackId()] = signature::kNeutrinoSignature; 
        }
        for (const auto& [type, sig] : pattern) {
            for (const auto& particle : sig) {
                truth_particles_U[particle->TrackId()] = type;
                truth_particles_V[particle->TrackId()] = type;
                truth_particles_W[particle->TrackId()] = type;
            }
        }
        constexpr float angle = 1.04719758034f;
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
                    if (new_row < 0 || new_row >= static_cast<int>(img.height()))
                        continue;
                    if (new_col < 0 || new_col >= static_cast<int>(img.width()))
                        continue;
                    if (img.get(new_row, new_col) != 0.0f) 
                        img.set(new_row, new_col, value, false);
                }
            }
        };
        for (const auto& hit : hits) {
            auto ch_id = hit->Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue;
            geo::View_t view = hit->View();
            double x = det_props.ConvertTicksToX(hit->PeakTime(), hit->WireID().planeID());
            TVector3 wire_center = geom->Wire(hit->WireID()).GetCenter();
            float transformed_coord = computeCoord(view, wire_center);
            for (auto& img : truth) {
                if (img.view() != view) continue;
                size_t row = img.row(x);
                size_t col = img.col(transformed_coord);
                if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1))
                    continue;
                size_t ker = img.kernel();
                applyKernel(img, row, col, ker, static_cast<float>(signature::kCosmicSignature));
            }
        }
        for (const auto& hit : hits) {
            auto ch_id = hit->Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue;
            try {
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
                if (!valid_hit) 
                    continue;
                geo::View_t hit_view = hit->View();
                signature::SignatureType label = signature::kEmptySignature;
                bool found = false;
                if (hit_view == geo::kU) {
                    auto it = truth_particles_U.find(track_id);
                    if (it == truth_particles_U.end()) 
                        continue;
                    label = it->second;
                    found = true;
                } else if (hit_view == geo::kV) {
                    auto it = truth_particles_V.find(track_id);
                    if (it == truth_particles_V.end()) 
                        continue;
                    label = it->second;
                    found = true;
                } else if (hit_view == geo::kW) {
                    auto it = truth_particles_W.find(track_id);
                    if (it == truth_particles_W.end()) 
                        continue;
                    label = it->second;
                    found = true;
                }
                if (!found) 
                    continue;
                double x = det_props.ConvertTicksToX(hit->PeakTime(), hit->WireID().planeID());
                TVector3 wire_center = geom->Wire(hit->WireID()).GetCenter();
                geo::View_t view = hit->View();
                float transformed_coord = computeCoord(view, wire_center);
                for (auto& img : truth) {
                    if (img.view() != view) 
                        continue;
                    size_t row = img.row(x);
                    size_t col = img.col(transformed_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1))
                        continue;
                    size_t ker = img.kernel();
                    applyKernel(img, row, col, ker, static_cast<float>(label));
                }
            } catch (const std::out_of_range&) {
                continue; 
            }
        }
        return truth;
    }

    std::vector<Image> constructInputHitImages(const art::Event& e,
                                               const std::vector<ImageProperties>& properties,
                                               const std::vector<art::Ptr<recob::Hit>>& hits,
                                               const geo::GeometryCore* geom) {
        std::vector<Image> input;
        for (const auto& prop : properties) {
            input.emplace_back(prop);
        }
        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        for (const auto& hit : hits) {
            geo::View_t view = hit->View();
            double x = det_props.ConvertTicksToX(hit->PeakTime(), hit->WireID().planeID());
            TVector3 wire_center = geom->Wire(hit->WireID()).GetCenter();
            float coord = (view == geo::kW) ? wire_center.Z() :
                          (view == geo::kU) ? (wire_center.Z() * cos(1.04719758034) - wire_center.Y() * sin(1.04719758034)) :
                                              (wire_center.Z() * cos(-1.04719758034) - wire_center.Y() * sin(-1.04719758034));
            for (auto& img : input) {
                if (img.view() == view) {
                    size_t row = img.row(x);
                    size_t col = img.col(coord);
                    if (row != static_cast<size_t>(-1) && col != static_cast<size_t>(-1)) {
                        img.set(row, col, hit->Integral(), true); 
                    }
                }
            }
        }
        return input;
    }*/

    std::vector<std::vector<float>> extractImages(const std::vector<image::Image>& images) {
        std::vector<std::vector<float>> data;
        data.reserve(images.size());
        for (const auto& img : images) {
            data.push_back(img.data());
        }
        return data;
    }
} 

#endif