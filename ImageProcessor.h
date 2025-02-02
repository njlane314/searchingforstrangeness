// ImageProcessor.h
#ifndef IMAGEPROCESSOR_H
#define IMAGEPROCESSOR_H

#include <vector>
#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <TFile.h>
#include <TTree.h>
#include <TObject.h>
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace image {

enum SignatureType {
    SignatureEmpty = 0,
    SignatureNoise,
    SignaturePrimaryMuon,
    SignatureChargedKaon,
    SignatureKaonShort,
    SignatureLambda,
    SignatureChargedSigma
};

class ImageMeta {
public:
    ImageMeta(double center_x, double center_y, size_t height, size_t width, double pixel_w, double pixel_h, geo::View_t view)
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

class Image : public TObject {
public:
    Image() = default;
    Image(const ImageMeta& meta)
        : meta_(meta), pixels_(meta.height() * meta.width(), 0.0f) {}

    void set(size_t row, size_t col, float value, bool accumulate = true) {
        size_t idx = meta_.index(row, col);
        if (idx != static_cast<size_t>(-1)) {
            if (accumulate)
                pixels_[idx] += value;
            else
                pixels_[idx] = value;
        }
    }

    void copy_col(size_t col, const std::vector<float>& data) {
        if (col >= meta_.width()) return;
        for (size_t row = 0; row < meta_.height(); ++row) {
            size_t idx = meta_.index(row, col);
            if (idx != static_cast<size_t>(-1)) {
                pixels_[idx] = data[row];
            }
        }
    }

    float get(size_t row, size_t col) const {
        size_t idx = meta_.index(row, col);
        return (idx != static_cast<size_t>(-1)) ? pixels_[idx] : 0.0f;
    }

    void clear(float value = 0.0f) {
        std::fill(pixels_.begin(), pixels_.end(), value);
    }

    const std::vector<float>& data() const {
        return pixels_;
    }

private:
    ImageMeta meta_;
    std::vector<float> pixels_;
};

std::vector<Image> SimChannelToImages(const std::vector<ImageMeta>& metas, const std::vector<sim::SimChannel>& channels, const geo::GeometryCore& geo) {
    std::unordered_map<int, SignatureType> track_signatures;
    std::unordered_map<size_t, float> pixel_energy;

    std::vector<Image> images;
    for (const auto& meta : metas) images.emplace_back(meta);

    for (const auto& channel : channels) {
        auto ch_id = channel.Channel();
        auto wire_ids = geo.ChannelToWire(ch_id);
        if (wire_ids.empty()) continue;

        auto wire_id = wire_ids.front();
        auto view = geo.View(wire_id);

        for (const auto& tdc_ide_pair : channel.TDCIDEMap()) {
            int tdc = tdc_ide_pair.first;
            for (const auto& ide : tdc_ide_pair.second) {
                double x = geo.ConvertTicksToX(tdc, wire_id.Plane, wire_id.TPC, wire_id.Cryostat);
                TVector3 wire_center = geo.Cryostat(wire_id.Cryostat).TPC(wire_id.TPC).Plane(wire_id.Plane).Wire(wire_id.Wire).GetCenter();

                float coord = (view == geo::kW) ? wire_center.Z() :
                              (wire_center.Z() * std::cos((view == geo::kU ? 1 : -1) * 1.04719758034)) - 
                              (wire_center.Y() * std::sin((view == geo::kU ? 1 : -1) * 1.04719758034));

                SignatureType sig = SignatureEmpty;
                if (ide.trackID > 0) {
                    if (track_signatures.find(ide.trackID) == track_signatures.end()) {
                        track_signatures[ide.trackID] = SignatureNoise;
                    }
                    sig = track_signatures[ide.trackID];
                }

                for (auto& img : images) {
                    if (img.meta_.view() == view) {
                        size_t col = img.meta_.col(coord);
                        size_t row = img.meta_.row(x);
                        size_t idx = img.meta_.index(row, col);
                        if (idx != static_cast<size_t>(-1)) {
                            if (pixel_energy[idx] < ide.energy) {
                                pixel_energy[idx] = ide.energy;
                                img.set(row, col, static_cast<float>(sig), false);
                            }
                        }
                    }
                }
            }
        }
    }
    return images;
}

std::vector<Image> WiresToImages(const std::vector<ImageMeta>& metas, const std::vector<recob::Wire>& wires, const geo::GeometryCore& geo) {
    std::vector<Image> images;
    for (const auto& meta : metas) images.emplace_back(meta);

    for (const auto& wire : wires) {
        auto ch_id = wire.Channel();
        auto wire_ids = geo.ChannelToWire(ch_id);
        if (wire_ids.empty()) continue;

        auto view = geo.View(wire_ids.front());
        for (const auto& range : wire.SignalROI().get_ranges()) {
            const auto& adcs = range.data();
            int start_idx = range.begin_index();

            for (size_t idx = 0; idx < adcs.size(); ++idx) {
                double x = geo.ConvertTicksToX(start_idx + idx, wire_ids.front().Plane, wire_ids.front().TPC, wire_ids.front().Cryostat);
                TVector3 wire_center = geo.Cryostat(wire_ids.front().Cryostat).TPC(wire_ids.front().TPC).Plane(wire_ids.front().Plane).Wire(wire_ids.front().Wire).GetCenter();

                float coord = (view == geo::kW) ? wire_center.Z() :
                              (wire_center.Z() * std::cos((view == geo::kU ? 1 : -1) * 1.04719758034)) - 
                              (wire_center.Y() * std::sin((view == geo::kU ? 1 : -1) * 1.04719758034));

                for (auto& img : images) {
                    if (img.meta_.view() == view) {
                        size_t row = img.meta_.row(x);
                        size_t col = img.meta_.col(coord);
                        img.set(row, col, adcs[idx]);
                    }
                }
            }
        }
    }
    return images;
}

class ImageManager {
public:
    ImageManager(TTree* tree, const geo::GeometryCore& geo)
        : tree_(tree), geo_(geo) {
        tree_->Branch("run", &run_);
        tree_->Branch("subrun", &subrun_);
        tree_->Branch("event", &event_);
        tree_->Branch("signal", &signal_);
        tree_->Branch("wires", &wires_);
        tree_->Branch("sim_channels", &sim_channels_);
    }

    void reset() {
        run_ = subrun_ = event_ = 0;
        signal_ = false;
        wires_.clear();
        sim_channels_.clear();
    }

    void add(int run, int subrun, int event, bool signal,
             const std::vector<recob::Wire>& wires,
             const std::vector<sim::SimChannel>& channels,
             const std::vector<ImageMeta>& metas) {
        run_ = run;
        subrun_ = subrun;
        event_ = event;
        signal_ = signal;

        wires_ = WiresToImages(metas, wires, geo_);
        sim_channels_ = SimChannelToImages(metas, channels, geo_);

        tree_->Fill();
    }

private:
    TTree* tree_;
    const geo::GeometryCore& geo_;
    int run_, subrun_, event_;
    bool signal_;
    std::vector<Image> wires_;
    std::vector<Image> sim_channels_;
};

} // namespace image

#endif // IMAGEPROCESSOR_H
