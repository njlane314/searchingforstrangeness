#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

#include <optional>
#include <array>
#include "Imaging/Image.h"
#include "Imaging/SemanticClassifier.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfo/DetectorPropertiesData.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <TVector3.h>
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace image {

struct PixelImageOptions {
    bool is_data{false};
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
    } producers;
    float adc_threshold{4.0f};
    const std::set<unsigned int>* bad_channels{nullptr};
    image::SemanticClassifier* semantic{nullptr};
};

struct CalibrationContext {
    calo::CalorimetryAlg* calo{nullptr};
    detinfo::DetectorClocksData const* clock{nullptr};
    detinfo::DetectorPropertiesData const* detprop{nullptr};
    double T0_ticks{0.0};
    bool enabled() const { return calo && clock && detprop; }
};

class PixelImageBuilder {
public:
    PixelImageBuilder(geo::GeometryCore const& geo,
                      PixelImageOptions const& opts)
        : geo_{&geo}, opts_{opts} {}

    void build(const art::Event &event,
               const std::vector<art::Ptr<recob::Hit>> &hits,
               const std::vector<ImageProperties> &properties,
               std::vector<Image<float>> &detector_images,
               std::vector<Image<int>> &semantic_images,
               const detinfo::DetectorProperties *detp,
               std::optional<CalibrationContext> const& cal = std::nullopt) const;

private:
    struct BuildContext {
        const std::vector<ImageProperties> &properties;
        std::vector<Image<float>> &detector_images;
        std::vector<Image<int>> &semantic_images;
        const art::FindManyP<recob::Hit> &wire_hit_assoc;
        const std::set<art::Ptr<recob::Hit>> &hit_set;
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc;
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector;
        const std::map<int, size_t> &trackid_to_index;
        const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector;
        bool is_data;
        bool has_mcps;
        const geo::GeometryCore *geo;
        const detinfo::DetectorProperties *detp;
        float adc_image_threshold;
        const std::set<unsigned int> *bad_channels;
        calo::CalorimetryAlg* calo_alg;
        detinfo::DetectorClocksData const* clock_data;
        detinfo::DetectorPropertiesData const* detprop_data;
        double T0_ticks{0.0};
    };

    static SemanticClassifier::SemanticLabel labelSemanticPixels(
        const art::Ptr<recob::Hit> &matched_hit,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
        bool has_mcps);

    static void fillDetectorImage(
        const recob::Wire &wire,
        size_t wire_idx,
        const std::set<art::Ptr<recob::Hit>> &hit_set,
        BuildContext const& ctx);

    geo::GeometryCore const* geo_{nullptr};
    PixelImageOptions opts_;
};

}

inline SemanticClassifier::SemanticLabel image::PixelImageBuilder::labelSemanticPixels(
    const art::Ptr<recob::Hit> &matched_hit,
    const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
        &mcp_bkth_assoc,
    const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
    const std::map<int, size_t> &trackid_to_index,
    const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
    bool has_mcps) {
    SemanticClassifier::SemanticLabel semantic_pixel_label =
        SemanticClassifier::SemanticLabel::Cosmic;

    if (!has_mcps || !mcp_vector.isValid()) {
        return semantic_pixel_label;
    }

    std::vector<art::Ptr<simb::MCParticle>> mcp_particles_ass_to_hit;
    std::vector<anab::BackTrackerHitMatchingData const *> bkth_data_ass_to_hit;
    mcp_bkth_assoc.get(matched_hit.key(), mcp_particles_ass_to_hit,
                       bkth_data_ass_to_hit);
    if (!bkth_data_ass_to_hit.empty()) {
        float max_ide_fraction = -1.0f;
        int best_match_track_id = -1;

        for (size_t i_bkth = 0; i_bkth < bkth_data_ass_to_hit.size(); ++i_bkth) {
            if (bkth_data_ass_to_hit[i_bkth] &&
                bkth_data_ass_to_hit[i_bkth]->ideFraction > max_ide_fraction) {
                max_ide_fraction = bkth_data_ass_to_hit[i_bkth]->ideFraction;
                best_match_track_id = mcp_particles_ass_to_hit[i_bkth]->TrackId();
            }
        }

        if (best_match_track_id != -1) {
            auto it_trackid = trackid_to_index.find(best_match_track_id);
            if (it_trackid != trackid_to_index.end()) {
                size_t particle_mcp_idx = it_trackid->second;

                if (particle_mcp_idx < semantic_label_vector.size()) {
                    if (max_ide_fraction <= 0.5f) {
                        semantic_pixel_label =
                            SemanticClassifier::SemanticLabel::Ambiguous;
                    } else {
                        semantic_pixel_label =
                            semantic_label_vector[particle_mcp_idx];
                    }
                }
            }
        }
    }

    return semantic_pixel_label;
}

inline void image::PixelImageBuilder::fillDetectorImage(
    const recob::Wire &wire,
    size_t wire_idx,
    const std::set<art::Ptr<recob::Hit>> &hit_set,
    BuildContext const& ctx) {
    auto ch_id = wire.Channel();
    if (ctx.bad_channels && ctx.bad_channels->count(ch_id)) {
        return;
    }

    std::vector<geo::WireID> wire_ids = ctx.geo->ChannelToWire(ch_id);
    if (wire_ids.empty()) {
        return;
    }

    geo::View_t view = ctx.geo->View(wire_ids.front().planeID());
    size_t view_idx = static_cast<size_t>(view);
    const geo::WireGeo *wire_geo = ctx.geo->WirePtr(wire_ids.front());
    TVector3 center = wire_geo->GetCenter();
    TVector3 wire_center(center.X(), center.Y(), center.Z());

    constexpr double plus60  =  1.04719758034;
    constexpr double minus60 = -1.04719758034;
    double wire_coord =
        (view == geo::kW)
            ? wire_center.Z()
        : (view == geo::kU)
            ? (wire_center.Z() * std::cos(plus60)  -
               wire_center.Y() * std::sin(plus60))
            : (wire_center.Z() * std::cos(minus60) -
               wire_center.Y() * std::sin(minus60));

    auto hits_for_wire = ctx.wire_hit_assoc.at(wire_idx);
    std::vector<art::Ptr<recob::Hit>> filtered_hits;
    filtered_hits.reserve(hits_for_wire.size());
    for (const auto &hit : hits_for_wire) {
        if (hit_set.count(hit)) {
            filtered_hits.push_back(hit);
        }
    }

    std::sort(filtered_hits.begin(), filtered_hits.end(),
              [](const art::Ptr<recob::Hit> &a, const art::Ptr<recob::Hit> &b) {
                  return a->StartTick() < b->StartTick();
              });

    size_t hit_index = 0;
    for (const auto &range : wire.SignalROI().get_ranges()) {
        const auto &adcs = range.data();
        int start_tick = range.begin_index();

        for (size_t adc_index = 0; adc_index < adcs.size(); ++adc_index) {
            int tick = start_tick + static_cast<int>(adc_index);
            double x_drift_pos = ctx.detp->ConvertTicksToX(static_cast<double>(tick),
                                                           wire_ids.front().planeID());

            auto row = ctx.properties[view_idx].row(x_drift_pos);
            auto col = ctx.properties[view_idx].col(wire_coord);
            if (!row || !col) {
                continue;
            }

            while (hit_index < filtered_hits.size() &&
                   filtered_hits[hit_index]->EndTick() <= tick) {
                ++hit_index;
            }

            if (hit_index < filtered_hits.size() &&
                filtered_hits[hit_index]->StartTick() <= tick &&
                tick < filtered_hits[hit_index]->EndTick()) {
                if (adcs[adc_index] > ctx.adc_image_threshold) {
                    float out_val = adcs[adc_index];
                    if (ctx.calo_alg && ctx.clock_data && ctx.detprop_data) {
                        const double pitch_cm_drift = ctx.properties[view_idx].pixel_w();
                        if (pitch_cm_drift > 0.) {
                            const double dQdx_ADC = static_cast<double>(adcs[adc_index]) / pitch_cm_drift;
                            const unsigned plane = wire_ids.front().planeID().Plane;
                            const double time_tick = static_cast<double>(tick);
                            const double dEdx_MeV_cm =
                                ctx.calo_alg->dEdx_AMP(*ctx.clock_data, *ctx.detprop_data,
                                                       dQdx_ADC, time_tick, plane,
                                                       ctx.T0_ticks);
                            out_val = static_cast<float>(dEdx_MeV_cm);
                        }
                    }
                    ctx.detector_images[view_idx].set(*row, *col, out_val);

                    if (!ctx.is_data) {
                        auto semantic_pixel_label = labelSemanticPixels(
                            filtered_hits[hit_index], ctx.mcp_bkth_assoc, ctx.mcp_vector,
                            ctx.trackid_to_index, ctx.semantic_label_vector, ctx.has_mcps);

                        ctx.semantic_images[view_idx].set(
                            *row, *col,
                            static_cast<int>(semantic_pixel_label), false);
                    }
                }
            }
        }
    }
}

inline void image::PixelImageBuilder::build(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<ImageProperties> &properties,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images,
    const detinfo::DetectorProperties *detp,
    std::optional<CalibrationContext> const& cal) const {
    detector_images.clear();
    semantic_images.clear();
    for (const auto &prop : properties) {
        Image<float> detector_image(prop);
        detector_image.clear(0.0f);
        detector_images.push_back(std::move(detector_image));
        Image<int> semantic_image(prop);
        semantic_image.clear(
            static_cast<int>(SemanticClassifier::SemanticLabel::Empty));
        semantic_images.push_back(std::move(semantic_image));
    }

    auto wire_vector =
        event.getValidHandle<std::vector<recob::Wire>>(opts_.producers.wire);
    auto hit_vector =
        event.getValidHandle<std::vector<recob::Hit>>(opts_.producers.hit);
    art::Handle<std::vector<simb::MCParticle>> mcp_vector;
    bool has_mcps = event.getByLabel(opts_.producers.mcp, mcp_vector);
    art::FindManyP<recob::Hit> wire_hit_assoc(wire_vector, event, opts_.producers.hit);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
        mcp_bkth_assoc(hit_vector, event, opts_.producers.bkt);
    std::vector<SemanticClassifier::SemanticLabel> semantic_label_vector;
    if (!opts_.is_data && has_mcps && mcp_vector.isValid() && opts_.semantic) {
        semantic_label_vector = opts_.semantic->classifyParticles(event);
    }
    std::map<int, size_t> trackid_to_index;
    if (!opts_.is_data && has_mcps && mcp_vector.isValid()) {
        for (size_t i = 0; i < mcp_vector->size(); ++i) {
            trackid_to_index[mcp_vector->at(i).TrackId()] = i;
        }
    }

    std::set<art::Ptr<recob::Hit>> hit_set(hits.begin(), hits.end());
    calo::CalorimetryAlg* calo_alg = (cal && cal->enabled()) ? cal->calo : nullptr;
    detinfo::DetectorClocksData const* clock_data = (cal && cal->enabled()) ? cal->clock : nullptr;
    detinfo::DetectorPropertiesData const* detprop_data = (cal && cal->enabled()) ? cal->detprop : nullptr;
    double T0_ticks = (cal && cal->enabled()) ? cal->T0_ticks : 0.0;

    BuildContext ctx{
        properties,
        detector_images,
        semantic_images,
        wire_hit_assoc,
        hit_set,
        mcp_bkth_assoc,
        mcp_vector,
        trackid_to_index,
        semantic_label_vector,
        opts_.is_data,
        has_mcps,
        geo_,
        detp,
        opts_.adc_threshold,
        opts_.bad_channels,
        calo_alg,
        clock_data,
        detprop_data,
        T0_ticks
    };

    for (size_t wire_idx = 0; wire_idx < wire_vector->size(); ++wire_idx) {
        fillDetectorImage(wire_vector->at(wire_idx), wire_idx, hit_set, ctx);
    }
}

#endif
