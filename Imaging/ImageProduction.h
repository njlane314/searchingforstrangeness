#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

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

class PixelImageBuilder {
public:
    static void constructPixelImages(const art::Event &event,
                                     const std::vector<art::Ptr<recob::Hit>> &hits,
                                     const std::vector<ImageProperties> &properties,
                                     std::vector<Image<float>> &detector_images,
                                     std::vector<Image<int>> &semantic_images,
                                     bool is_data,
                                     const art::InputTag &wire_producer,
                                     const art::InputTag &hit_producer,
                                     const art::InputTag &mcp_producer,
                                     const art::InputTag &bkt_producer,
                                     const geo::GeometryCore *geo,
                                     const detinfo::DetectorProperties *detp,
                                     float adc_image_threshold,
                                     SemanticClassifier *semantic_classifier,
                                     const std::set<unsigned int> &bad_channels);

    static void constructPixelImages(const art::Event &event,
                                     const std::vector<art::Ptr<recob::Hit>> &hits,
                                     const std::vector<ImageProperties> &properties,
                                     std::vector<Image<float>> &detector_images,
                                     std::vector<Image<int>> &semantic_images,
                                     bool is_data,
                                     const art::InputTag &wire_producer,
                                     const art::InputTag &hit_producer,
                                     const art::InputTag &mcp_producer,
                                     const art::InputTag &bkt_producer,
                                     const geo::GeometryCore *geo,
                                     const detinfo::DetectorProperties *detp,
                                     float adc_image_threshold,
                                     SemanticClassifier *semantic_classifier,
                                     const std::set<unsigned int> &bad_channels,
                                     calo::CalorimetryAlg* calo_alg,
                                     detinfo::DetectorClocksData const* clock_data,
                                     detinfo::DetectorPropertiesData const* detprop_data,
                                     double T0_ticks);

private:
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
        const std::vector<ImageProperties> &properties,
        const art::FindManyP<recob::Hit> &wire_hit_assoc,
        const std::set<art::Ptr<recob::Hit>> &hit_set,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
        std::vector<Image<float>> &detector_images,
        std::vector<Image<int>> &semantic_images,
        bool is_data,
        bool has_mcps,
        const geo::GeometryCore *geo,
        const detinfo::DetectorProperties *detp,
        float adc_image_threshold,
        const std::set<unsigned int> &bad_channels,
        calo::CalorimetryAlg* calo_alg,
        detinfo::DetectorClocksData const* clock_data,
        detinfo::DetectorPropertiesData const* detprop_data,
        double T0_ticks);
};

} // namespace image

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
    const std::vector<ImageProperties> &properties,
    const art::FindManyP<recob::Hit> &wire_hit_assoc,
    const std::set<art::Ptr<recob::Hit>> &hit_set,
    const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
        &mcp_bkth_assoc,
    const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
    const std::map<int, size_t> &trackid_to_index,
    const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images,
    bool is_data,
    bool has_mcps,
    const geo::GeometryCore *geo,
    const detinfo::DetectorProperties *detp,
    float adc_image_threshold,
    const std::set<unsigned int> &bad_channels,
    calo::CalorimetryAlg* calo_alg,
    detinfo::DetectorClocksData const* clock_data,
    detinfo::DetectorPropertiesData const* detprop_data,
    double T0_ticks) {
    auto ch_id = wire.Channel();
    if (bad_channels.count(ch_id)) {
        return;
    }

    std::vector<geo::WireID> wire_ids = geo->ChannelToWire(ch_id);
    if (wire_ids.empty()) {
        return;
    }

    geo::View_t view = geo->View(wire_ids.front().planeID());
    size_t view_idx = static_cast<size_t>(view);
    const geo::WireGeo *wire_geo = geo->WirePtr(wire_ids.front());
    TVector3 center = wire_geo->GetCenter();
    TVector3 wire_center(center.X(), center.Y(), center.Z());

    double wire_coord =
        (view == geo::kW)
            ? wire_center.Z()
        : (view == geo::kU)
            ? (wire_center.Z() * std::cos(1.04719758034) -
               wire_center.Y() * std::sin(1.04719758034))
            : (wire_center.Z() * std::cos(-1.04719758034) -
               wire_center.Y() * std::sin(-1.04719758034));

    auto hits_for_wire = wire_hit_assoc.at(wire_idx);
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
            double x_drift_pos = detp->ConvertTicksToX(static_cast<double>(tick),
                                                       wire_ids.front().planeID());

            auto row = properties[view_idx].row(x_drift_pos);
            auto col = properties[view_idx].col(wire_coord);
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
                if (adcs[adc_index] > adc_image_threshold) {
                    float out_val = adcs[adc_index];
                    if (calo_alg && clock_data && detprop_data) {
                        const double pitch_cm_drift = properties[view_idx].pixel_w();
                        if (pitch_cm_drift > 0.) {
                            const double dQdx_ADC = static_cast<double>(adcs[adc_index]) / pitch_cm_drift;
                            const unsigned plane = wire_ids.front().planeID().Plane;
                            const double time_tick = static_cast<double>(tick);
                            const double dEdx_MeV_cm =
                                calo_alg->dEdx_AMP(*clock_data, *detprop_data,
                                                   dQdx_ADC, time_tick, plane,
                                                   /*T0=*/T0_ticks);
                            out_val = static_cast<float>(dEdx_MeV_cm);
                        }
                    }
                    detector_images[view_idx].set(*row, *col, out_val);

                    if (!is_data) {
                        auto semantic_pixel_label = labelSemanticPixels(
                            filtered_hits[hit_index], mcp_bkth_assoc, mcp_vector,
                            trackid_to_index, semantic_label_vector, has_mcps);

                        semantic_images[view_idx].set(
                            *row, *col,
                            static_cast<int>(semantic_pixel_label), false);
                    }
                }
            }
        }
    }
}

inline void image::PixelImageBuilder::constructPixelImages(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<ImageProperties> &properties,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images,
    bool is_data,
    const art::InputTag &wire_producer,
    const art::InputTag &hit_producer,
    const art::InputTag &mcp_producer,
    const art::InputTag &bkt_producer,
    const geo::GeometryCore *geo,
    const detinfo::DetectorProperties *detp,
    float adc_image_threshold,
    SemanticClassifier *semantic_classifier,
    const std::set<unsigned int> &bad_channels) {
    constructPixelImages(event, hits, properties, detector_images, semantic_images,
                         is_data, wire_producer, hit_producer, mcp_producer, bkt_producer,
                         geo, detp, adc_image_threshold, semantic_classifier, bad_channels,
                         /*calo_alg*/nullptr, /*clock_data*/nullptr, /*detprop_data*/nullptr,
                         /*T0_ticks*/0.0);
}

inline void image::PixelImageBuilder::constructPixelImages(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<ImageProperties> &properties,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images,
    bool is_data,
    const art::InputTag &wire_producer,
    const art::InputTag &hit_producer,
    const art::InputTag &mcp_producer,
    const art::InputTag &bkt_producer,
    const geo::GeometryCore *geo,
    const detinfo::DetectorProperties *detp,
    float adc_image_threshold,
    SemanticClassifier *semantic_classifier,
    const std::set<unsigned int> &bad_channels,
    calo::CalorimetryAlg* calo_alg,
    detinfo::DetectorClocksData const* clock_data,
    detinfo::DetectorPropertiesData const* detprop_data,
    double T0_ticks) {
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
        event.getValidHandle<std::vector<recob::Wire>>(wire_producer);
    auto hit_vector =
        event.getValidHandle<std::vector<recob::Hit>>(hit_producer);
    art::Handle<std::vector<simb::MCParticle>> mcp_vector;
    bool has_mcps = event.getByLabel(mcp_producer, mcp_vector);
    art::FindManyP<recob::Hit> wire_hit_assoc(wire_vector, event, hit_producer);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
        mcp_bkth_assoc(hit_vector, event, bkt_producer);
    std::vector<SemanticClassifier::SemanticLabel> semantic_label_vector;
    if (!is_data && has_mcps && mcp_vector.isValid() && semantic_classifier) {
        semantic_label_vector = semantic_classifier->classifyParticles(event);
    }
    std::map<int, size_t> trackid_to_index;
    if (!is_data && has_mcps && mcp_vector.isValid()) {
        for (size_t i = 0; i < mcp_vector->size(); ++i) {
            trackid_to_index[mcp_vector->at(i).TrackId()] = i;
        }
    }

    std::set<art::Ptr<recob::Hit>> hit_set(hits.begin(), hits.end());
    for (size_t wire_idx = 0; wire_idx < wire_vector->size(); ++wire_idx) {
        fillDetectorImage(wire_vector->at(wire_idx), wire_idx, properties,
                          wire_hit_assoc, hit_set, mcp_bkth_assoc, mcp_vector,
                          trackid_to_index, semantic_label_vector,
                          detector_images, semantic_images, is_data, has_mcps,
                          geo, detp, adc_image_threshold, bad_channels,
                          calo_alg, clock_data, detprop_data, T0_ticks);
    }
}


#endif // IMAGEPRODUCTION_H
