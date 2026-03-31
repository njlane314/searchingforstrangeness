#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

#include <algorithm>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <vector>

#include "Helpers/PandoraUtilities.h"
#include "ImagePipeline/Image.h"
#include "ImagePipeline/SemanticClassifier.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataalg/DetectorInfo/DetectorProperties.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include <TVector3.h>
#include <cetlib_except/exception.h>

namespace image {

struct PixelImageOptions {
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
    } producers;
    image::SemanticClassifier *semantic{nullptr};
};

class ImageProduction {
  public:
    ImageProduction(geo::GeometryCore const &geo, PixelImageOptions const &opts) : geo_{&geo}, opts_{opts} {}

    void build(const art::Event &event, const std::vector<art::Ptr<recob::Hit>> &hits,
               const std::vector<art::Ptr<recob::Hit>> &slice_hits, const std::vector<ImageProperties> &properties,
               std::vector<Image<float>> &detector_images, std::vector<Image<int>> &semantic_images,
               std::vector<Image<uint8_t>> &slice_images,
               detinfo::DetectorProperties const *detprop) const {
        detector_images.clear();
        semantic_images.clear();
        slice_images.clear();

        for (auto const &p : properties) {
            Image<float> det(p);
            det.clear(0.0f);
            detector_images.push_back(std::move(det));

            Image<int> sem(p);
            sem.clear(static_cast<int>(image::SemanticClassifier::SemanticLabel::Empty));
            semantic_images.push_back(std::move(sem));

            Image<uint8_t> slice(p);
            slice.clear(0);
            slice_images.push_back(std::move(slice));
        }

        auto wires = event.getValidHandle<std::vector<recob::Wire>>(opts_.producers.wire);
        auto all_hits = event.getValidHandle<std::vector<recob::Hit>>(opts_.producers.hit);

        art::Handle<std::vector<simb::MCParticle>> mcps;
        bool has_mcps = event.getByLabel(opts_.producers.mcp, mcps);

        art::FindManyP<recob::Hit> wire_hit_assoc(wires, event, opts_.producers.hit);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(all_hits, event, opts_.producers.bkt);

        std::vector<image::SemanticClassifier::SemanticLabel> sem_labels;
        if (has_mcps && mcps.isValid() && opts_.semantic)
            sem_labels = opts_.semantic->classifyParticles(event);

        std::map<int, std::size_t> trkid_to_idx;
        if (has_mcps && mcps.isValid()) {
            for (std::size_t i = 0; i < mcps->size(); ++i)
                trkid_to_idx[mcps->at(i).TrackId()] = i;
        }

        std::map<art::Ptr<recob::Hit>, std::size_t> hit_to_key;
        for (auto const &ph : hits)
            hit_to_key.emplace(ph, static_cast<std::size_t>(ph.key()));

        std::set<std::size_t> slice_hit_keys;
        for (auto const &ph : slice_hits) {
            if (ph)
                slice_hit_keys.insert(static_cast<std::size_t>(ph.key()));
        }

        BuildContext ctx{properties,   detector_images, semantic_images, slice_images, wire_hit_assoc, hit_to_key,
                         slice_hit_keys, mcp_bkth_assoc, mcps, trkid_to_idx, sem_labels, has_mcps, geo_,
                         kAdcThreshold, detprop};

        if (!ctx.detprop) {
            throw cet::exception("ImageProduction") << "ImageProduction::build called with null DetectorProperties pointer.";
        }

        for (std::size_t wi = 0; wi < wires->size(); ++wi) {
            fillImagesForWire(wires->at(wi), wi, ctx);
        }
    }

  private:
    static constexpr float kAdcThreshold = 0.0f;

    struct BuildContext {
        const std::vector<ImageProperties> &properties;
        std::vector<Image<float>> &detector_images;
        std::vector<Image<int>> &semantic_images;
        std::vector<Image<uint8_t>> &slice_images;

        const art::FindManyP<recob::Hit> &wire_hit_assoc;
        const std::map<art::Ptr<recob::Hit>, std::size_t> &hit_to_key;
        const std::set<std::size_t> &slice_hit_keys;

        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc;
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector;
        const std::map<int, std::size_t> &trackid_to_index;
        const std::vector<image::SemanticClassifier::SemanticLabel> &semantic_label_vector;

        bool has_mcps{false};

        const geo::GeometryCore *geo{nullptr};

        float adc_image_threshold{0.f};

        detinfo::DetectorProperties const *detprop{nullptr};
    };

    static image::SemanticClassifier::SemanticLabel
    labelSemanticPixels(const art::Ptr<recob::Hit> &matched_hit, std::size_t matched_hit_key,
                        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
                        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector, const std::map<int, std::size_t> &trackid_to_index,
                        const std::vector<image::SemanticClassifier::SemanticLabel> &semantic_label_vector, bool has_mcps) {
        image::SemanticClassifier::SemanticLabel out = image::SemanticClassifier::SemanticLabel::Cosmic;
        if (!has_mcps || !mcp_vector.isValid())
            return out;

        std::vector<art::Ptr<simb::MCParticle>> parts;
        std::vector<anab::BackTrackerHitMatchingData const *> bkd;
        mcp_bkth_assoc.get(matched_hit_key, parts, bkd);
        if (bkd.empty())
            return out;

        float max_ide_fraction = -1.0f;
        int best_trkid = -1;
        for (std::size_t i = 0; i < bkd.size(); ++i) {
            if (bkd[i] && bkd[i]->ideFraction > max_ide_fraction) {
                max_ide_fraction = bkd[i]->ideFraction;
                best_trkid = parts[i]->TrackId();
            }
        }
        if (best_trkid == -1)
            return out;

        auto it = trackid_to_index.find(best_trkid);
        if (it == trackid_to_index.end())
            return out;
        const std::size_t idx = it->second;
        if (idx >= semantic_label_vector.size())
            return out;

        if (max_ide_fraction <= 0.5f)
            return image::SemanticClassifier::SemanticLabel::Ambiguous;
        return semantic_label_vector[idx];
    }

    static void mergeSemanticPixel(Image<int> &image,
                                   size_t row, size_t col,
                                   image::SemanticClassifier::SemanticLabel label) {
        int const empty =
            static_cast<int>(image::SemanticClassifier::SemanticLabel::Empty);
        int const ambiguous =
            static_cast<int>(image::SemanticClassifier::SemanticLabel::Ambiguous);
        int const next = static_cast<int>(label);
        int const current = image.get(row, col);

        if (current == empty) {
            image.set(row, col, next, false);
            return;
        }
        if (current != next && current != ambiguous) {
            image.set(row, col, ambiguous, false);
        }
    }

    struct WirePrep {
        geo::PlaneID planeID;
        std::size_t view_idx;
        TVector3 wire_center;
        double wire_coord{0.0};
        std::vector<art::Ptr<recob::Hit>> hits_filtered;
    };

    static std::optional<WirePrep> prepareWire(const recob::Wire &wire, std::size_t wire_idx, BuildContext const &ctx) {
        auto hits_for_wire = ctx.wire_hit_assoc.at(wire_idx);

        std::vector<art::Ptr<recob::Hit>> hits_filtered;
        hits_filtered.reserve(hits_for_wire.size());
        for (auto const &ph : hits_for_wire) {
            if (ctx.hit_to_key.find(ph) != ctx.hit_to_key.end())
                hits_filtered.push_back(ph);
        }
        if (hits_filtered.empty())
            return std::nullopt;

        recob::Hit const &firstHit = *hits_filtered.front();
        geo::WireID const &wid = firstHit.WireID();
        geo::PlaneID const planeID = wid.planeID();

        geo::View_t const view = ctx.geo->View(planeID);

        std::size_t view_idx = 0;
        bool found_view_idx = false;
        for (std::size_t i = 0; i < ctx.properties.size(); ++i) {
            if (ctx.properties[i].view() == view) {
                view_idx = i;
                found_view_idx = true;
                break;
            }
        }
        if (!found_view_idx)
            return std::nullopt;

        geo::WireGeo const &wire_geo = ctx.geo->WireIDToWireGeo(wid);
        TVector3 const wire_center = wire_geo.GetCenter();

        common::PandoraView pandora_view;
        switch (view) {
        case geo::kU:
            pandora_view = common::TPC_VIEW_U;
            break;
        case geo::kV:
            pandora_view = common::TPC_VIEW_V;
            break;
        case geo::kW:
        case geo::kY:
            pandora_view = common::TPC_VIEW_W;
            break;
        default:
            throw cet::exception("ImageProduction") << "Unsupported geo::View_t: " << static_cast<int>(view);
        }

        TVector3 const proj = common::ProjectToWireView(static_cast<float>(wire_center.X()), static_cast<float>(wire_center.Y()),
                                                        static_cast<float>(wire_center.Z()), pandora_view);

        return WirePrep{planeID, view_idx, wire_center, proj.Z(), std::move(hits_filtered)};
    }

    static void fillImagesForWire(const recob::Wire &wire, std::size_t wire_idx, BuildContext const &ctx) {
        auto prep = prepareWire(wire, wire_idx, ctx);
        if (!prep)
            return;
        const auto &w = *prep;

        auto const col = ctx.properties[w.view_idx].col(w.wire_coord);
        if (!col)
            return;

        int min_tick = std::numeric_limits<int>::max();
        int max_tick = std::numeric_limits<int>::min();
        for (auto const &ph : w.hits_filtered) {
            if (!ph)
                continue;
            min_tick = std::min(min_tick, ph->StartTick());
            max_tick = std::max(max_tick, ph->EndTick());
        }
        if (min_tick >= max_tick)
            return;

        auto const nticks = static_cast<std::size_t>(max_tick - min_tick);
        std::vector<unsigned char> tick_active(nticks, 0);
        std::vector<unsigned char> tick_in_slice(nticks, 0);

        int const empty_label =
            static_cast<int>(image::SemanticClassifier::SemanticLabel::Empty);
        int const ambiguous_label =
            static_cast<int>(image::SemanticClassifier::SemanticLabel::Ambiguous);
        std::vector<int> tick_semantic;
        if (ctx.has_mcps)
            tick_semantic.assign(nticks, empty_label);

        for (auto const &ph : w.hits_filtered) {
            const recob::Hit &hit = *ph;
            bool const hit_in_slice =
                ctx.slice_hit_keys.count(static_cast<std::size_t>(ph.key())) > 0;

            image::SemanticClassifier::SemanticLabel sem = image::SemanticClassifier::SemanticLabel::Cosmic;
            if (ctx.has_mcps) {
                auto hit_it = ctx.hit_to_key.find(ph);
                if (hit_it == ctx.hit_to_key.end())
                    continue;
                sem = labelSemanticPixels(ph, hit_it->second, ctx.mcp_bkth_assoc, ctx.mcp_vector, ctx.trackid_to_index,
                                          ctx.semantic_label_vector, ctx.has_mcps);
            }

            int const tick_start = std::max(hit.StartTick(), min_tick);
            int const tick_end = std::min(hit.EndTick(), max_tick);
            if (tick_start >= tick_end)
                continue;

            for (int t = tick_start; t < tick_end; ++t) {
                auto const ti = static_cast<std::size_t>(t - min_tick);
                tick_active[ti] = 1;
                if (hit_in_slice)
                    tick_in_slice[ti] = 1;
                if (ctx.has_mcps) {
                    int &slot = tick_semantic[ti];
                    int const next = static_cast<int>(sem);
                    if (slot == empty_label)
                        slot = next;
                    else if (slot != next)
                        slot = ambiguous_label;
                }
            }
        }

        for (auto const &rr : wire.SignalROI().get_ranges()) {
            const int rbeg = rr.begin_index();
            const auto &adcs = rr.data();
            const int rend = rbeg + static_cast<int>(adcs.size());
            const int s = std::max(rbeg, min_tick);
            const int e = std::min(rend, max_tick);
            if (s >= e)
                continue;

            for (int t = s; t < e; ++t) {
                auto const ti = static_cast<std::size_t>(t - min_tick);
                if (!tick_active[ti])
                    continue;

                const float a = adcs[t - rbeg];
                if (a <= ctx.adc_image_threshold)
                    continue;

                double const x = ctx.detprop->ConvertTicksToX(static_cast<double>(t), w.planeID);
                auto row = ctx.properties[w.view_idx].row(x);
                if (!row)
                    continue;

                ctx.detector_images[w.view_idx].set(*row, *col, a, true);
                if (tick_in_slice[ti])
                    ctx.slice_images[w.view_idx].set(*row, *col, static_cast<uint8_t>(1), false);
                if (ctx.has_mcps && tick_semantic[ti] != empty_label) {
                    mergeSemanticPixel(
                        ctx.semantic_images[w.view_idx], *row, *col,
                        static_cast<image::SemanticClassifier::SemanticLabel>(tick_semantic[ti]));
                }
            }
        }
    }

    geo::GeometryCore const *geo_{nullptr};
    PixelImageOptions opts_;
};

} // namespace image

#endif
