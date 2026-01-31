#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

#include <algorithm>
#include <map>
#include <optional>
#include <vector>

#include "Common/PandoraUtilities.h"
#include "Imaging/Image.h"
#include "Imaging/SemanticClassifier.h"

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
    sem::SemanticClassifier *semantic{nullptr};
};

class ImageProduction {
  public:
    ImageProduction(geo::GeometryCore const &geo, PixelImageOptions const &opts) : geo_{&geo}, opts_{opts} {}

    void build(const art::Event &event, const std::vector<art::Ptr<recob::Hit>> &hits, const std::vector<ImageProperties> &properties,
               std::vector<Image<float>> &detector_images, std::vector<Image<int>> &semantic_images,
               detinfo::DetectorProperties const *detprop) const {
        detector_images.clear();
        semantic_images.clear();

        for (auto const &p : properties) {
            Image<float> det(p);
            det.clear(0.0f);
            detector_images.push_back(std::move(det));

            Image<int> sem(p);
            sem.clear(static_cast<int>(sem::SemanticClassifier::SemanticLabel::Empty));
            semantic_images.push_back(std::move(sem));
        }

        auto wires = event.getValidHandle<std::vector<recob::Wire>>(opts_.producers.wire);
        auto all_hits = event.getValidHandle<std::vector<recob::Hit>>(opts_.producers.hit);

        art::Handle<std::vector<simb::MCParticle>> mcps;
        bool has_mcps = event.getByLabel(opts_.producers.mcp, mcps);

        art::FindManyP<recob::Hit> wire_hit_assoc(wires, event, opts_.producers.hit);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(all_hits, event, opts_.producers.bkt);

        std::vector<sem::SemanticClassifier::SemanticLabel> sem_labels;
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

        BuildContext ctx{properties,   detector_images, semantic_images, wire_hit_assoc, hit_to_key,    mcp_bkth_assoc, mcps,
                         trkid_to_idx, sem_labels,      has_mcps,        geo_,           kAdcThreshold, detprop};

        if (!ctx.detprop) {
            throw cet::exception("ImageProduction") << "ImageProduction::build called with null DetectorProperties pointer.";
        }

        for (std::size_t wi = 0; wi < wires->size(); ++wi) {
            fillImagesForWire(wires->at(wi), wi, ctx);
        }
    }

  private:
    static constexpr float kAdcThreshold = 4.0f;

    struct BuildContext {
        const std::vector<ImageProperties> &properties;
        std::vector<Image<float>> &detector_images;
        std::vector<Image<int>> &semantic_images;

        const art::FindManyP<recob::Hit> &wire_hit_assoc;
        const std::map<art::Ptr<recob::Hit>, std::size_t> &hit_to_key;

        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc;
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector;
        const std::map<int, std::size_t> &trackid_to_index;
        const std::vector<sem::SemanticClassifier::SemanticLabel> &semantic_label_vector;

        bool has_mcps{false};

        const geo::GeometryCore *geo{nullptr};

        float adc_image_threshold{0.f};

        detinfo::DetectorProperties const *detprop{nullptr};
    };

    static sem::SemanticClassifier::SemanticLabel
    labelSemanticPixels(const art::Ptr<recob::Hit> &matched_hit, std::size_t matched_hit_key,
                        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
                        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector, const std::map<int, std::size_t> &trackid_to_index,
                        const std::vector<sem::SemanticClassifier::SemanticLabel> &semantic_label_vector, bool has_mcps) {
        sem::SemanticClassifier::SemanticLabel out = sem::SemanticClassifier::SemanticLabel::Cosmic;
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
            return sem::SemanticClassifier::SemanticLabel::Ambiguous;
        return semantic_label_vector[idx];
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

        for (auto const &ph : w.hits_filtered) {
            const recob::Hit &hit = *ph;

            sem::SemanticClassifier::SemanticLabel sem = sem::SemanticClassifier::SemanticLabel::Cosmic;
            if (ctx.has_mcps) {
                auto hit_it = ctx.hit_to_key.find(ph);
                if (hit_it == ctx.hit_to_key.end())
                    continue;
                sem = labelSemanticPixels(ph, hit_it->second, ctx.mcp_bkth_assoc, ctx.mcp_vector, ctx.trackid_to_index,
                                          ctx.semantic_label_vector, ctx.has_mcps);
            }

            const int tick_start = hit.StartTick();
            const int tick_end = hit.EndTick();

            for (auto const &rr : wire.SignalROI().get_ranges()) {
                const int rbeg = rr.begin_index();
                const auto &adcs = rr.data();
                const int rend = rbeg + static_cast<int>(adcs.size());
                const int s = std::max(tick_start, rbeg);
                const int e = std::min(tick_end, rend);
                if (s >= e)
                    continue;

                for (int t = s; t < e; ++t) {
                    const float a = adcs[t - rbeg];
                    if (a <= ctx.adc_image_threshold)
                        continue;

                    double const x = ctx.detprop->ConvertTicksToX(static_cast<double>(t), w.planeID);
                    auto row = ctx.properties[w.view_idx].row(x);
                    if (!row)
                        continue;

                    ctx.detector_images[w.view_idx].set(*row, *col, a, true);
                    if (ctx.has_mcps)
                        ctx.semantic_images[w.view_idx].set(*row, *col, static_cast<int>(sem), false);
                }
            }
        }
    }

    geo::GeometryCore const *geo_{nullptr};
    PixelImageOptions opts_;
};

} // namespace image

#endif
