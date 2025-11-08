#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

#include <optional>
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

#include "Imaging/Image.h"
#include "Imaging/SemanticClassifier.h"
#include "Imaging/ImageCorrections.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

// Use services in the header; we only need pointers to the Data types here.
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
namespace detinfo { class DetectorClocksData; class DetectorPropertiesData; class DetectorProperties; }

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationServices/TPCEnergyCalibService.h"
#include "larevt/CalibrationServices/ChannelStatusService.h"

#include <TVector3.h>

namespace image {

struct PixelImageOptions {
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
    } producers;
    sem::SemanticClassifier* semantic{nullptr};
};

struct CalibrationContext {
    calo::CalorimetryAlg* calo{nullptr};
    detinfo::DetectorClocksData const* clocks{nullptr};
    detinfo::DetectorPropertiesData const* detprop{nullptr};

    lariov::TPCEnergyCalib const* tpcCalib{nullptr};
    spacecharge::SpaceCharge const* sce{nullptr};
    lariov::ChannelStatusProvider const* chanStatus{nullptr};
    double T0_ticks{0.0};

    bool enabled() const { return calo && clocks && detprop; }
};

class ImageProduction {
public:
    ImageProduction(geo::GeometryCore const& geo,
                    PixelImageOptions const& opts)
        : geo_{&geo}, opts_{opts} {}

    void build(const art::Event &event,
               const std::vector<art::Ptr<recob::Hit>> &hits,
               const std::vector<ImageProperties> &properties,
               std::vector<Image<float>> &detector_images,
               std::vector<Image<int>> &semantic_images,
               const detinfo::DetectorProperties *,
               std::optional<CalibrationContext> const& cal = std::nullopt) const
    {
        detector_images.clear();
        semantic_images.clear();

        for (auto const& p : properties) {
            Image<float> det(p); det.clear(0.0f); detector_images.push_back(std::move(det));
            Image<int>   sem(p); sem.clear(static_cast<int>(sem::SemanticClassifier::SemanticLabel::Empty));
            semantic_images.push_back(std::move(sem));
        }

        auto wires = event.getValidHandle<std::vector<recob::Wire>>(opts_.producers.wire);
        auto all_hits = event.getValidHandle<std::vector<recob::Hit>>(opts_.producers.hit);

        art::Handle<std::vector<simb::MCParticle>> mcps;
        bool has_mcps = event.getByLabel(opts_.producers.mcp, mcps);

        art::FindManyP<recob::Hit> wire_hit_assoc(wires, event, opts_.producers.hit);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
            mcp_bkth_assoc(all_hits, event, opts_.producers.bkt);

        std::vector<sem::SemanticClassifier::SemanticLabel> sem_labels;
        if (has_mcps && mcps.isValid() && opts_.semantic)
            sem_labels = opts_.semantic->classifyParticles(event);

        std::map<int, size_t> trkid_to_idx;
        if (has_mcps && mcps.isValid())
            for (size_t i = 0; i < mcps->size(); ++i) trkid_to_idx[mcps->at(i).TrackId()] = i;


        std::map<art::Ptr<recob::Hit>, std::size_t> hit_to_key;
        for (auto const& ph : hits)
            hit_to_key.emplace(ph, static_cast<std::size_t>(ph.key()));
        BuildContext ctx{
            properties,
            detector_images,
            semantic_images,
            wire_hit_assoc,
            hit_to_key,
            mcp_bkth_assoc,
            mcps,
            trkid_to_idx,
            sem_labels,
            has_mcps,
            geo_,
            kAdcThreshold,
            (cal && cal->enabled()) ? cal->calo      : nullptr,
            (cal && cal->enabled()) ? cal->clocks    : nullptr,
            (cal && cal->enabled()) ? cal->detprop   : nullptr,
            (cal && cal->enabled()) ? cal->tpcCalib  : nullptr,
            (cal && cal->enabled()) ? cal->sce       : nullptr,
            (cal && cal->enabled()) ? cal->chanStatus: nullptr,
            (cal && cal->enabled()) ? cal->T0_ticks  : 0.0
        };


        for (size_t wi = 0; wi < wires->size(); ++wi) {
            fillImagesForWire(wires->at(wi), wi, ctx);
        }

        smoothDetectorImages(ctx);
    }

private:
    static constexpr float kAdcThreshold    = 4.0f;
    static constexpr float kGaussianSigmaPx = 1.0f;  

    struct BuildContext {
        const std::vector<ImageProperties> &properties;
        std::vector<Image<float>> &detector_images;
        std::vector<Image<int>> &semantic_images;

        const art::FindManyP<recob::Hit> &wire_hit_assoc;
        const std::map<art::Ptr<recob::Hit>, std::size_t> &hit_to_key;

        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc;
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector;
        const std::map<int, size_t> &trackid_to_index;
        const std::vector<sem::SemanticClassifier::SemanticLabel> &semantic_label_vector;

        bool has_mcps;

        const geo::GeometryCore *geo;

        float adc_image_threshold;

        calo::CalorimetryAlg* calo_alg;
        detinfo::DetectorClocksData const* clocks;
        detinfo::DetectorPropertiesData const* detprop_data;
        lariov::TPCEnergyCalib const* tpcCalib;
        spacecharge::SpaceCharge const* sce;
        lariov::ChannelStatusProvider const* chanStatus;

        double T0_ticks{0.0};
    };

    static sem::SemanticClassifier::SemanticLabel labelSemanticPixels(
        const art::Ptr<recob::Hit> &matched_hit,
        std::size_t matched_hit_key,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<sem::SemanticClassifier::SemanticLabel> &semantic_label_vector,
        bool has_mcps)
    {
        sem::SemanticClassifier::SemanticLabel out = sem::SemanticClassifier::SemanticLabel::Cosmic;
        if (!has_mcps || !mcp_vector.isValid()) return out;

        std::vector<art::Ptr<simb::MCParticle>> parts;
        std::vector<anab::BackTrackerHitMatchingData const*> bkd;
        mcp_bkth_assoc.get(matched_hit_key, parts, bkd);
        if (bkd.empty()) return out;

        float max_ide_fraction = -1.0f;
        int best_trkid = -1;
        for (size_t i = 0; i < bkd.size(); ++i) {
            if (bkd[i] && bkd[i]->ideFraction > max_ide_fraction) {
                max_ide_fraction = bkd[i]->ideFraction;
                best_trkid = parts[i]->TrackId();
            }
        }
        if (best_trkid == -1) return out;

        auto it = trackid_to_index.find(best_trkid);
        if (it == trackid_to_index.end()) return out;
        const size_t idx = it->second;
        if (idx >= semantic_label_vector.size()) return out;

        if (max_ide_fraction <= 0.5f) return sem::SemanticClassifier::SemanticLabel::Ambiguous;
        return semantic_label_vector[idx];
    }

    struct WirePrep {
        geo::PlaneID planeID;
        size_t view_idx;
        TVector3 wire_center;
        std::vector<art::Ptr<recob::Hit>> hits_filtered;
    };

    static std::optional<WirePrep> prepareWire(
        const recob::Wire &wire,
        size_t wire_idx,
        BuildContext const& ctx)
    {
        const unsigned ch = wire.Channel();
        if (ctx.chanStatus && ctx.chanStatus->IsBad(ch)) return std::nullopt;

        auto wire_ids = ctx.geo->ChannelToWire(ch);
        if (wire_ids.empty()) return std::nullopt;

        const auto planeID = wire_ids.front().planeID();
        const geo::View_t view = ctx.geo->View(planeID);
        const size_t view_idx = static_cast<size_t>(view);
        if (view_idx >= ctx.properties.size()) return std::nullopt;

        const geo::WireGeo* wire_geo = ctx.geo->WirePtr(wire_ids.front());
        const TVector3 wire_center = wire_geo->GetCenter();

        auto hits_for_wire = ctx.wire_hit_assoc.at(wire_idx);
        std::vector<art::Ptr<recob::Hit>> hits_filtered;
        hits_filtered.reserve(hits_for_wire.size());
        for (auto const& ph : hits_for_wire)
            if (ctx.hit_to_key.find(ph) != ctx.hit_to_key.end()) hits_filtered.push_back(ph);
        if (hits_filtered.empty()) return std::nullopt;

        return WirePrep{planeID, view_idx, wire_center, std::move(hits_filtered)};
    }

    static void fillImagesForWire(
        const recob::Wire &wire,
        size_t wire_idx,
        BuildContext const& ctx)
    {
        auto prep = prepareWire(wire, wire_idx, ctx);
        if (!prep) return;
        const auto& w = *prep;

        for (auto const& ph : w.hits_filtered) {
            const recob::Hit& hit = *ph;
            const unsigned plane = w.planeID.Plane;
            const int tick_c = static_cast<int>(hit.PeakTime());

            auto geo_res = cal::applyGeometry(ctx.detprop_data, ctx.sce, w.planeID,
                                              tick_c, w.wire_center, ctx.properties[w.view_idx]);
            if (!geo_res.col) continue;

            double pitch_cm = 1.0;
            if (ctx.geo) pitch_cm = std::max(1e-6, ctx.geo->Plane(w.planeID).WirePitch());

            auto calo_res = cal::applyCalorimetry(hit, plane, geo_res.p_corr, pitch_cm,
                                                  ctx.calo_alg, ctx.clocks, ctx.detprop_data,
                                                  ctx.tpcCalib, ctx.sce, ctx.T0_ticks);

            sem::SemanticClassifier::SemanticLabel sem = sem::SemanticClassifier::SemanticLabel::Cosmic;
            if (ctx.has_mcps) {
                auto hit_it = ctx.hit_to_key.find(ph);
                if (hit_it == ctx.hit_to_key.end()) continue;
                sem = labelSemanticPixels(ph, hit_it->second, ctx.mcp_bkth_assoc, ctx.mcp_vector,
                                          ctx.trackid_to_index, ctx.semantic_label_vector,
                                          ctx.has_mcps);
            }

            const int tick_start = hit.StartTick();
            const int tick_end   = hit.EndTick();

            struct PixelContribution {
                int row;
                float adc;
            };
            std::vector<PixelContribution> contributions;
            contributions.reserve(std::max(0, tick_end - tick_start));

            double sumw = 0.0;

            for (auto const& rr : wire.SignalROI().get_ranges()) {
                const int rbeg = rr.begin_index();
                const auto& adcs = rr.data();
                const int rend = rbeg + static_cast<int>(adcs.size());
                const int s = std::max(tick_start, rbeg);
                const int e = std::min(tick_end,   rend);
                if (s >= e) continue;

                for (int t = s; t < e; ++t) {
                    const float a = adcs[t - rbeg];
                    if (a <= ctx.adc_image_threshold) continue;

                    auto row = geo_res.row(t);
                    if (!row) continue;

                    contributions.push_back({*row, a});
                    sumw += static_cast<double>(a);
                }
            }

            if (contributions.empty()) continue;

            if (sumw <= 0.0)
                sumw = static_cast<double>(std::max(1, tick_end - tick_start));

            for (auto const& contrib : contributions) {
                const double wgt = (sumw > 0.0) ? (static_cast<double>(contrib.adc) / sumw) : 0.0;
                const float E_pix = static_cast<float>(calo_res.E_hit_MeV * wgt);

                ctx.detector_images[w.view_idx].set(contrib.row, *geo_res.col, E_pix, true);
                if (ctx.has_mcps)
                    ctx.semantic_images[w.view_idx].set(contrib.row, *geo_res.col, static_cast<int>(sem), false);
            }
        }
    }

    static void smoothDetectorImages(BuildContext const& ctx)
    {
        for (size_t view_idx = 0; view_idx < ctx.detector_images.size(); ++view_idx) {
            ctx.detector_images[view_idx].blur(kGaussianSigmaPx);
            auto& sem_img = ctx.semantic_images[view_idx];
            auto& det_img = ctx.detector_images[view_idx];
            const int empty =
                static_cast<int>(sem::SemanticClassifier::SemanticLabel::Empty);
            sem_img.dilate(det_img,
                           kGaussianSigmaPx,
                           empty);
        }
    }

    geo::GeometryCore const* geo_{nullptr};
    PixelImageOptions opts_;
};

}

#endif
