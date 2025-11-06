#ifndef IMAGEPRODUCTION_H
#define IMAGEPRODUCTION_H

#include <optional>
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "Imaging/Image.h"
#include "Imaging/SemanticClassifier.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "lardata/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfo/DetectorPropertiesData.h"
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationServices/TPCEnergyCalibService.h"
#include "larevt/CalibrationServices/ChannelStatusService.h"

#include <TVector3.h>

namespace image {

// Minimal, option-free inputs.
struct PixelImageOptions {
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
    } producers;
    image::SemanticClassifier* semantic{nullptr}; // optional
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

namespace cal {

inline geo::Point_t nominal_point_from_tick(detinfo::DetectorProperties const* detp_api,
                                            detinfo::DetectorPropertiesData const* detprop,
                                            geo::PlaneID const& planeID,
                                            double tick,
                                            TVector3 const& wire_center)
{
    (void)detp_api; 
    double x_nom = detprop->ConvertTicksToX(tick, planeID);
    return geo::Point_t{x_nom, wire_center.Y(), wire_center.Z()};
}

inline geo::Point_t correct_spacepoint(spacecharge::SpaceCharge const* sce,
                                       geo::Point_t p_in,
                                       bool always_apply)
{
    if (sce && sce->EnableCalSpatialSCE() && always_apply) {
        auto off = sce->GetCalPosOffsets(p_in, 0);
        return geo::Point_t{ p_in.X() - off.X(), p_in.Y() + off.Y(), p_in.Z() + off.Z() };
    }
    return p_in;
}

inline std::optional<size_t> tick_to_corrected_row(detinfo::DetectorProperties const* detp_api,
                                                   detinfo::DetectorPropertiesData const* detprop,
                                                   spacecharge::SpaceCharge const* sce,
                                                   geo::PlaneID const& planeID,
                                                   int tick,
                                                   TVector3 const& wire_center,
                                                   bool always_apply_sce,
                                                   ImageProperties const& prop)
{
    auto p_nom = nominal_point_from_tick(detp_api, detprop, planeID, static_cast<double>(tick), wire_center);
    auto p = correct_spacepoint(sce, p_nom, always_apply_sce);
    return prop.row(p.X());
}

inline double wire_column_coordinate(geo::View_t view, TVector3 const& wire_center)
{
    constexpr double plus60  =  1.04719758034;
    constexpr double minus60 = -1.04719758034;
    if (view == geo::kW) return wire_center.Z();
    if (view == geo::kU) return wire_center.Z() * std::cos(plus60)  - wire_center.Y() * std::sin(plus60);
    return                    wire_center.Z() * std::cos(minus60) - wire_center.Y() * std::sin(minus60);
}

inline double yz_gain_corr(lariov::TPCEnergyCalib const* tpcCalib,
                           unsigned plane,
                           geo::Point_t const& p,
                           bool always_apply)
{
    if (tpcCalib && always_apply) return tpcCalib->YZdqdxCorrection(plane, p.Y(), p.Z());
    return 1.0;
}

inline double local_efield_kV_cm(spacecharge::SpaceCharge const* sce,
                                 detinfo::DetectorPropertiesData const* detprop,
                                 geo::Point_t const& p,
                                 bool always_apply)
{
    double E = detprop->Efield();
    if (sce && sce->EnableCalEfieldSCE() && always_apply) {
        auto fo = sce->GetCalEfieldOffsets(p, 0);
        E *= std::hypot(1.0 + fo.X(), fo.Y(), fo.Z());
    }
    return E;
}

struct RangeRef { int begin; const std::vector<float>* adcs; };

inline double sum_adc_weights_in_window(std::vector<RangeRef> const& ranges,
                                        int tick_start, int tick_end,
                                        float adc_threshold)
{
    double sumw = 0.0;
    for (auto const& rr : ranges) {
        const int rbeg = rr.begin;
        const int rend = rbeg + static_cast<int>(rr.adcs->size());
        const int s = std::max(tick_start, rbeg);
        const int e = std::min(tick_end,   rend);
        if (s >= e) continue;
        for (int t = s; t < e; ++t) {
            const float a = (*rr.adcs)[t - rbeg];
            if (a > adc_threshold) sumw += static_cast<double>(a);
        }
    }
    if (sumw <= 0.0) sumw = static_cast<double>(std::max(1, tick_end - tick_start));
    return sumw;
}

} 

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
               const detinfo::DetectorProperties *detp,
               std::optional<CalibrationContext> const& cal = std::nullopt) const;

private:
    static constexpr int   kCollectionPlane = 2;   // minimal, fixed choice
    static constexpr float kAdcThreshold    = 4.0f; // minimal, fixed threshold

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

        bool has_mcps;

        const geo::GeometryCore *geo;
        const detinfo::DetectorProperties *detp_api;

        float adc_image_threshold;

        calo::CalorimetryAlg* calo_alg;
        detinfo::DetectorClocksData const* clocks;
        detinfo::DetectorPropertiesData const* detprop_data;
        lariov::TPCEnergyCalib const* tpcCalib;
        spacecharge::SpaceCharge const* sce;
        lariov::ChannelStatusProvider const* chan_status;

        double T0_ticks{0.0};
    };

    static SemanticClassifier::SemanticLabel labelSemanticPixels(
        const art::Ptr<recob::Hit> &matched_hit,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
        bool has_mcps);

    static void paint_hit_energy(
        recob::Hit const& hit,
        recob::Wire const& wire,
        geo::PlaneID const& planeID,
        geo::View_t view,
        TVector3 const& wire_center,
        size_t view_idx,
        std::vector<cal::RangeRef> const& roi_ranges,
        BuildContext const& ctx);

    static void fillDetectorImage(
        const recob::Wire &wire,
        size_t wire_idx,
        const std::set<art::Ptr<recob::Hit>> &hit_set,
        BuildContext const& ctx);

    geo::GeometryCore const* geo_{nullptr};
    PixelImageOptions opts_;
};

inline SemanticClassifier::SemanticLabel
ImageProduction::labelSemanticPixels(
    const art::Ptr<recob::Hit> &matched_hit,
    const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
    const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
    const std::map<int, size_t> &trackid_to_index,
    const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
    bool has_mcps)
{
    SemanticClassifier::SemanticLabel out = SemanticClassifier::SemanticLabel::Cosmic;
    if (!has_mcps || !mcp_vector.isValid()) return out;

    std::vector<art::Ptr<simb::MCParticle>> parts;
    std::vector<anab::BackTrackerHitMatchingData const*> bkd;
    mcp_bkth_assoc.get(matched_hit.key(), parts, bkd);
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

    if (max_ide_fraction <= 0.5f) return SemanticClassifier::SemanticLabel::Ambiguous;
    return semantic_label_vector[idx];
}

inline void ImageProduction::paint_hit_energy(
    recob::Hit const& hit,
    recob::Wire const& wire,
    geo::PlaneID const& planeID,
    geo::View_t view,
    TVector3 const& wire_center,
    size_t view_idx,
    std::vector<cal::RangeRef> const& roi_ranges,
    BuildContext const& ctx)
{
    const int tick_start = hit.StartTick();
    const int tick_end   = hit.EndTick();
    const double tick_c = static_cast<double>(hit.PeakTime());
    auto p_nom = cal::nominal_point_from_tick(ctx.detp_api, ctx.detprop_data, planeID, tick_c, wire_center);
    auto p = cal::correct_spacepoint(ctx.sce, p_nom, /*always_apply=*/true);
    const unsigned plane = planeID.Plane;
    const double yz = cal::yz_gain_corr(ctx.tpcCalib, plane, p, /*always_apply=*/true);
    double pitch_cm = 1.0;
    if (ctx.geo) pitch_cm = std::max(1e-6, ctx.geo->Plane(planeID).WirePitch());
    const double E_loc = cal::local_efield_kV_cm(ctx.sce, ctx.detprop_data, p, /*always_apply=*/true);
    const double T0_ns = ctx.detprop_data ? (ctx.detprop_data->SamplingRate() * ctx.T0_ticks) : 0.0;
    const double phi = 0.0;
    double dEdx_MeV_cm = 0.0;
    if (ctx.calo_alg && ctx.detprop_data && ctx.clocks && pitch_cm > 0.0) {
        const double dQdx_adc_cm = (hit.Integral() / pitch_cm) * yz;
        dEdx_MeV_cm = ctx.calo_alg->dEdx_AREA(
            *ctx.clocks, *ctx.detprop_data, dQdx_adc_cm, tick_c, plane, T0_ns, E_loc, phi);
    }
    const double E_hit_MeV = dEdx_MeV_cm * pitch_cm;
    const double sumw = cal::sum_adc_weights_in_window(roi_ranges, tick_start, tick_end, ctx.adc_image_threshold);
    const double wire_coord = cal::wire_column_coordinate(view, wire_center);

    for (auto const& rr : roi_ranges) {
        const int rbeg = rr.begin;
        const int rend = rbeg + static_cast<int>(rr.adcs->size());
        const int s = std::max(tick_start, rbeg);
        const int e = std::min(tick_end,   rend);
        if (s >= e) continue;

        for (int t = s; t < e; ++t) {
            const float a = (*rr.adcs)[t - rbeg];
            if (a <= ctx.adc_image_threshold) continue;

            
            auto row = cal::tick_to_corrected_row(ctx.detp_api, ctx.detprop_data, ctx.sce,
                                                  planeID, t, wire_center, /*always_apply_sce=*/true,
                                                  ctx.properties[view_idx]);
            auto col = ctx.properties[view_idx].col(wire_coord);
            if (!row || !col) continue;

            const double w = (sumw > 0.0) ? (static_cast<double>(a) / sumw) : 0.0;
            const float E_pix = static_cast<float>(E_hit_MeV * w);

            
            const float prev = ctx.detector_images[view_idx].get(*row, *col);
            ctx.detector_images[view_idx].set(*row, *col, prev + E_pix);

            
            if (ctx.has_mcps) {




            }
        }
    }
}

inline void ImageProduction::fillDetectorImage(
    const recob::Wire &wire,
    size_t wire_idx,
    const std::set<art::Ptr<recob::Hit>> &hit_set,
    BuildContext const& ctx)
{
    const unsigned ch = wire.Channel();
    if (ctx.chan_status && ctx.chan_status->IsBad(ch)) return;

    auto wire_ids = ctx.geo->ChannelToWire(ch);
    if (wire_ids.empty()) return;

    const auto planeID = wire_ids.front().planeID();
    const unsigned plane = planeID.Plane;
    const geo::View_t view = ctx.geo->View(planeID);
    const size_t view_idx = static_cast<size_t>(view);

    const geo::WireGeo* wire_geo = ctx.geo->WirePtr(wire_ids.front());
    const TVector3 wire_center = wire_geo->GetCenter();

    
    auto hits_for_wire = ctx.wire_hit_assoc.at(wire_idx);
    std::vector<art::Ptr<recob::Hit>> hits_filtered;
    hits_filtered.reserve(hits_for_wire.size());
    for (auto const& ph : hits_for_wire) if (hit_set.count(ph)) hits_filtered.push_back(ph);
    if (hits_filtered.empty()) return;

    
    std::vector<cal::RangeRef> ranges;
    ranges.reserve(wire.SignalROI().get_ranges().size());
    for (auto const& r : wire.SignalROI().get_ranges())
        ranges.push_back(cal::RangeRef{ r.begin_index(), &r.data() });

    
    for (auto const& ph : hits_filtered) {
        if (static_cast<int>(plane) != kCollectionPlane)
            continue;

        
        paint_hit_energy(*ph, wire, planeID, view, wire_center, view_idx, ranges, ctx);

        
        if (ctx.has_mcps) {
            
            const int tick_start = ph->StartTick();
            const int tick_end   = ph->EndTick();
            const double wire_coord = cal::wire_column_coordinate(view, wire_center);

            for (auto const& rr : ranges) {
                const int rbeg = rr.begin;
                const int rend = rbeg + static_cast<int>(rr.adcs->size());
                const int s = std::max(tick_start, rbeg);
                const int e = std::min(tick_end,   rend);
                if (s >= e) continue;

                for (int t = s; t < e; ++t) {
                    const float a = (*rr.adcs)[t - rbeg];
                    if (a <= ctx.adc_image_threshold) continue;

                    auto row = cal::tick_to_corrected_row(ctx.detp_api, ctx.detprop_data, ctx.sce,
                                                          planeID, t, wire_center, /*always_apply_sce=*/true,
                                                          ctx.properties[view_idx]);
                    auto col = ctx.properties[view_idx].col(wire_coord);
                    if (!row || !col) continue;

                    auto sem = labelSemanticPixels(ph, ctx.mcp_bkth_assoc, ctx.mcp_vector,
                                                   ctx.trackid_to_index, ctx.semantic_label_vector,
                                                   ctx.has_mcps);
                    ctx.semantic_images[view_idx].set(*row, *col, static_cast<int>(sem), false);
                }
            }
        }
    }
}

inline void ImageProduction::build(
    const art::Event &event,
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const std::vector<ImageProperties> &properties,
    std::vector<Image<float>> &detector_images,
    std::vector<Image<int>> &semantic_images,
    const detinfo::DetectorProperties *detp,
    std::optional<CalibrationContext> const& cal) const
{
    detector_images.clear();
    semantic_images.clear();

    for (auto const& p : properties) {
        Image<float> det(p); det.clear(0.0f); detector_images.push_back(std::move(det));
        Image<int>   sem(p); sem.clear(static_cast<int>(SemanticClassifier::SemanticLabel::Empty));
        semantic_images.push_back(std::move(sem));
    }

    auto wires = event.getValidHandle<std::vector<recob::Wire>>(opts_.producers.wire);
    auto all_hits = event.getValidHandle<std::vector<recob::Hit>>(opts_.producers.hit);

    art::Handle<std::vector<simb::MCParticle>> mcps;
    bool has_mcps = event.getByLabel(opts_.producers.mcp, mcps);

    art::FindManyP<recob::Hit> wire_hit_assoc(wires, event, opts_.producers.hit);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>
        mcp_bkth_assoc(all_hits, event, opts_.producers.bkt);

    std::vector<SemanticClassifier::SemanticLabel> sem_labels;
    if (has_mcps && mcps.isValid() && opts_.semantic)
        sem_labels = opts_.semantic->classifyParticles(event);

    std::map<int, size_t> trkid_to_idx;
    if (has_mcps && mcps.isValid())
        for (size_t i = 0; i < mcps->size(); ++i) trkid_to_idx[mcps->at(i).TrackId()] = i;

    
    std::set<art::Ptr<recob::Hit>> hit_set(hits.begin(), hits.end());
    BuildContext ctx{
        properties,
        detector_images,
        semantic_images,
        wire_hit_assoc,
        hit_set,
        mcp_bkth_assoc,
        mcps,
        trkid_to_idx,
        sem_labels,
        has_mcps,
        geo_,
        detp,
        kAdcThreshold,
        (cal && cal->enabled()) ? cal->calo      : nullptr,
        (cal && cal->enabled()) ? cal->clocks    : nullptr,
        (cal && cal->enabled()) ? cal->detprop   : nullptr,
        (cal && cal->enabled()) ? cal->tpcCalib  : nullptr,
        (cal && cal->enabled()) ? cal->sce       : nullptr,
        (cal && cal->enabled()) ? cal->chanStatus: nullptr,
        (cal && cal->enabled()) ? cal->T0_ticks  : 0.0
    };

    
    for (size_t wi = 0; wi < wires->size(); ++wi)
        fillDetectorImage(wires->at(wi), wi, hit_set, ctx);
}

} 

#endif
