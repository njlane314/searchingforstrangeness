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

struct PixelImageOptions {
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
    } producers;
    image::SemanticClassifier* semantic{nullptr};
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

struct GeometryResult {
    geo::Point_t p_corr;
    double wire_coord{0.0};
    std::optional<size_t> col;

    detinfo::DetectorPropertiesData const* detprop{nullptr};
    spacecharge::SpaceCharge const* sce{nullptr};
    geo::PlaneID planeID;
    TVector3 wire_center;
    ImageProperties const* prop{nullptr};

    inline std::optional<size_t> row(int tick) const
    {
        double x_nom = detprop->ConvertTicksToX(static_cast<double>(tick), planeID);
        geo::Point_t p{x_nom, wire_center.Y(), wire_center.Z()};
        if (sce && sce->EnableCalSpatialSCE()) {
            auto off = sce->GetCalPosOffsets(p, 0);
            p = geo::Point_t{ p.X() - off.X(), p.Y() + off.Y(), p.Z() + off.Z() };
        }
        return prop->row(p.X());
    }
};

inline GeometryResult applyGeometry(detinfo::DetectorPropertiesData const* detprop,
                                    spacecharge::SpaceCharge const* sce,
                                    geo::PlaneID const& planeID,
                                    int tick_center,
                                    TVector3 const& wire_center,
                                    geo::View_t view,
                                    ImageProperties const& prop)
{
    double x_nom = detprop->ConvertTicksToX(static_cast<double>(tick_center), planeID);
    geo::Point_t p_corr{x_nom, wire_center.Y(), wire_center.Z()};
    if (sce && sce->EnableCalSpatialSCE()) {
        auto off = sce->GetCalPosOffsets(p_corr, 0);
        p_corr = geo::Point_t{ p_corr.X() - off.X(), p_corr.Y() + off.Y(), p_corr.Z() + off.Z() };
    }

    constexpr double plus60  =  1.04719758034;
    constexpr double minus60 = -1.04719758034;
    double wire_coord = (view == geo::kW) ? wire_center.Z()
                      : (view == geo::kU)
                        ? (wire_center.Z() * std::cos(plus60) - wire_center.Y() * std::sin(plus60))
                        : (wire_center.Z() * std::cos(minus60) - wire_center.Y() * std::sin(minus60));

    GeometryResult out;
    out.p_corr = p_corr;
    out.wire_coord = wire_coord;
    out.col = prop.col(wire_coord);
    out.detprop = detprop;
    out.sce = sce;
    out.planeID = planeID;
    out.wire_center = wire_center;
    out.prop = &prop;
    return out;
}

struct CaloResult {
    double E_hit_MeV{0.0};
    double dEdx_MeV_cm{0.0};
    double yz_corr{1.0};
    double E_loc_kV_cm{0.0};
};

inline CaloResult applyCalorimetry(recob::Hit const& hit,
                                   unsigned plane,
                                   geo::Point_t const& p_corr,
                                   double pitch_cm,
                                   calo::CalorimetryAlg* calo_alg,
                                   detinfo::DetectorClocksData const* clocks,
                                   detinfo::DetectorPropertiesData const* detprop,
                                   lariov::TPCEnergyCalib const* tpcCalib,
                                   spacecharge::SpaceCharge const* sce,
                                   double T0_ticks)
{
    CaloResult out;
    out.yz_corr = (tpcCalib ? tpcCalib->YZdqdxCorrection(plane, p_corr.Y(), p_corr.Z()) : 1.0);
    out.E_loc_kV_cm = detprop->Efield();
    if (sce && sce->EnableCalEfieldSCE()) {
        auto fo = sce->GetCalEfieldOffsets(p_corr, 0);
        out.E_loc_kV_cm *= std::hypot(1.0 + fo.X(), fo.Y(), fo.Z());
    }
    if (calo_alg && clocks && detprop && pitch_cm > 0.0) {
        const double T0_ns = detprop->SamplingRate() * T0_ticks;
        const double dQdx_adc_cm = (hit.Integral() / pitch_cm) * out.yz_corr;
        out.dEdx_MeV_cm = calo_alg->dEdx_AREA(*clocks, *detprop,
                                              dQdx_adc_cm,
                                              static_cast<double>(hit.PeakTime()),
                                              plane, T0_ns, out.E_loc_kV_cm,
                                              0.0);
        out.E_hit_MeV = out.dEdx_MeV_cm * pitch_cm;
    }
    return out;
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

inline std::vector<float> gaussianKernel1D(float sigma_px)
{
    if (!(sigma_px > 0.f)) return {1.f};
    const int R = std::max(1, static_cast<int>(std::ceil(3.0f * sigma_px)));
    std::vector<float> k(2 * R + 1);
    double s = 0.0;
    for (int i = -R; i <= R; ++i) {
        const double v = std::exp(-(static_cast<double>(i) * static_cast<double>(i)) / (2.0 * sigma_px * sigma_px));
        k[i + R] = static_cast<float>(v);
        s += v;
    }
    const float invs = static_cast<float>(1.0 / s);
    for (auto& v : k) v *= invs;
    return k;
}

inline void convolveSeparableUnitSum(Image<float>& img,
                                     ImageProperties const& prop,
                                     std::vector<float> const& k_row,
                                     std::vector<float> const& k_col)
{
    const int Rr = static_cast<int>(k_row.size() / 2);
    const int Rc = static_cast<int>(k_col.size() / 2);

    const size_t H = prop.rows();
    const size_t W = prop.cols();

    Image<float> tmp(prop);
    tmp.clear(0.0f);

    for (size_t r = 0; r < H; ++r) {
        for (size_t c = 0; c < W; ++c) {
            double acc = 0.0;
            for (int dr = -Rr; dr <= Rr; ++dr) {
                const int rr = std::min<int>(std::max<int>(0, static_cast<int>(r) + dr), static_cast<int>(H) - 1);
                acc += k_row[dr + Rr] * img.get(rr, c);
            }
            tmp.set(r, c, static_cast<float>(acc));
        }
    }

    Image<float> out(prop);
    out.clear(0.0f);
    for (size_t r = 0; r < H; ++r) {
        for (size_t c = 0; c < W; ++c) {
            double acc = 0.0;
            for (int dc = -Rc; dc <= Rc; ++dc) {
                const int cc = std::min<int>(std::max<int>(0, static_cast<int>(c) + dc), static_cast<int>(W) - 1);
                acc += k_col[dc + Rc] * tmp.get(r, cc);
            }
            out.set(r, c, static_cast<float>(acc));
        }
    }
    img = std::move(out);
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
               const detinfo::DetectorProperties *,
               std::optional<CalibrationContext> const& cal = std::nullopt) const;

private:
    static constexpr int   kCollectionPlane = 2;
    static constexpr float kAdcThreshold    = 4.0f;

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

        float adc_image_threshold;

        calo::CalorimetryAlg* calo_alg;
        detinfo::DetectorClocksData const* clocks;
        detinfo::DetectorPropertiesData const* detprop_data;
        lariov::TPCEnergyCalib const* tpcCalib;
        spacecharge::SpaceCharge const* sce;
        lariov::ChannelStatusProvider const* chanStatus;

        double T0_ticks{0.0};
    };

    static SemanticClassifier::SemanticLabel labelSemanticPixels(
        const art::Ptr<recob::Hit> &matched_hit,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> &mcp_bkth_assoc,
        const art::Handle<std::vector<simb::MCParticle>> &mcp_vector,
        const std::map<int, size_t> &trackid_to_index,
        const std::vector<SemanticClassifier::SemanticLabel> &semantic_label_vector,
        bool has_mcps);

    static void rasterizeHitEnergy(
        recob::Hit const& hit,
        recob::Wire const&,
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

inline void ImageProduction::rasterizeHitEnergy(
    recob::Hit const& hit,
    recob::Wire const&,
    geo::PlaneID const& planeID,
    geo::View_t view,
    TVector3 const& wire_center,
    size_t view_idx,
    std::vector<cal::RangeRef> const& roi_ranges,
    BuildContext const& ctx)
{
    const unsigned plane = planeID.Plane;
    const int tick_c = static_cast<int>(hit.PeakTime());

    auto geoRes = cal::applyGeometry(ctx.detprop_data, ctx.sce, planeID,
                                     tick_c, wire_center, view, ctx.properties[view_idx]);
    if (!geoRes.col) return;

    double pitch_cm = 1.0;
    if (ctx.geo) pitch_cm = std::max(1e-6, ctx.geo->Plane(planeID).WirePitch());

    auto calRes = cal::applyCalorimetry(hit, plane, geoRes.p_corr, pitch_cm,
                                        ctx.calo_alg, ctx.clocks, ctx.detprop_data,
                                        ctx.tpcCalib, ctx.sce, ctx.T0_ticks);

    const int tick_start = hit.StartTick();
    const int tick_end   = hit.EndTick();
    const double sumw = cal::sum_adc_weights_in_window(roi_ranges, tick_start, tick_end, ctx.adc_image_threshold);

    for (auto const& rr : roi_ranges) {
        const int rbeg = rr.begin;
        const int rend = rbeg + static_cast<int>(rr.adcs->size());
        const int s = std::max(tick_start, rbeg);
        const int e = std::min(tick_end,   rend);
        if (s >= e) continue;

        for (int t = s; t < e; ++t) {
            const float a = (*rr.adcs)[t - rbeg];
            if (a <= ctx.adc_image_threshold) continue;

            auto row = geoRes.row(t);
            if (!row) continue;

            const double w = (sumw > 0.0) ? (static_cast<double>(a) / sumw) : 0.0;
            const float E_pix = static_cast<float>(calRes.E_hit_MeV * w);

            const float prev = ctx.detector_images[view_idx].get(*row, *geoRes.col);
            ctx.detector_images[view_idx].set(*row, *geoRes.col, prev + E_pix);
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
    if (ctx.chanStatus && ctx.chanStatus->IsBad(ch)) return;

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

        rasterizeHitEnergy(*ph, wire, planeID, view, wire_center, view_idx, ranges, ctx);

        if (ctx.has_mcps) {
            const int tick_c = static_cast<int>(ph->PeakTime());
            auto geoRes = cal::applyGeometry(ctx.detprop_data, ctx.sce, planeID,
                                             tick_c, wire_center, view, ctx.properties[view_idx]);
            if (!geoRes.col) continue;

            const int tick_start = ph->StartTick();
            const int tick_end   = ph->EndTick();

            for (auto const& rr : ranges) {
                const int rbeg = rr.begin;
                const int rend = rbeg + static_cast<int>(rr.adcs->size());
                const int s = std::max(tick_start, rbeg);
                const int e = std::min(tick_end,   rend);
                if (s >= e) continue;

                for (int t = s; t < e; ++t) {
                    const float a = (*rr.adcs)[t - rbeg];
                    if (a <= ctx.adc_image_threshold) continue;

                    auto row = geoRes.row(t);
                    if (!row) continue;

                    auto sem = labelSemanticPixels(ph, ctx.mcp_bkth_assoc, ctx.mcp_vector,
                                                   ctx.trackid_to_index, ctx.semantic_label_vector,
                                                   ctx.has_mcps);
                    ctx.semantic_images[view_idx].set(*row, *geoRes.col, static_cast<int>(sem), false);
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
    const detinfo::DetectorProperties *,
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

    {
        constexpr float kSigmaPx = 1.0f;
        const auto k = cal::gaussianKernel1D(kSigmaPx);
        for (size_t v = 0; v < detector_images.size(); ++v) {
            cal::convolveSeparableUnitSum(
                detector_images[v], properties[v],
                k, k
            );
        }
    }
}

}

#endif
