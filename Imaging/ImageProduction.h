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

#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace detinfo {
  using DetectorPropertiesData = detinfo::DetectorProperties;
}

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationServices/TPCEnergyCalibService.h"
#include "larevt/CalibrationServices/UBElectronLifetimeService.h"
#include "larevt/CalibrationServices/ChannelStatusService.h"

#include <TVector3.h>

namespace image {

// ========================== Options & Context ==========================

struct PixelImageOptions {
    bool is_data{false};
    struct Producers {
        art::InputTag wire;
        art::InputTag hit;
        art::InputTag mcp;
        art::InputTag bkt;
        art::InputTag bad; // optional per-event bad-channel product (vector<int>)
    } producers;

    float adc_threshold{4.0f};
    const std::set<unsigned int>* bad_channels{nullptr}; // optional static mask
    image::SemanticClassifier* semantic{nullptr};

    // Calorimetry & corrections toggles (match BlipRecoAlg)
    int  calo_plane{2};           // collection by default
    bool apply_yz{true};
    bool apply_lifetime{true};
    bool apply_sce_spatial{true};
    bool apply_sce_efield{true};
    bool veto_noisy{false};

    // Modified Box parameters and assumed dE/dx (same as blips)
    float modbox_A{0.93f};
    float modbox_B{0.212f};
    float assumed_dEdx{2.8f}; // MeV/cm

    enum class ValueMode { ADC, dEdx, Energy };
    ValueMode value_mode{ValueMode::Energy};
};

struct CalibrationContext {
    calo::CalorimetryAlg* calo{nullptr};
    detinfo::DetectorPropertiesData const* detprop{nullptr};

    lariov::TPCEnergyCalib const* tpcCalib{nullptr};
    lariov::UBElectronLifetime const* lifetime{nullptr};
    spacecharge::SpaceCharge const* sce{nullptr};
    lariov::ChannelStatusProvider const* chanStatus{nullptr};

    double T0_ticks{0.0};

    bool enabled() const { return calo && detprop; }
};

// ============================== Helpers ===============================

namespace cal {

// --- geometry mapping ---

inline geo::Point_t nominal_point_from_tick(detinfo::DetectorProperties const* detp_api,
                                            detinfo::DetectorPropertiesData const* detprop,
                                            geo::PlaneID const& planeID,
                                            double tick,
                                            TVector3 const& wire_center)
{
    (void)detp_api; // kept for signature parity with your original interfaces
    double x_nom = detprop->ConvertTicksToX(tick, planeID);
    return geo::Point_t{x_nom, wire_center.Y(), wire_center.Z()};
}

inline geo::Point_t correct_spacepoint(spacecharge::SpaceCharge const* sce,
                                       geo::Point_t p_in,
                                       bool apply_sce_spatial)
{
    if (sce && apply_sce_spatial && sce->EnableCalSpatialSCE()) {
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
                                                   bool apply_sce_spatial,
                                                   ImageProperties const& prop)
{
    auto p_nom = nominal_point_from_tick(detp_api, detprop, planeID, static_cast<double>(tick), wire_center);
    auto p = correct_spacepoint(sce, p_nom, apply_sce_spatial);
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

// --- calibrations and physics ---

inline double yz_gain_corr(lariov::TPCEnergyCalib const* tpcCalib,
                           unsigned plane,
                           geo::Point_t const& p,
                           bool apply_yz)
{
    if (tpcCalib && apply_yz) return tpcCalib->YZdqdxCorrection(plane, p.Y(), p.Z());
    return 1.0;
}

inline double lifetime_corr(lariov::UBElectronLifetime const* lifetime,
                            detinfo::DetectorPropertiesData const* detprop,
                            double T0_ticks,
                            double tick,
                            bool apply_lifetime)
{
    if (!lifetime || !detprop || !apply_lifetime) return 1.0;
    const double tick_us = detprop->SamplingRate() * 1e-3; // ns -> us
    const double t_ms    = (tick - T0_ticks) * tick_us * 1e-3;
    const double tau_ms  = lifetime->Lifetime();
    if (tau_ms <= 0.0) return 1.0;
    return std::exp(t_ms / tau_ms);
}

inline double electrons_from_adc_area(calo::CalorimetryAlg* calo,
                                      recob::Hit const& hit,
                                      unsigned plane)
{
    if (!calo) return 0.0;
    return static_cast<double>(calo->ElectronsFromADCArea(hit.Integral(), plane));
}

inline double local_efield_kV_cm(spacecharge::SpaceCharge const* sce,
                                 detinfo::DetectorPropertiesData const* detprop,
                                 geo::Point_t const& p,
                                 bool apply_sce_efield)
{
    double E = detprop->Efield();
    if (sce && apply_sce_efield && sce->EnableCalEfieldSCE()) {
        auto fo = sce->GetCalEfieldOffsets(p, 0);
        E *= std::hypot(1.0 + fo.X(), fo.Y(), fo.Z());
    }
    return E;
}

inline double modbox_recomb(double dEdx_MeV_cm,
                            double Efield_kV_cm,
                            double LArDensity_g_cm3,
                            double A, double B)
{
    // same as your BlipRecoAlg
    const double Xi = (B * dEdx_MeV_cm) / (Efield_kV_cm * LArDensity_g_cm3);
    return std::log(A + Xi) / Xi; // fraction of charge that survives
}

inline double energy_from_electrons(double Q_electrons,
                                    double recomb,
                                    double Wion_MeV_per_e)
{
    if (recomb <= 0.0) return 0.0;
    return Q_electrons * (1.0 / recomb) * Wion_MeV_per_e;
}

// --- ROI utilities ---

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

} // namespace cal

// ============================ Builder ============================

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
        const detinfo::DetectorProperties *detp_api;

        float adc_image_threshold;

        // Combined bad channels (static + per-event)
        std::set<unsigned int> merged_bad;

        // Calibration & services
        calo::CalorimetryAlg* calo_alg;
        detinfo::DetectorPropertiesData const* detprop_data;
        lariov::TPCEnergyCalib const* tpcCalib;
        lariov::UBElectronLifetime const* lifetime;
        spacecharge::SpaceCharge const* sce;
        lariov::ChannelStatusProvider const* chan_status;

        // Knobs
        double T0_ticks{0.0};
        int  calo_plane{2};
        bool apply_yz{true};
        bool apply_lifetime{true};
        bool apply_sce_spatial{true};
        bool apply_sce_efield{true};
        bool veto_noisy{false};

        float modbox_A{0.93f};
        float modbox_B{0.212f};
        float assumed_dEdx{2.8f};

        // Material + constants
        double LArDensity{1.396};         // overwritten from detprop when present
        double WionMeVPerElectron{2.36e-5};

        PixelImageOptions::ValueMode value_mode{PixelImageOptions::ValueMode::Energy};
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

// ----------------- Implementation -----------------

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

    // Representative tick (peak) for corrections
    const double tick_c = static_cast<double>(hit.PeakTime());

    // Build nominal point at tick_c, then SCE spatial correction
    auto p_nom = cal::nominal_point_from_tick(ctx.detp_api, ctx.detprop_data, planeID, tick_c, wire_center);
    auto p = cal::correct_spacepoint(ctx.sce, p_nom, ctx.apply_sce_spatial);

    // Per-hit calibrations
    const unsigned plane = planeID.Plane;
    const double yz    = cal::yz_gain_corr(ctx.tpcCalib, plane, p, ctx.apply_yz);
    const double life  = cal::lifetime_corr(ctx.lifetime, ctx.detprop_data, ctx.T0_ticks, tick_c, ctx.apply_lifetime);
    double Q_e         = cal::electrons_from_adc_area(ctx.calo_alg, hit, plane) * yz * life;

    const double E_loc = cal::local_efield_kV_cm(ctx.sce, ctx.detprop_data, p, ctx.apply_sce_efield);
    const double recomb = cal::modbox_recomb(static_cast<double>(ctx.assumed_dEdx),
                                             E_loc, ctx.LArDensity,
                                             static_cast<double>(ctx.modbox_A),
                                             static_cast<double>(ctx.modbox_B));

    const double E_hit_MeV = cal::energy_from_electrons(Q_e, recomb, ctx.WionMeVPerElectron);

    // Sum ADC weights across this hit window
    const double sumw = cal::sum_adc_weights_in_window(roi_ranges, tick_start, tick_end, ctx.adc_image_threshold);

    // Column coordinate for this wire (fixed)
    const double wire_coord = cal::wire_column_coordinate(view, wire_center);

    // Paint energy across ticks ∝ ADC
    for (auto const& rr : roi_ranges) {
        const int rbeg = rr.begin;
        const int rend = rbeg + static_cast<int>(rr.adcs->size());
        const int s = std::max(tick_start, rbeg);
        const int e = std::min(tick_end,   rend);
        if (s >= e) continue;

        for (int t = s; t < e; ++t) {
            const float a = (*rr.adcs)[t - rbeg];
            if (a <= ctx.adc_image_threshold) continue;

            // Row position from per-tick (SCE-corrected) X
            auto row = cal::tick_to_corrected_row(ctx.detp_api, ctx.detprop_data, ctx.sce,
                                                  planeID, t, wire_center, ctx.apply_sce_spatial,
                                                  ctx.properties[view_idx]);
            auto col = ctx.properties[view_idx].col(wire_coord);
            if (!row || !col) continue;

            const double w = (sumw > 0.0) ? (static_cast<double>(a) / sumw) : 0.0;
            const float E_pix = static_cast<float>(E_hit_MeV * w);

            // Accumulate; if your Image<T> lacks get(r,c), replace this with your accessor
            const float prev = ctx.detector_images[view_idx].get(*row, *col);
            ctx.detector_images[view_idx].set(*row, *col, prev + E_pix);

            // Optional semantic label (do not overwrite)
            if (!ctx.is_data) {
                // We need a Ptr<recob::Hit> for semantics; this function receives a const&.
                // Callers invoke paint_hit_energy from within a loop over Ptr<recob::Hit>
                // if they also want semantics. When called from fillDetectorImage (below),
                // semantics are set there where Ptr<> is available.
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

    // Static + per-event channel veto
    if (ctx.merged_bad.count(ch)) return;
    if (ctx.chan_status) {
        if (ctx.chan_status->IsBad(ch)) return;
        if (ctx.veto_noisy && ctx.chan_status->IsNoisy(ch)) return;
    }

    auto wire_ids = ctx.geo->ChannelToWire(ch);
    if (wire_ids.empty()) return;

    const auto planeID = wire_ids.front().planeID();
    const unsigned plane = planeID.Plane;
    const geo::View_t view = ctx.geo->View(planeID);
    const size_t view_idx = static_cast<size_t>(view);

    const geo::WireGeo* wire_geo = ctx.geo->WirePtr(wire_ids.front());
    const TVector3 wire_center = wire_geo->GetCenter();

    // Gather hits for this wire & filter to the requested subset
    auto hits_for_wire = ctx.wire_hit_assoc.at(wire_idx);
    std::vector<art::Ptr<recob::Hit>> hits_filtered;
    hits_filtered.reserve(hits_for_wire.size());
    for (auto const& ph : hits_for_wire) if (hit_set.count(ph)) hits_filtered.push_back(ph);
    if (hits_filtered.empty()) return;

    // Cache ROI ranges for this wire
    std::vector<cal::RangeRef> ranges;
    ranges.reserve(wire.SignalROI().get_ranges().size());
    for (auto const& r : wire.SignalROI().get_ranges())
        ranges.push_back(cal::RangeRef{ r.begin_index(), &r.data() });

    // If writing energy (MeV per pixel), mirror blips and restrict to calo plane
    const bool restrict_to_calo_plane = (ctx.value_mode == PixelImageOptions::ValueMode::Energy);

    for (auto const& ph : hits_filtered) {
        if (restrict_to_calo_plane && static_cast<int>(plane) != ctx.calo_plane)
            continue;

        // Paint energy for this hit
        paint_hit_energy(*ph, wire, planeID, view, wire_center, view_idx, ranges, ctx);

        // Optional semantic label pass (per pixel) while we have Ptr<recob::Hit>
        if (!ctx.is_data) {
            // Re-walk the pixels touched by this hit to set labels (cheap; small windows)
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
                                                          planeID, t, wire_center, ctx.apply_sce_spatial,
                                                          ctx.properties[view_idx]);
                    auto col = ctx.properties[view_idx].col(wire_coord);
                    if (!row || !col) continue;

                    auto sem = labelSemanticPixels(ph, ctx.mcp_bkth_assoc, ctx.mcp_vector,
                                                   ctx.trackid_to_index, ctx.semantic_label_vector,
                                                   ctx.has_mcps);
                    ctx.semantic_images[view_idx].set(*row, *col, static_cast<int>(sem), /*overwrite=*/false);
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
    if (!opts_.is_data && has_mcps && mcps.isValid() && opts_.semantic)
        sem_labels = opts_.semantic->classifyParticles(event);

    std::map<int, size_t> trkid_to_idx;
    if (!opts_.is_data && has_mcps && mcps.isValid())
        for (size_t i = 0; i < mcps->size(); ++i) trkid_to_idx[mcps->at(i).TrackId()] = i;

    // Hits of interest
    std::set<art::Ptr<recob::Hit>> hit_set(hits.begin(), hits.end());

    // Merge static + per-event bad channels
    std::set<unsigned int> merged_bad;
    if (opts_.bad_channels) merged_bad.insert(opts_.bad_channels->begin(), opts_.bad_channels->end());
    if (opts_.producers.bad != art::InputTag{}) {
        art::Handle<std::vector<int>> evtBad;
        if (event.getByLabel(opts_.producers.bad, evtBad) && evtBad.isValid())
            for (int ch : *evtBad) merged_bad.insert(static_cast<unsigned>(ch));
    }

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
        opts_.is_data,
        has_mcps,
        geo_,
        detp,
        opts_.adc_threshold,
        {}, // merged_bad populated below
        (cal && cal->enabled()) ? cal->calo      : nullptr,
        (cal && cal->enabled()) ? cal->detprop   : nullptr,
        (cal && cal->enabled()) ? cal->tpcCalib  : nullptr,
        (cal && cal->enabled()) ? cal->lifetime  : nullptr,
        (cal && cal->enabled()) ? cal->sce       : nullptr,
        (cal && cal->enabled()) ? cal->chanStatus: nullptr,
        (cal && cal->enabled()) ? cal->T0_ticks  : 0.0,
        opts_.calo_plane,
        opts_.apply_yz,
        opts_.apply_lifetime,
        opts_.apply_sce_spatial,
        opts_.apply_sce_efield,
        opts_.veto_noisy,
        opts_.modbox_A,
        opts_.modbox_B,
        opts_.assumed_dEdx,
        1.396,
        2.36e-5,
        opts_.value_mode
    };
    ctx.merged_bad = std::move(merged_bad);
    if (ctx.detprop_data) ctx.LArDensity = ctx.detprop_data->Density();

    // Loop all wires
    for (size_t wi = 0; wi < wires->size(); ++wi)
        fillDetectorImage(wires->at(wi), wi, hit_set, ctx);
}

} // namespace image

#endif
