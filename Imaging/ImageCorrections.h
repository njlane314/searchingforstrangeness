#ifndef IMAGECORRECTIONS_H
#define IMAGECORRECTIONS_H

#include <algorithm>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include "Imaging/Image.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "lardata/DetectorInfo/DetectorClocksData.h"
#include "lardata/DetectorInfo/DetectorPropertiesData.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larevt/CalibrationServices/TPCEnergyCalibService.h"

#include <TVector3.h>

namespace image {
namespace cal {

inline geo::Point_t correctedPointFromTick(detinfo::DetectorPropertiesData const* detprop,
                                           spacecharge::SpaceCharge const* sce,
                                           geo::PlaneID const& planeID,
                                           TVector3 const& wire_center,
                                           int tick)
{
    const double x_nom = detprop->ConvertTicksToX(static_cast<double>(tick), planeID);
    geo::Point_t p{x_nom, wire_center.Y(), wire_center.Z()};
    if (sce && sce->EnableCalSpatialSCE()) {
        const unsigned tpcID = planeID.TPC;
        auto off = sce->GetCalPosOffsets(p, tpcID);
        p = geo::Point_t{ p.X() - off.X(), p.Y() + off.Y(), p.Z() + off.Z() };
    }
    return p;
}

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
        auto p = correctedPointFromTick(detprop, sce, planeID, wire_center, tick);
        return prop->row(p.X());
    }
};

inline GeometryResult applyGeometry(detinfo::DetectorPropertiesData const* detprop,
                                    spacecharge::SpaceCharge const* sce,
                                    geo::PlaneID const& planeID,
                                    int tick_center,
                                    TVector3 const& wire_center,
                                    ImageProperties const& prop)
{
    geo::Point_t p_corr = correctedPointFromTick(detprop, sce, planeID, wire_center, tick_center);

    constexpr double plus60  =  1.04719758034;
    constexpr double minus60 = -1.04719758034;
    const double y = p_corr.Y();
    const double z = p_corr.Z();
    const geo::View_t view = prop.view();
    double wire_coord = (view == geo::kW) ? z
                      : (view == geo::kU)
                        ? (z * std::cos(plus60) - y * std::sin(plus60))
                        : (z * std::cos(minus60) - y * std::sin(minus60));

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

    const size_t H = prop.height();
    const size_t W = prop.width();

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
}

#endif
