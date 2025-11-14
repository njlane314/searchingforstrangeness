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

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <TVector3.h>

namespace image {
namespace cal {

// Geometry corrections -------------------------------------------------------

inline geo::Point_t correctedPointFromTick(detinfo::DetectorProperties const* detprop,
                                           spacecharge::SpaceCharge const* sce,
                                           geo::PlaneID const& plane,
                                           TVector3 const& wire_center,
                                           int tick)
{
    const double x_nom = detprop->ConvertTicksToX(static_cast<double>(tick), plane);
    geo::Point_t p{x_nom, wire_center.Y(), wire_center.Z()};

    if (sce && sce->EnableSimSpatialSCE()) {
        auto off = sce->GetPosOffsets(p_meas);  // cm
        return geo::Point_t{ p_meas.X() - off.X(),
                             p_meas.Y() + off.Y(),
                             p_meas.Z() + off.Z() };
    }
    return p;
}

struct GeometryResult {
    geo::Point_t p_corr;
    double wire_coord{0.0};
    std::optional<size_t> col;

    detinfo::DetectorProperties const* detprop{nullptr};
    spacecharge::SpaceCharge const* sce{nullptr};
    geo::PlaneID plane;
    TVector3 wire_center;
    ImageProperties const* prop{nullptr};

    inline std::optional<size_t> row(int tick) const
    {
        auto p = correctedPointFromTick(detprop, sce, plane, wire_center, tick);
        return prop->row(p.X());
    }
};

inline GeometryResult applyGeometry(detinfo::DetectorProperties const* detprop,
                                    spacecharge::SpaceCharge const* sce,
                                    geo::PlaneID const& plane,
                                    int tick_center,
                                    TVector3 const& wire_center,
                                    ImageProperties const& prop)
{
    geo::Point_t p_corr = correctedPointFromTick(detprop, sce, plane, wire_center, tick_center);

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
    out.plane = plane;
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

// Calorimetry corrections ----------------------------------------------------

inline CaloResult applyCalorimetry(recob::Hit const& hit,
                                   unsigned /*plane*/,
                                   geo::Point_t const& p_corr,
                                   double pitch_cm,
                                   calo::CalorimetryAlg* calo_alg,
                                   detinfo::DetectorClocks const* /*clocks*/,
                                   detinfo::DetectorProperties const* detprop,
                                   spacecharge::SpaceCharge const* sce,
                                   double T0_ticks)
{
    CaloResult out;
    out.E_loc_kV_cm = detprop ? detprop->Efield() : 0.0;
    if (sce && sce->EnableCalEfieldSCE()) {
        auto fo = sce->GetCalEfieldOffsets(p_corr); 
        if (detprop) {
            const double ex = 1.0 + fo.X();
            out.E_loc_kV_cm *= std::sqrt(ex * ex + fo.Y() * fo.Y() + fo.Z() * fo.Z());
        }
    }
    if (calo_alg && detprop && pitch_cm > 0.0) {
        const double T0_ns = detprop->SamplingRate() * T0_ticks;
        out.dEdx_MeV_cm = calo_alg->dEdx_AREA(hit, pitch_cm, T0_ns);
        out.E_hit_MeV = out.dEdx_MeV_cm * pitch_cm;
    }
    return out;
}

}
}

#endif
