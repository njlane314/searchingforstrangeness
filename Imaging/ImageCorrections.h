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

#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardataobj/RecoBase/Hit.h"

// Space charge service (brings in spacecharge::SpaceCharge)
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <TVector3.h>

namespace image {
namespace cal {

// Geometry corrections -------------------------------------------------------

inline geo::Point_t correctedPointFromTick(detinfo::DetectorProperties const* detprop,
                                           spacecharge::SpaceCharge const* sce,
                                           geo::PlaneID const& planeID,
                                           TVector3 const& wire_center,
                                           int tick)
{
    const double x_nom = detprop->ConvertTicksToX(static_cast<double>(tick), planeID);
    geo::Point_t p{x_nom, wire_center.Y(), wire_center.Z()};
    if (sce && sce->EnableCalSpatialSCE()) {
        // v08-era interface takes a single point argument
        auto off = sce->GetCalPosOffsets(p);
        p = geo::Point_t{ p.X() - off.X(), p.Y() + off.Y(), p.Z() + off.Z() };
    }
    return p;
}

struct GeometryResult {
    geo::Point_t p_corr;
    double wire_coord{0.0};
    std::optional<size_t> col;

    detinfo::DetectorProperties const* detprop{nullptr};
    spacecharge::SpaceCharge const* sce{nullptr};
    geo::PlaneID planeID;
    TVector3 wire_center;
    ImageProperties const* prop{nullptr};

    // Convenience helper using cached context.
    inline std::optional<size_t> row(int tick) const
    {
        auto p = correctedPointFromTick(detprop, sce, planeID, wire_center, tick);
        return prop->row(p.X());
    }
};

inline GeometryResult applyGeometry(detinfo::DetectorProperties const* detprop,
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
    double yz_corr{1.0}; // kept for completeness; set to 1.0 w/o DB correction
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
    // No DB correction dependency in this build; leave yz_corr at 1.0
    out.E_loc_kV_cm = detprop ? detprop->Efield() : 0.0;
    if (sce && sce->EnableCalEfieldSCE()) {
        auto fo = sce->GetCalEfieldOffsets(p_corr); // one-arg signature in v08
        if (detprop) {
            const double ex = 1.0 + fo.X();
            out.E_loc_kV_cm *= std::sqrt(ex * ex + fo.Y() * fo.Y() + fo.Z() * fo.Z());
        }
    }
    if (calo_alg && detprop && pitch_cm > 0.0) {
        const double T0_ns = detprop->SamplingRate() * T0_ticks;
        // Use legacy CalorimetryAlg signature for v08
        out.dEdx_MeV_cm = calo_alg->dEdx_AREA(hit, pitch_cm, T0_ns);
        out.E_hit_MeV = out.dEdx_MeV_cm * pitch_cm;
    }
    return out;
}

}
}

#endif
