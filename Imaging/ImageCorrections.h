#ifndef IMAGECORRECTIONS_H
#define IMAGECORRECTIONS_H

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Imaging/Image.h"

#include "Common/PandoraUtilities.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

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
        auto off = sce->GetPosOffsets(p);  // cm
        return geo::Point_t{ p.X() - off.X(),
                             p.Y() + off.Y(),
                             p.Z() + off.Z() };
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

    const geo::View_t view = prop.view();
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
            throw std::runtime_error("applyGeometry: Unsupported geo::View_t: " + std::to_string(view));
    }

    const TVector3 proj = common::ProjectToWireView(static_cast<float>(p_corr.X()),
                                                    static_cast<float>(p_corr.Y()),
                                                    static_cast<float>(p_corr.Z()),
                                                    pandora_view);
    const double wire_coord = static_cast<double>(proj.Z());

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
                                   geo::Point_t const& p_corr,
                                   double pitch_cm,
                                   calo::CalorimetryAlg* calo_alg,
                                   detinfo::DetectorProperties const* detprop,
                                   spacecharge::SpaceCharge const* sce,
                                   double T0_ticks)
{
    CaloResult out;
    out.E_loc_kV_cm = detprop ? detprop->Efield() : 0.0;
    if (sce && sce->EnableSimEfieldSCE()) {
        auto fo = sce->GetEfieldOffsets(p_corr); 
        if (detprop) {
            const double ex = 1.0 + fo.X();
            out.E_loc_kV_cm *= std::sqrt(ex * ex + fo.Y() * fo.Y() + fo.Z() * fo.Z());
        }
    }
    if (calo_alg && detprop && pitch_cm > 0.0) {
        const double T0_ns = 1000.0 * detprop->SamplingRate() * T0_ticks;
        out.dEdx_MeV_cm = calo_alg->dEdx_AREA(hit, pitch_cm, T0_ns);
        out.E_hit_MeV = out.dEdx_MeV_cm * pitch_cm;
    }
    return out;
}

}
}

#endif
