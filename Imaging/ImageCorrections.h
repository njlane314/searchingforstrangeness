#ifndef IMAGECORRECTIONS_H
#define IMAGECORRECTIONS_H

#include <algorithm>
#include <cmath>
#include <optional>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Imaging/Image.h"

#include "Common/PandoraUtilities.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <TVector3.h>
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace image {
namespace cal {

namespace {
constexpr bool kImageCorrectionsCout = false;

template <typename F>
inline void print(F&& fill)
{
    std::ostringstream oss;
    fill(oss);
    const auto msg = oss.str();

    if constexpr (kImageCorrectionsCout) {
        std::cout << "[ImageCorrections] " << msg << std::endl;
    } else {
        mf::LogInfo("ImageCorrections") << msg;
    }
}
}

inline geo::Point_t correctedPointFromTick(detinfo::DetectorProperties const* detprop,
                                           spacecharge::SpaceCharge const* sce,
                                           geo::PlaneID const& plane,
                                           TVector3 const& wire_center,
                                           int tick)
{
    const double x_nom = detprop->ConvertTicksToX(static_cast<double>(tick), plane);
    geo::Point_t p{x_nom, wire_center.Y(), wire_center.Z()};
    geo::Point_t p_corr = p;

    if (sce && sce->EnableSimSpatialSCE()) {
        auto off = sce->GetPosOffsets(p);
        std::cout << "[ImageCorrections] Spatial SCE correction offsets: ("
                  << off.X() << "," << off.Y() << "," << off.Z() << ")" << std::endl;
        p_corr = geo::Point_t{ p.X() - off.X(),
                               p.Y() + off.Y(),
                               p.Z() + off.Z() };
    } else {
        std::cout << "[ImageCorrections] Spatial SCE correction offsets: (0,0,0)" << std::endl;
    }

    print([&](std::ostream& os) {
        os << "correctedPointFromTick:"
           << " tick=" << tick
           << " plane=(" << plane.Cryostat << "," << plane.TPC << "," << plane.Plane << ")"
           << " x_nom=" << x_nom
           << " wire_center=(" << wire_center.X() << "," << wire_center.Y() << "," << wire_center.Z() << ")"
           << " sce=" << (sce && sce->EnableSimSpatialSCE() ? "on" : "off")
           << " p_corr=(" << p_corr.X() << "," << p_corr.Y() << "," << p_corr.Z() << ")";
    });

    return p_corr;
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
        auto r = prop->row(p.X());

        print([&](std::ostream& os) {
            os << "GeometryResult::row:"
               << " tick=" << tick
               << " x_corr=" << p.X()
               << " -> row=" << (r ? std::to_string(*r) : std::string{"n/a"});
        });

        return r;
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

    print([&](std::ostream& os) {
        os << "applyGeometry:"
           << " view=" << static_cast<int>(view)
           << " pandora_view=" << static_cast<int>(pandora_view)
           << " tick_center=" << tick_center
           << " p_corr=(" << p_corr.X() << "," << p_corr.Y() << "," << p_corr.Z() << ")"
           << " wire_coord=" << wire_coord
           << " col=" << (out.col ? std::to_string(*out.col) : std::string{"n/a"})
           << " img_width=" << prop.width()
           << " img_height=" << prop.height();
    });

    return out;
}

struct CaloResult {
    double E_hit_MeV{0.0};
    double dEdx_MeV_cm{0.0};
    double yz_corr{1.0};
    double E_loc_kV_cm{0.0};
};

inline CaloResult applyCalorimetry(recob::Hit const& hit,
                                   geo::Point_t const& p_corr,
                                   double pitch_cm,
                                   calo::CalorimetryAlg* calo_alg,
                                   detinfo::DetectorClocks const* clocks,
                                   detinfo::DetectorProperties const* detprop,
                                   spacecharge::SpaceCharge const* sce,
                                   double T0_ns)
{
    CaloResult out;
    out.E_loc_kV_cm = detprop ? detprop->Efield() : 0.0;
    double nominal_E_kV_cm = out.E_loc_kV_cm;
    const bool use_efield_sce = sce && sce->EnableSimEfieldSCE();
    if (use_efield_sce) {
        auto fo = sce->GetEfieldOffsets(p_corr);
        if (detprop) {
            const double ex = 1.0 + fo.X();
            out.E_loc_kV_cm *= std::sqrt(ex * ex + fo.Y() * fo.Y() + fo.Z() * fo.Z());
        }
        std::cout << "[ImageCorrections] E-field: nominal=" << nominal_E_kV_cm
                  << " kV/cm sce_offsets=(" << fo.X() << "," << fo.Y() << "," << fo.Z() << ")"
                  << " -> corrected=" << out.E_loc_kV_cm << " kV/cm" << std::endl;
        print([&](std::ostream& os) {
            os << "applyCalorimetry: E-field SCE offsets at p_corr=("
               << p_corr.X() << "," << p_corr.Y() << "," << p_corr.Z() << ") = ("
               << fo.X() << "," << fo.Y() << "," << fo.Z() << ")"
               << " -> E_loc_kV_cm=" << out.E_loc_kV_cm;
        });
    } else {
        std::cout << "[ImageCorrections] E-field: nominal=" << nominal_E_kV_cm
                  << " kV/cm (no SCE offsets)" << std::endl;
        print([&](std::ostream& os) {
            os << "applyCalorimetry: no E-field SCE at p_corr=("
               << p_corr.X() << "," << p_corr.Y() << "," << p_corr.Z() << ")"
               << " E_loc_kV_cm=" << out.E_loc_kV_cm;
        });
    }
    if (calo_alg && clocks && pitch_cm > 0.0) {
        // For debugging: approximate drift time in microseconds
        double const tick_period_us = clocks->TPCClock().TickPeriod();
        double const drift_time_us  = hit.PeakTime() * tick_period_us - T0_ns * 1.0e-3;

        // Always emit drift time to std::cout for quick inspection
        std::cout << "[ImageCorrections] drift_time_us≈" << drift_time_us
                  << " (T0_ns=" << T0_ns << ", tick_period_us=" << tick_period_us << ")"
                  << std::endl;

        out.dEdx_MeV_cm = calo_alg->dEdx_AREA(hit, pitch_cm, T0_ns);
        out.E_hit_MeV = out.dEdx_MeV_cm * pitch_cm;
        print([&](std::ostream& os) {
            os << "applyCalorimetry: channel=" << hit.Channel()
               << " peakTime=" << hit.PeakTime()
               << " pitch_cm=" << pitch_cm
               << " T0_ns=" << T0_ns
               << " drift_time_us≈" << drift_time_us
               << " dEdx=" << out.dEdx_MeV_cm << " MeV/cm"
               << " E_hit=" << out.E_hit_MeV << " MeV";
        });
    } else {
        print([&](std::ostream& os) {
            os << "applyCalorimetry: skipping dEdx/E_hit calc:"
               << " calo_alg=" << (calo_alg ? "set" : "null")
               << " detprop=" << (detprop ? "set" : "null")
               << " clocks=" << (clocks ? "set" : "null")
               << " pitch_cm=" << pitch_cm;
        });
    }
    return out;
}

}
}

#endif
