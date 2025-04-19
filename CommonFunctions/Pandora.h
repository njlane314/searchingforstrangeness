#ifndef PANDORAFUNCS_H
#define PANDORAFUNCS_H

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraInterface/LArPandoraGeometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TVector3.h"
#include <cmath>
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larcore/Geometry/WireReadout.h"

namespace common
{
    enum PandoraView {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

    constexpr size_t N_VIEWS = 3;

    PandoraView GetPandoraView(const art::Ptr<recob::Hit> &hit)
    {
        const geo::WireID hit_wire(hit->WireID());
        const geo::View_t hit_view(hit->View());
        const geo::View_t pandora_view(lar_pandora::LArPandoraGeometry::GetGlobalView(hit_wire.Cryostat, hit_wire.TPC, hit_view));

        if (pandora_view == geo::kW || pandora_view == geo::kY)
            return TPC_VIEW_W;
        else if (pandora_view == geo::kU)
            return TPC_VIEW_U;
        else if (pandora_view == geo::kV)
            return TPC_VIEW_V;
        else
            throw cet::exception("PandoraFuncs") << "wire view not recognised";
    }

    float YZtoU(const float y_coord, const float z_coord)
    {
        const float m_uWireAngle = 1.04719758034;
        return (z_coord * std::cos(m_uWireAngle)) - (y_coord * std::sin(m_uWireAngle));
    }

    float YZtoV(const float y_coord, const float z_coord)
    {
        const float m_vWireAngle = -1.04719758034;
        return (z_coord * std::cos(m_vWireAngle)) - (y_coord * std::sin(m_vWireAngle));
    }

    float YZtoW(const float y_coord, const float z_coord)
    {
        const float m_wWireAngle = 0.0;
        return (z_coord * std::cos(m_wWireAngle)) - (y_coord * std::sin(m_wWireAngle));
    }

    TVector3 ProjectToWireView(const float input_x, const float input_y, const float input_z, const PandoraView pandora_view)
    {
        const float x_coord = input_x;
        const float y_coord = input_y;
        const float z_coord = input_z;

        return TVector3(x_coord, 0.f, pandora_view == TPC_VIEW_U ? YZtoU(y_coord, z_coord) : pandora_view == TPC_VIEW_V ? YZtoV(y_coord, z_coord) : YZtoW(y_coord, z_coord));
    }

    TVector3 GetPandoraHitPosition(const art::Event& e, art::Ptr<recob::Hit> hit, PandoraView view) 
    {
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();
        auto const detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);

        const geo::WireID hit_wire = hit->WireID();
        const double hit_time = hit->PeakTime();

        double x_coord = detp.ConvertTicksToX(hit_time, hit_wire.Plane, hit_wire.TPC, hit_wire.Cryostat);

        auto const wire_geo = wireReadout.Wire(hit_wire);
        auto const center = wire_geo.GetCenter();

        TVector3 xyz(center.X(), center.Y(), center.Z());

        xyz.SetX(x_coord);
        return xyz;
    }
} 

#endif