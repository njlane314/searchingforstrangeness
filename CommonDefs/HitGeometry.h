#ifndef HITGEOMETRY_H
#define HITGEOMETRY_H

#include <array>
#include <cmath>
#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"

namespace common
{
    geo::View_t GetTPCView(const art::Ptr<recob::Hit>& hit) {
        geo::View_t view = hit->View();
        if (view == geo::kU || view == geo::kV || view == geo::kW || view == geo::kY) {
            return view;
        } else {
            throw cet::exception("TPCFuncs") << "Wire view not recognized";
        }
    }

    float ProjectYZToWire(float y, float z, geo::View_t view) {
        static constexpr std::array<float, 3> angles = {1.04719758034f, -1.04719758034f, 0.0f};
        size_t index;
        switch (view) {
            case geo::kU: index = 0; break;
            case geo::kV: index = 1; break;
            case geo::kW: case geo::kY: index = 2; break;
            default: throw cet::exception("TPCFuncs") << "Invalid view for projection";
        }
        float angle = angles[index];
        return z * std::cos(angle) - y * std::sin(angle);
    }

    TVector3 ProjectToWireView(float x, float y, float z, geo::View_t view) {
        return TVector3(x, 0.f, ProjectYZToWire(y, z, view));
    }

    TVector3 GetTPCHitPosition(const art::Event& e, const art::Ptr<recob::Hit>& hit, geo::View_t view) {
        art::ServiceHandle<geo::Geometry> geo;
        auto const* det = lar::providerFrom<detinfo::DetectorPropertiesService>();
        
        geo::WireID wire = hit->WireID();
        double x = det->ConvertTicksToX(hit->PeakTime() + 0.5, wire.Plane, wire.TPC, wire.Cryostat);
        TVector3 center = geo->Cryostat(wire.Cryostat).TPC(wire.TPC).Plane(wire.Plane).Wire(wire.Wire).GetCenter();
        
        return ProjectToWireView(x, center.Y(), center.Z(), view);
    }
}

#endif