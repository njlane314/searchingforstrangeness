#ifndef CENTROID_H
#define CENTROID_H

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "detinfo/DetectorPropertiesService.h"
#include "TVector3.h"

namespace common {

TVector3 calculateChargeCentroid(
    const art::Event& e,
    const std::vector<art::Ptr<recob::Hit>>& hits
) {
    double total_charge = 0.0;
    double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
    auto det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
    auto geom = art::ServiceHandle<geo::Geometry>()->provider();
    for (const auto& hit : hits) {
        double charge = hit->Integral();
        total_charge += charge;
        double x = det_props.ConvertTicksToX(hit->PeakTime(), hit->WireID().Plane);
        geo::WireID wire_id(hit->WireID());
        auto wire_pos = geom->WireIDToWireGeo(wire_id).GetCenter();
        double y = wire_pos.Y();
        double z = wire_pos.Z();
        x_sum += charge * x;
        y_sum += charge * y;
        z_sum += charge * z;
    }
    if (total_charge <= 0) {
        return TVector3(-999, -999, -999);
    }
    return TVector3(x_sum / total_charge, y_sum / total_charge, z_sum / total_charge);
}

}

#endif