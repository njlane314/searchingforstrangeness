#ifndef POSITIONTOWIRE_H
#define POSITIONTOWIRE_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TMatrixDSymEigen.h" 

namespace common
{
  /**
   * @brief go from 3D coordinates to 2D coordinates in order to compare 3D reco with hits
   * @input pt3d -> 3d point to be projected
   * @input pl -> which plane are we on?
   * @input wire2cm -> wire 2 cm conversion
   * @input time2cm -> time 2 cm conversion
   * @output wirecm -> pt3d wire coordinate in cm
   * @output timecm -> pt3d time coordinate in cm
   */
  void Project3Dto2D(const TVector3& pt3d, const int& pl,
		     const float& wire2cm, const float& time2cm,
		     float& wirecm, float& timecm) {

    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    
    wirecm = geom->WireCoordinate(pt3d[1],pt3d[2],geo::PlaneID(0,0,pl)) * wire2cm;
    timecm = pt3d[0];

    return;
  }

  /**
   * @brief get hit wire/time in cm
   * @input recob::hit 
   * @output hitwire -> hit wire coordinate in cm
   * @output hittime -> hit time coordinate in cm
   */
  void GetHitWireTime(const art::Ptr<recob::Hit> &hit, 
		      const float& wire2cm, const float& time2cm,		      
		      float& hitwire, float& hittime) {

    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

    hitwire = hit->WireID().Wire * wire2cm;
    hittime = (hit->PeakTime() - detp->TriggerOffset())  * time2cm;

    return;
  }

  /**
   * @brief given a 3D pt and a 2D hit get their distance in 2D on the plane
   * @input pt3d -> 3d point to be projected
   * @input hit -> hit
   * @input wire2cm -> wire 2 cm conversion
   * @input time2cm -> time 2 cm conversion
   * @return 2d distance [cm]
   */
  float HitPtDistance(const TVector3& pt3d, const art::Ptr<recob::Hit> &hit,
		      const float& wire2cm, const float& time2cm) {

    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    // what plane are we on?
    auto pl = hit->WireID().Plane;

    float ptwire, pttime;
    Project3Dto2D(pt3d,pl,wire2cm,time2cm,ptwire,pttime);
    
    float hitwire = hit->WireID().Wire * wire2cm;
    float hittime = (hit->PeakTime() - detp->TriggerOffset())  * time2cm;
    
    float distance = sqrt( (ptwire-hitwire)*(ptwire-hitwire) + (pttime-hittime)*(pttime-hittime) );

    return distance;
  }

}

#endif