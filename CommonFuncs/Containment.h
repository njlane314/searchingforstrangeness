#ifndef TRUTHCONTAINMENT_H
#define CONTAINMENT_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"

#include "lardata/Utilities/FindManyInChainP.h"

namespace common
{

    bool isFiducial(const double x[3], 
            const double fid_x_start, const double fid_y_start, const double fid_z_start,
            const double fid_x_end, const double fid_y_end, const double fid_z_end)
    {
        art::ServiceHandle<geo::Geometry> geo;
        geo::TPCGeo const &thisTPC = geo->TPC();
        geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
        std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
        bool is_x = x[0] > (bnd[0] + fid_x_start) && x[0] < (bnd[1] - fid_x_end);
        bool is_y = x[1] > (bnd[2] + fid_y_start) && x[1] < (bnd[3] - fid_y_end);
        bool is_z = x[2] > (bnd[4] + fid_z_start) && x[2] < (bnd[5] - fid_z_end);
        
        return is_x && is_y && is_z;
    }
  
} 

#endif