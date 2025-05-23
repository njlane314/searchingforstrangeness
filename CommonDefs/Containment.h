#ifndef CONTAINMENT_H
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

    const double FVxmin = 0.0;
    const double FVxmax = 256.35;
    const double FVymin = -115.53;
    const double FVymax = 117.47;
    const double FVzmin = 0.1;
    const double FVzmax = 1036.9;

    bool point_inside_fv(const double x[3], 
        const double fid_x_start=FVxmin, const double fid_y_start=FVymin, const double fid_z_start=FVzmin,
        const double fid_x_end=FVxmax, const double fid_y_end=FVymax, const double fid_z_end=FVzmax)
    {
        if(x[0] > FVxmax || x[0] < FVxmin) return false;
        if(x[1] > FVymax || x[1] < FVymin) return false;
        if(x[2] > FVzmax || x[2] < FVzmin) return false;
        return true;
    }

    /*bool point_inside_fv(const double x[3], 
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
    }*/

    bool truth_contained(const float& FVxS, const float& FVyS, const float& FVzS,
		      const float& FVxE, const float& FVyE, const float& FVzE,
		      const std::vector<sim::MCShower> &inputMCShower,
		      const std::vector<sim::MCTrack> &inputMCTrack ) 
    {
        // require truth-containment by
        // (1) requiring the vertex is in the FV
        // (2) require all MCTracks are contained within the FV
        // (3) require all MCShowers to deposit > some fraction of energy in the FV
        
        for (auto mcs : inputMCShower) {
            if (mcs.Process() == "primary" || (mcs.MotherPdgCode() == 111 && mcs.Process() == "Decay" && mcs.MotherProcess() == "primary") ) {
                float contained = mcs.DetProfile().E() / mcs.Start().E();
                float edep = mcs.DetProfile().E();
                if ( (contained < 0.6) && (edep < 100.) ) {
                    return false;
                }
            }// if primary
        }
        
        for (auto mct : inputMCTrack) {
        
            if (mct.Process() == "primary") {
                //is the start point in the FV?
                sim::MCStep mc_step_track_start = mct.Start();
                sim::MCStep mc_step_track_end   = mct.End();
                
                double start[3];
                start[0] = mc_step_track_start.X();
                start[1] = mc_step_track_start.Y();
                start[2] = mc_step_track_start.Z();
                
                double end[3];
                end[0] = mc_step_track_end.X();
                end[1] = mc_step_track_end.Y();
                end[2] = mc_step_track_end.Z();
                
                
                if (point_inside_fv(start,FVxS,FVyS,FVzS,FVxE,FVyE,FVzE) == false) {
                    return false;
                }
                if (point_inside_fv(end  ,FVxS,FVyS,FVzS,FVxE,FVyE,FVzE) == false) {
                    return false;
                }
            }
        }
        
        return true;
    } 

} 

#endif