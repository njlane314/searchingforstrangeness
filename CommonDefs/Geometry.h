#ifndef GEOMETRYFUNCS_H
#define GEOMETRYFUNCS_H

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "TMatrixDSymEigen.h" 


namespace common
{
    float distance2d(const float& x1, const float& y1,
                    const float& x2, const float& y2)
    {
        return sqrt((x1-x2)*(x1-x2) +
                    (y1-y2)*(y1-y2));
    }

    float distance3d(const float& x1, const float& y1, const float& z1,
                    const float& x2, const float& y2, const float& z2)
    {
        return sqrt((x1-x2)*(x1-x2) +
                    (y1-y2)*(y1-y2) +
                    (z1-z2)*(z1-z2));
    }

    double distance3d(const double& x1, const double& y1, const double& z1,
                    const double& x2, const double& y2, const double& z2)
    {
        return sqrt((x1-x2)*(x1-x2) +
                    (y1-y2)*(y1-y2) +
                    (z1-z2)*(z1-z2));
    }

    float distance3d(const float& x1, const float& y1, const float& z1,
                    const double& x2, const double& y2, const double& z2)
    {
        return sqrt((x1-x2)*(x1-x2) +
                    (y1-y2)*(y1-y2) +
                    (z1-z2)*(z1-z2));
    }

    float distance3d(const double& x1, const double& y1, const double& z1,
                    const float& x2, const float& y2, const float& z2)
    {
        return sqrt((x1-x2)*(x1-x2) +
                    (y1-y2)*(y1-y2) +
                    (z1-z2)*(z1-z2));
    }

    float YZtoPlanecoordinate(const float y, const float z, const int plane)
    {
        auto const* geom = ::lar::providerFrom<geo::Geometry>();
        double _wire2cm = geom->WirePitch(0, 0, 0);
        return geom->WireCoordinate(y, z, geo::PlaneID(0, 0, plane)) * _wire2cm;
    }

    float getPitch(float dir_y, float dir_z, int plane)
    {
        float aux_cos = 1.;
        if (plane == 0)
            aux_cos = dir_y * (-sqrt(3)/2) + dir_z * (1/2);
        if (plane == 1)
            aux_cos = dir_y * (sqrt(3)/2) + dir_z * (1/2);
        if (plane == 2)
            aux_cos = dir_z;

        return 0.3/aux_cos;
    }

    void TrkDirectionAtXYZ(const recob::Track trk, const double x, const double y, const double z, float out[3])
    {
        float min_dist = 100;
        size_t i_min = -1;
        for(size_t i=0; i < trk.NumberTrajectoryPoints(); i++)
        {
            if (trk.HasValidPoint(i))
            { // check this point is valid
                auto point_i = trk.LocationAtPoint(i);
                float distance = common::distance3d((double)point_i.X(), (double)point_i.Y(), (double)point_i.Z(),
                        x, y, z);
                if (distance < min_dist)
                {
                min_dist = distance;
                i_min = i;
                }
            }// if point is valid
        }// for all track points

        auto direction = trk.DirectionAtPoint(i_min);
        out[0] = (float)direction.X();
        out[1] = (float)direction.Y();
        out[2] = (float)direction.Z();

        float norm;
        norm = out[0]*out[0] + out[1]*out[1] + out[2]*out[2];
        if (fabs(norm -1) > 0.001)
            {
            std::cout << "i_min = " << i_min << std::endl;
            std::cout << "minimum distance = " << min_dist << std::endl;
            std::cout << "out[0], out[1], out[2] = " << out[0] << " , " << out[1] << " , " << out[2] << std::endl;
            std::cout << "norm = " << norm << std::endl;
        }
    }

    std::vector<float> polarAngles(float dir_x, float dir_y, float dir_z, size_t axis, size_t plane)
    {
        float dir_y_prime, dir_z_prime;
        if (plane == 0)
        {
            dir_y_prime = dir_y * (1/2) + dir_z * (sqrt(3)/2);
            dir_z_prime = dir_y * (-sqrt(3)/2) + dir_z * (1/2);
        }

        if (plane == 1)
        {
            dir_y_prime = dir_y * (1/2) + dir_z * (-sqrt(3)/2);
            dir_z_prime = dir_y * (sqrt(3)/2) + dir_z * (1/2);
        }

        if (plane == 2)
        {
            dir_y_prime = dir_y;
            dir_z_prime = dir_z;
        }

        std::vector<float> abs_angle;

        if (axis == 0)
        {
            abs_angle.push_back(acos( abs(dir_x)));
            abs_angle.push_back(atan2( abs(dir_y_prime), abs(dir_z_prime)));
        }

        if (axis == 1)
        {
            abs_angle.push_back(acos(abs(dir_y_prime)));
            abs_angle.push_back(atan2(abs(dir_x), abs(dir_z_prime)));
        }

        if (axis == 2)
        {
            abs_angle.push_back(acos(abs(dir_z_prime)));
            abs_angle.push_back(atan2(abs(dir_y_prime), abs(dir_x)));
        }

        return abs_angle;
    }

    std::vector<std::vector<float>> polarAngles(std::vector<float> dir_x, std::vector<float> dir_y, std::vector<float> dir_z, size_t axis, size_t plane)
    {
        std::vector<float> aux_theta_v, aux_phi_v;
        for (size_t i = 0; i < dir_x.size(); i++)
        {
            std::vector<float> aux_angles = polarAngles(dir_x[i], dir_y[i], dir_z[i], axis, plane);
            aux_theta_v.push_back(aux_angles[0]);
            aux_phi_v.push_back(aux_angles[1]);
        }

        std::vector<std::vector<float>> out;
        out.push_back(aux_theta_v);
        out.push_back(aux_phi_v);

        return out;
    }

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
        
        wirecm = geom->WireCoordinate(pt3d[1], pt3d[2], geo::PlaneID(0,0,pl)) * wire2cm;
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
        Project3Dto2D(pt3d,pl, wire2cm, time2cm, ptwire, pttime);
        
        float hitwire = hit->WireID().Wire * wire2cm;
        float hittime = (hit->PeakTime() - detp->TriggerOffset())  * time2cm;
        
        float distance = sqrt( (ptwire - hitwire)*(ptwire - hitwire) + (pttime - hittime)*(pttime - hittime) );

        return distance;
    }

    
} 

#endif