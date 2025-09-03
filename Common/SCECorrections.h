#ifndef SCECORRECTIONSFUNCS_H
#define SCECORRECTIONSFUNCS_H

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

namespace common {

void ApplySCEMappingXYZ(float &x, float &y, float &z) {

    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableSimSpatialSCE() == true) {
        auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
        x -= offset.X();
        y += offset.Y();
        z += offset.Z();
    }
}

void ApplySCECorrectionXYZ(float &x, float &y, float &z) {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableCalSpatialSCE() == true) {

        auto offset = SCE->GetCalPosOffsets(geo::Point_t(x, y, z));
        x -= offset.X();
        y += offset.Y();
        z += offset.Z();
    }
}

float GetSCECorrTrackLength(const art::Ptr<recob::Track> &trk) {

    float SCElength = 0.;
    int previousvalidpoint = -1;

    for (size_t i = 0; i < trk->NumberTrajectoryPoints(); i++) {
        if (trk->HasValidPoint(i)) {

            if (previousvalidpoint >= 0) {
                auto point1 = trk->LocationAtPoint(i);
                auto point0 = trk->LocationAtPoint(previousvalidpoint);

                float point1X = point1.X();
                float point1Y = point1.Y();
                float point1Z = point1.Z();
                ApplySCECorrectionXYZ(point1X, point1Y, point1Z);
                float point0X = point0.X();
                float point0Y = point0.Y();
                float point0Z = point0.Z();
                ApplySCECorrectionXYZ(point0X, point0Y, point0Z);

                float distance3D = sqrt((point1X - point0X) * (point1X - point0X) + (point1Y - point0Y) * (point1Y - point0Y) + (point1Z - point0Z) * (point1Z - point0Z));

                SCElength += distance3D;
            }
            previousvalidpoint = i;
        }
    }

    return SCElength;
}

void True2RecoMappingXYZ(float &t, float &x, float &y, float &z) {
    ApplySCEMappingXYZ(x, y, z);

    auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
    float _xtimeoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);

    x += _xtimeoffset;
    x += 0.6;
}

void ApplySCEMappingXYZ(float x, float y, float z, float out[3]) {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableSimSpatialSCE() == true) {
        auto offset = SCE->GetPosOffsets(geo::Point_t(x, y, z));
        out[0] = (x - offset.X());
        out[1] = (y + offset.Y());
        out[2] = (z + offset.Z());
    }
}

void ApplySCECorrectionXYZ(float x, float y, float z, float out[3]) {
    auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    if (SCE->EnableCalSpatialSCE() == true) {
        auto offset = SCE->GetCalPosOffsets(geo::Point_t(x, y, z));
        out[0] = (x - offset.X());
        out[1] = (y + offset.Y());
        out[2] = (z + offset.Z());
    }
}

float x_offset(float t) {
    auto const &detProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    auto const &detClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    double g4Ticks = detClocks->TPCG4Time2Tick(t) + detProperties->GetXTicksOffset(0, 0, 0) - detProperties->TriggerOffset();
    float xoffset = detProperties->ConvertTicksToX(g4Ticks, 0, 0, 0);
    xoffset += 0.6;
    return xoffset;
}

void True2RecoMappingXYZ(float t, float x, float y, float z, float out[3]) {
    ApplySCEMappingXYZ(x, y, z, out);
    float _xoffset = x_offset(t);
    out[0] += _xoffset;
}

float GetLocalEFieldMag(const float x, const float y, const float z) {

    const detinfo::DetectorProperties *detprop = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
    auto const *sce = lar::providerFrom<spacecharge::SpaceChargeService>();

    double E_field_nominal = detprop->Efield();

    geo::Vector_t E_field_offsets = {0., 0., 0.};
    E_field_offsets = sce->GetCalEfieldOffsets(geo::Point_t{x, y, z});
    TVector3 E_field_vector = {E_field_nominal * (1 + E_field_offsets.X()), E_field_nominal * E_field_offsets.Y(), E_field_nominal * E_field_offsets.Z()};
    float E_field = E_field_vector.Mag();

    return E_field;
}

float GetdEdxfromdQdx(const float dqdx, const float x, const float y, const float z, const float dedxfixed, const float adctoe) {
    auto efield = common::GetLocalEFieldMag(x, y, z);
    float B = 0.212 / (1.383 * efield);
    float r = log(dedxfixed * B + 0.93) / (dedxfixed * B);
    return dqdx * adctoe * (23.6 / 1e6) / r;
}

std::vector<float> GetdEdxfromdQdx(const std::vector<float> dqdx_v,
                                   const std::vector<float> x_v,
                                   const std::vector<float> y_v,
                                   const std::vector<float> z_v,
                                   const float dedxfixed,
                                   const float adctoe) {

    std::vector<float> dedx_v;

    if ((x_v.size() < dqdx_v.size()) || (y_v.size() < dqdx_v.size()) || (z_v.size() < dqdx_v.size())) {
        std::cout << "ERROR. Vector size does not match in CalorimetryAnalysis_tool [common]" << std::endl;
        return dedx_v;
    }

    for (size_t i = 0; i < dqdx_v.size(); i++) {
        dedx_v.push_back(common::GetdEdxfromdQdx(dqdx_v[i], x_v[i], y_v[i], z_v[i], dedxfixed, adctoe));
    }

    return dedx_v;
}

} // namespace common

#endif