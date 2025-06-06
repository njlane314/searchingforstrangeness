#include "HitGeometry.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "TVector3.h"

#include <vector>
#include <utility>

namespace common
{
    void CalculateHitChargeCentroid(const art::Event& event,
                                    PandoraView view, 
                                    const std::vector<art::Ptr<recob::Hit>>& hits,
                                    std::pair<double, double>& charge_centroid)
    {
        double sum_charge = 0.0;
        double weighted_sum_drift = 0.0;
        double weighted_sum_wire = 0.0;

        for (const auto& hit : hits) {
            if (GetPandoraView(hit) != view) continue;

            const double charge = hit->Integral();
            const TVector3 pos = common::GetPandoraHitPosition(event, hit, view);

            sum_charge += charge;
            weighted_sum_drift += pos.X() * charge;
            weighted_sum_wire += pos.Z() * charge;
        }

        if (sum_charge > 1e-9) {
            charge_centroid.first = weighted_sum_drift / sum_charge;
            charge_centroid.second = weighted_sum_wire / sum_charge;
        } else {
            charge_centroid.first = -9999.0;
            charge_centroid.second = -9999.0;
        }
    }
}