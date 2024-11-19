#ifndef REGION_H
#define REGION_H

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "CommonFunctions/Types.h"

namespace common
{
    void addDaughters(const ProxyPfpElem_t &pfp_pxy,
                                           const ProxyPfpColl_t &pfp_pxy_col,
                                           std::vector<ProxyPfpElem_t> &slice_v)
    {
        std::map<unsigned int, unsigned int> pfp_map;

        unsigned int p = 0;
        for (const auto &pfp_pxy : pfp_pxy_col)
        {
            pfp_map[pfp_pxy->Self()] = p;
            p++;
        }
    
        auto daughters = pfp_pxy->Daughters();

        slice_v.push_back(pfp_pxy);

        std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;

        for (auto const &daughterid : daughters)
        {

            if (pfp_map.find(daughterid) == pfp_map.end())
                continue;

            auto pfp_pxy2 = pfp_pxy_col.begin();
            for (size_t j = 0; j < pfp_map.at(daughterid); ++j)
                ++pfp_pxy2;

            common::addDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);

        } 

        return;
    } 

    std::pair<std::vector<art::Ptr<recob::Hit>>, std::vector<ProxyPfpElem_t>> getNuSliceHits(const common::ProxyPfpColl_t& pfp_proxy, 
                                                 const common::ProxyClusColl_t& clus_proxy)
    {
        std::vector<art::Ptr<recob::Hit>> nu_slice_hits;
        std::vector<ProxyPfpElem_t> nu_slice;

        for (const ProxyPfpElem_t& pfp_pxy : pfp_proxy)
        {
            if (!pfp_pxy->IsPrimary()) continue;

            int pdg = abs(pfp_pxy->PdgCode());
            if (pdg == 12 || pdg == 14) 
            {
                nu_slice.clear();
                common::addDaughters(pfp_pxy, pfp_proxy, nu_slice); 
                break;  
            }
        }

        for (const ProxyPfpElem_t& pfp_pxy : nu_slice)
        {
            auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();

            for (auto ass_clus : clus_pxy_v)
            {
                const auto& clus = clus_proxy[ass_clus.key()];
                auto clus_hit_v = clus.get<recob::Hit>();

                nu_slice_hits.insert(nu_slice_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
            }
        }

        return {nu_slice_hits, nu_slice};
    }

    void initialiseChargeMap(
        std::map<common::PandoraView, std::array<float, 2>>& q_centre_map,
        std::map<common::PandoraView, float>& tot_q_map)
    {
        for (const auto& view : {common::TPC_VIEW_U, common::TPC_VIEW_V, common::TPC_VIEW_W}) 
        {
            q_centre_map[view] = {0.0f, 0.0f};  
            tot_q_map[view] = 0.0f;  
        }
    }

    std::tuple<float, unsigned int, unsigned int, unsigned int> getMaxDetectorLimits() 
    {
        const geo::GeometryCore* geom = lar::providerFrom<geo::Geometry>();
        const detinfo::DetectorProperties* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

        const float max_time = detprop->NumberTimeSamples(); 
        const unsigned int max_wire_u = geom->Nwires(geo::PlaneID(0, 0, 0)); 
        const unsigned int max_wire_v = geom->Nwires(geo::PlaneID(0, 1, 0));  
        const unsigned int max_wire_w = geom->Nwires(geo::PlaneID(0, 2, 0)); 

        return {max_time, max_wire_u, max_wire_v, max_wire_w};
    }

    unsigned int getMaxWires(common::PandoraView view, unsigned int max_wire_u, unsigned int max_wire_v, unsigned int max_wire_w)
    {
        switch (view)
        {
            case common::TPC_VIEW_U: return max_wire_u;
            case common::TPC_VIEW_V: return max_wire_v;
            case common::TPC_VIEW_W: return max_wire_w;
            default: return 0;
        }
    }

} 

#endif