#ifndef SELECTION_NEUTRALDECAY_CXX
#define SELECTION_NEUTRALDECAY_CXX

#include <iostream>
#include <vector> 

#include "SelectionToolBase.h"

#include "../CommonFuncs/ProximityClustering.h"
#include "../CommonFuncs/Geometry.h"
#include "../CommonFuncs/SpaceChargeCorrections.h"
#include "../CommonFuncs/PositionToWire.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

namespace selection
{
class NeutralDecaySelection : public SelectionToolBase {

public:

    NeutralDecaySelection(const fhicl::ParameterSet& pset);
    ~NeutralDecaySelection(){};
    
    void configure(fhicl::ParameterSet const & pset);

    bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_proxy);

    void setBranches(TTree* _tree){};

    void resetTTree(TTree* _tree){};

    using ProxyPfpColl_t = common::ProxyPfpColl_t;
    using ProxyPfpElem_t = common::ProxyPfpElem_t;
    
private:

    art::InputTag fPFPproducer;
    art::InputTag fCLSproducer; // cluster associated to PFP
    art::InputTag fSLCproducer; // slice associated to PFP
    art::InputTag fHITproducer; // hit associated to cluster
    art::InputTag fSHRproducer; // shower associated to PFP
    art::InputTag fVTXproducer; // vertex associated to PFP
    art::InputTag fTRKproducer; // track associated to PFP
    art::InputTag fPCAproducer; // PCAxis associated to PFP
    art::InputTag fMCTproducer;
    
    float _wire2cm, _time2cm;
};

NeutralDecaySelection::NeutralDecaySelection(const fhicl::ParameterSet& pset)
{
    fPFPproducer = pset.get<art::InputTag>("PFPproducer");
    fSHRproducer = pset.get<art::InputTag>("SHRproducer");
    fHITproducer = pset.get<art::InputTag>("HITproducer");
    fCLSproducer = pset.get<art::InputTag>("CLSproducer");
    fSLCproducer = pset.get<art::InputTag>("SLCproducer");
    fVTXproducer = pset.get<art::InputTag>("VTXproducer");
    fPCAproducer = pset.get<art::InputTag>("PCAproducer");
    fTRKproducer = pset.get<art::InputTag>("TRKproducer");
    fMCTproducer = pset.get<art::InputTag>("MCTproducer");

    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
    _wire2cm = geom->WirePitch(0,0,0);
    _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
}

void NeutralDecaySelection::configure(fhicl::ParameterSet const & pset)
{
}

bool NeutralDecaySelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_proxy_v)
{
    std::cout << "[NeutralDecaySelection] Starting selection..." << std::endl;
    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(fCLSproducer),
                                                        proxy::withAssociated<recob::Slice>(fSLCproducer),
                                                        proxy::withAssociated<recob::Track>(fTRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(fVTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(fPCAproducer),
                                                        proxy::withAssociated<recob::Shower>(fSHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(fPFPproducer));

    common::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster> >(e, fPFPproducer, proxy::withAssociated<recob::Hit>(fPFPproducer));

    double primary_vertex_xyz[3] = {};
    double primary_vertex_xyz_sce[3] = {};
    bool primary_vertex_found = false;
    std::vector<common::ProxyPfpElem_t> secondary_pfp_list;
    std::vector<std::vector<common::ProxyPfpElem_t>> vertex_groups;

    /*std::cout << "[NeutralDecaySelection] Number of pfp particles in the event is..." << pfp_proxy.size() << std::endl;
    std::cout << "[NeutralDecaySelection] Number of pfp particles in the neutrino slice is..." << pfp_proxy.size() << std::endl;

    std::cout << "[NeutralDecaySelection] Searching for primary vertex..." << std::endl;*/
    for (const ProxyPfpElem_t &pfp : pfp_proxy) 
    {
        if (pfp->IsPrimary()) 
        {
            auto vtx_v = pfp.get<recob::Vertex>();
            if (vtx_v.size() == 1) 
            {
                vtx_v.at(0)->XYZ(primary_vertex_xyz);
                /*std::cout << "[NeutralDecaySelection] Primary vertex found at (" 
                          << primary_vertex_xyz[0] << ", " 
                          << primary_vertex_xyz[1] << ", " 
                          << primary_vertex_xyz[2] << ")" << std::endl;*/

                float x = primary_vertex_xyz[0];
                float y = primary_vertex_xyz[1];
                float z = primary_vertex_xyz[2];
                common::ApplySCECorrectionXYZ(x, y, z);

                primary_vertex_xyz_sce[0] = x;
                primary_vertex_xyz_sce[1] = y;
                primary_vertex_xyz_sce[2] = z;

                primary_vertex_found = true;
                break;
            }
        }
    }

    if (!primary_vertex_found) {
        std::cout << "[NeutralDecaySelection] No primary vertex found, rejecting event." << std::endl;
        return false;
    }

    std::vector<art::Ptr<recob::Hit>> hits_around_primary_vertex;
    for (const ProxyPfpElem_t &pfp : pfp_proxy)
    {   
        if (pfp->IsPrimary())
            continue;

        auto vtx_v = pfp.get<recob::Vertex>();
        if (vtx_v.size() == 1) 
        {
            double vertex_xyz[3];
            vtx_v.at(0)->XYZ(vertex_xyz);

            float x = vertex_xyz[0];
            float y = vertex_xyz[1];
            float z = vertex_xyz[2];
            common::ApplySCECorrectionXYZ(x, y, z);

            vertex_xyz[0] = x;
            vertex_xyz[1] = y;
            vertex_xyz[2] = z;

            /*std::cout << "[NeutralDecaySelection] Secondary vertex found at (" 
                          << vertex_xyz[0] << ", " 
                          << vertex_xyz[1] << ", " 
                          << vertex_xyz[2] << ")" << std::endl;*/

            double distance_to_primary = common::distance3d(primary_vertex_xyz_sce[0], primary_vertex_xyz_sce[1], primary_vertex_xyz_sce[2],
                                                            vertex_xyz[0], vertex_xyz[1], vertex_xyz[2]);

            //std::cout << "[NeutralDecaySelection] Found secondary vertex at distance " << distance_to_primary << " from primary vertex." << std::endl;

            if (distance_to_primary > 0.4 && distance_to_primary < 18.0) 
            {
                secondary_pfp_list.push_back(pfp);
            }

            if (distance_to_primary < 18.0) 
            {   
                std::vector<art::Ptr<recob::Hit>> hit_v;
                auto clus_pxy_v = pfp.get<recob::Cluster>();
                if (clus_pxy_v.size() != 0) 
                {
                    for (auto ass_clus : clus_pxy_v) 
                    {
                        const auto& clus = clus_proxy[ass_clus.key()];
                        auto clus_hit_v = clus.get<recob::Hit>();
                        for (const auto& hit : clus_hit_v)
                            hit_v.push_back(hit);
                    }
                }

                hits_around_primary_vertex.insert(hits_around_primary_vertex.end(), hit_v.begin(), hit_v.end());
            }
        }   
    }

    //std::cout << "[NeutralDecaySelection] Number of secondary vertices: " << secondary_pfp_list.size() << std::endl;
    if (secondary_pfp_list.size() < 2) {
        std::cout << "[NeutralDecaySelection] Less than 2 secondary vertices found, rejecting event." << std::endl;
        return false;
    }

    std::vector<bool> vertex_used(secondary_pfp_list.size(), false);
    for (size_t i = 0; i < secondary_pfp_list.size(); ++i)
    {
        if (vertex_used[i])
            continue;

        std::vector<common::ProxyPfpElem_t> current_group;
        current_group.push_back(secondary_pfp_list[i]);
        vertex_used[i] = true;

        auto vertex1_list = secondary_pfp_list[i].get<recob::Vertex>(); 
        double vertex1_xyz[3];
        vertex1_list.at(0)->XYZ(vertex1_xyz);

        float x1 = vertex1_xyz[0];
        float y1 = vertex1_xyz[1];
        float z1 = vertex1_xyz[2];
        common::ApplySCECorrectionXYZ(x1, y1, z1);

        vertex1_xyz[0] = x1;
        vertex1_xyz[1] = y1;
        vertex1_xyz[2] = z1;

        for (size_t j = i + 1; j < secondary_pfp_list.size(); ++j)
        {
            if (vertex_used[j])
                continue; 

            auto vertex2_list = secondary_pfp_list[j].get<recob::Vertex>(); 
            double vertex2_xyz[3];
            vertex2_list.at(0)->XYZ(vertex2_xyz);

            float x2 = vertex2_xyz[0];
            float y2 = vertex2_xyz[1];
            float z2 = vertex2_xyz[2];
            common::ApplySCECorrectionXYZ(x2, y2, z2);

            vertex2_xyz[0] = x2;
            vertex2_xyz[1] = y2;
            vertex2_xyz[2] = z2;

            double distance_between_vertices = common::distance3d(vertex1_xyz[0], vertex1_xyz[1], vertex1_xyz[2],
                                                      vertex2_xyz[0], vertex2_xyz[1], vertex2_xyz[2]);

            //std::cout << "[NeutralDecaySelection] Distance between secondary vertices: " << distance_between_vertices << std::endl;

            if (distance_between_vertices < 1.6) {
                current_group.push_back(secondary_pfp_list[j]);
                vertex_used[j] = true;
            }
        }

        if (current_group.size() > 1)
            vertex_groups.push_back(current_group);
    }

    //std::cout << "[NeutralDecaySelection] Number of vertex groups: " << vertex_groups.size() << std::endl;

    std::vector<std::vector<unsigned int>> cluster_vector;
    float cell_size = 0.4;
    float clustering_radius = 0.5;

    if(!common::cluster(hits_around_primary_vertex, cluster_vector, cell_size, clustering_radius)) {
        std::cout << "[NeutralDecaySelection] Clustering failed, rejecting event." << std::endl;
        return false;
    }

    //std::cout << "[NeutralDecaySelection] Number of clusters formed " << cluster_vector.size() << std::endl;

    for (const std::vector<common::ProxyPfpElem_t>& group : vertex_groups)
    {
        for (const std::vector<unsigned int>& cluster : cluster_vector)
        {
            bool all_vertices_found_in_cluster = true;

            for (const common::ProxyPfpElem_t& pfp : group)
            {
                auto vertex_list = pfp.get<recob::Vertex>(); 
                if (vertex_list.size() != 1) 
                    continue;

                double vertex_xyz[3];
                vertex_list.at(0)->XYZ(vertex_xyz);

                bool vertex_found_in_cluster = false;
                for (const auto& hit_idx : cluster)
                {
                    const auto& hit = hits_around_primary_vertex[hit_idx];

                    TVector3 pos(vertex_xyz[0], vertex_xyz[1], vertex_xyz[2]);
                    float distance_to_vertex = common::HitPtDistance(pos, hit, _wire2cm, _time2cm);

                    if (distance_to_vertex < 1.6)
                    {
                        //std::cout << "[NeutralDecaySelection] Vertex found in cluster" << std::endl;
                        vertex_found_in_cluster = true;
                        break;
                    }
                }

                if (!vertex_found_in_cluster)
                {
                    all_vertices_found_in_cluster = false;
                    break;
                }
            }

            if (all_vertices_found_in_cluster)
            {
                bool cluster_is_isolated = true;

                for (const auto& hit_idx : cluster)
                {
                    const auto& hit = hits_around_primary_vertex[hit_idx];

                    TVector3 pos(primary_vertex_xyz[0], primary_vertex_xyz[1], primary_vertex_xyz[2]);
                    float distance_to_primary = common::HitPtDistance(pos, hit, _wire2cm, _time2cm);

                    //std::cout << "[NeutralDecaySelection] Distance from hit to primary interaction vertex " << distance_to_primary << std::endl;

                    if (distance_to_primary < 0.4)
                    {
                        cluster_is_isolated = false;
                        break;
                    }
                }

                if (cluster_is_isolated)
                {
                    std::cout << "[NeutralDecaySelection] Event passes selection." << std::endl;
                    return true;
                }
            }
        }
    }

    std::cout << "[NeutralDecaySelection] Event does not pass selection." << std::endl;
    return false;
}


DEFINE_ART_CLASS_TOOL(NeutralDecaySelection)
} 

#endif
