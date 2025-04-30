#ifndef SELECTION_STRANGENESS_CXX
#define SELECTION_STRANGENESS_CXX

#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "../CommonDefs/Corrections.h"
#include "../CommonDefs/Image.h"
#include <vector>
#include <limits>

namespace selection 
{
    class StrangenessSelection : public SelectionToolBase {
    public:
        StrangenessSelection(const fhicl::ParameterSet& pset);
        ~StrangenessSelection() {}

        void configure(fhicl::ParameterSet const& pset);
        bool selectEvent(art::Event const& e, 
                         const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                         const std::vector<image::Image<float>>& calo_images, 
                         const std::vector<image::Image<int>>& reco_images, 
                         const std::vector<image::Image<int>>& label_images);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);

    private:
        double m_fidvol_xstart, m_fidvol_xend;
        double m_fidvol_ystart, m_fidvol_yend;
        double m_fidvol_zstart, m_fidvol_zend;
        art::ServiceHandle<geo::Geometry> m_geometry;
        art::ServiceHandle<detinfo::DetectorPropertiesService> m_detprop;
        double m_reco_nu_vtx_x; 
        double m_reco_nu_vtx_y;
        double m_reco_nu_vtx_z;
        bool m_in_fiducial;

        art::InputTag m_CLSproducer;
        int m_nhits_u, m_nhits_v, m_nhits_w;
        float m_charge_u, m_charge_v, m_charge_w;
        float m_wirerange_u, m_wirerange_v, m_wirerange_w;
        float m_timerange_u, m_timerange_v, m_timerange_w;

        std::vector<float> calo_pixels_u;
        std::vector<float> calo_pixels_v;
        std::vector<float> calo_pixels_w;
        std::vector<int> reco_pixels_u;
        std::vector<int> reco_pixels_v;
        std::vector<int> reco_pixels_w;
        std::vector<int> label_pixels_u;
        std::vector<int> label_pixels_v;
        std::vector<int> label_pixels_w;

        bool isFiducial(const double x[3]) const;
        void computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
    };

    StrangenessSelection::StrangenessSelection(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void StrangenessSelection::configure(fhicl::ParameterSet const& pset) {
        m_fidvol_xstart = pset.get<double>("fidvolXstart");
        m_fidvol_xend = pset.get<double>("fidvolXend");
        m_fidvol_ystart = pset.get<double>("fidvolYstart");
        m_fidvol_yend = pset.get<double>("fidvolYend");
        m_fidvol_zstart = pset.get<double>("fidvolZstart");
        m_fidvol_zend = pset.get<double>("fidvolZend");

        m_CLSproducer = pset.get<art::InputTag>("CLSproducer");
    }

    bool StrangenessSelection::isFiducial(const double x[3]) const {
        geo::TPCGeo const &thisTPC = m_geometry->TPC();
        geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
        std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
        bool is_x = x[0] > (bnd[0] + m_fidvol_xstart) && x[0] < (bnd[1] - m_fidvol_xend);
        bool is_y = x[1] > (bnd[2] + m_fidvol_ystart) && x[1] < (bnd[3] - m_fidvol_yend);
        bool is_z = x[2] > (bnd[4] + m_fidvol_zstart) && x[2] < (bnd[5] - m_fidvol_zend);
        return is_x && is_y && is_z;
    }

    void StrangenessSelection::computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        auto cluster_handle = e.getValidHandle<std::vector<recob::Cluster>>(m_CLSproducer);
        art::FindManyP<recob::Hit> cluster_hits(cluster_handle, e, m_CLSproducer);
        
        std::map<int, std::vector<art::Ptr<recob::Hit>>> slicehits_plane;
        for (const auto& pfp : pfp_pxy_v) {
            auto clusters = pfp.get<recob::Cluster>();
            for (const auto& cluster : clusters) {
                const std::vector<art::Ptr<recob::Hit>>& hits = cluster_hits.at(cluster.key());
                for (const auto& hit : hits) {
                    int plane = hit->WireID().planeID().Plane;
                    if (plane < 0 || plane > 2) continue;  
                    slicehits_plane[plane].push_back(hit);
                }
            }
        }

        for (int plane = 0; plane < 3; ++plane) {
            auto& hits = slicehits_plane[plane];
            float min_wire = std::numeric_limits<float>::max(), max_wire = std::numeric_limits<float>::lowest();
            float min_time = std::numeric_limits<float>::max(), max_time = std::numeric_limits<float>::lowest();
            float charge_sum = 0.0f;

            for (const auto& hit : hits) {
                float wire = hit->WireID().Wire, time = hit->PeakTime();
                min_wire = std::min(min_wire, wire); max_wire = std::max(max_wire, wire);
                min_time = std::min(min_time, time); max_time = std::max(max_time, time);
                charge_sum += hit->Integral();
            }

            int nhits = hits.size();
            float wirerange = nhits ? max_wire - min_wire : 0.0f;
            float timerange = nhits ? max_time - min_time : 0.0f;

            if (plane == 0) {
                m_nhits_u = nhits; m_charge_u = charge_sum; m_wirerange_u = wirerange; m_timerange_u = timerange;
            } else if (plane == 1) {
                m_nhits_v = nhits; m_charge_v = charge_sum; m_wirerange_v = wirerange; m_timerange_v = timerange;
            } else {
                m_nhits_w = nhits; m_charge_w = charge_sum; m_wirerange_w = wirerange; m_timerange_w = timerange;
            }
        }
    }

    bool StrangenessSelection::selectEvent(art::Event const& e, 
                                           const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                           const std::vector<image::Image<float>>& calo_images, 
                                           const std::vector<image::Image<int>>& reco_images, 
                                           const std::vector<image::Image<int>>& label_images) {
        m_in_fiducial = false;
        m_reco_nu_vtx_x = std::numeric_limits<double>::min();
        m_reco_nu_vtx_y = std::numeric_limits<double>::min();
        m_reco_nu_vtx_z = std::numeric_limits<double>::min();

        for (const auto& pfp_pxy : pfp_pxy_v) {
            if (pfp_pxy->IsPrimary()) {
                auto vertices = pfp_pxy.get<recob::Vertex>();
                if (vertices.size() == 1) {
                    const auto& vertex = vertices[0];
                    double xyz[3] = {};
                    vertex->XYZ(xyz);

                    float x = static_cast<float>(xyz[0]);
                    float y = static_cast<float>(xyz[1]);
                    float z = static_cast<float>(xyz[2]);

                    common::ApplySCECorrectionXYZ(x, y, z);

                    xyz[0] = static_cast<double>(x);
                    xyz[1] = static_cast<double>(y);
                    xyz[2] = static_cast<double>(z);

                    m_reco_nu_vtx_x = xyz[0];
                    m_reco_nu_vtx_y = xyz[1];
                    m_reco_nu_vtx_z = xyz[2];

                    m_in_fiducial = this->isFiducial(xyz);
                    break;
                }
            }
        }

        this->computeFeatures(e, pfp_pxy_v);

        if (calo_images.size() == 3 && reco_images.size() == 3 && label_images.size() == 3) {
            for (const auto& img : calo_images) {
                if (img.view() == geo::kU) {
                    calo_pixels_u = img.data();
                } else if (img.view() == geo::kV) {
                    calo_pixels_v = img.data();
                } else if (img.view() == geo::kW) {
                    calo_pixels_w = img.data();
                }
            }

            for (const auto& img : reco_images) {
                if (img.view() == geo::kU) {
                    reco_pixels_u = img.data();
                } else if (img.view() == geo::kV) {
                    reco_pixels_v = img.data();
                } else if (img.view() == geo::kW) {
                    reco_pixels_w = img.data();
                }
            }

            for (const auto& img : label_images) {
                if (img.view() == geo::kU) {
                    label_pixels_u = img.data();
                } else if (img.view() == geo::kV) {
                    label_pixels_v = img.data();
                } else if (img.view() == geo::kW) {
                    label_pixels_w = img.data();
                }
            }
        }

        return true;
    }

    void StrangenessSelection::setBranches(TTree* _tree) {
        _tree->Branch("in_fiducial", &m_in_fiducial, "in_fiducial/O");
        _tree->Branch("nhits_u", &m_nhits_u, "nhits_u/I");
        _tree->Branch("nhits_v", &m_nhits_v, "nhits_v/I");
        _tree->Branch("nhits_w", &m_nhits_w, "nhits_w/I");
        _tree->Branch("charge_u", &m_charge_u, "charge_u/F");
        _tree->Branch("charge_v", &m_charge_v, "charge_v/F");
        _tree->Branch("charge_w", &m_charge_w, "charge_w/F");
        _tree->Branch("wirerange_u", &m_wirerange_u, "wirerange_u/F");
        _tree->Branch("wirerange_v", &m_wirerange_v, "wirerange_v/F");
        _tree->Branch("wirerange_w", &m_wirerange_w, "wirerange_w/F");
        _tree->Branch("timerange_u", &m_timerange_u, "timerange_u/F");
        _tree->Branch("timerange_v", &m_timerange_v, "timerange_v/F");
        _tree->Branch("timerange_w", &m_timerange_w, "timerange_w/F");
        _tree->Branch("calo_pixels_u", &calo_pixels_u);
        _tree->Branch("calo_pixels_v", &calo_pixels_v);
        _tree->Branch("calo_pixels_w", &calo_pixels_w);
        _tree->Branch("reco_pixels_u", &reco_pixels_u);
        _tree->Branch("reco_pixels_v", &reco_pixels_v);
        _tree->Branch("reco_pixels_w", &reco_pixels_w);
        _tree->Branch("label_pixels_u", &label_pixels_u);
        _tree->Branch("label_pixels_v", &label_pixels_v);
        _tree->Branch("label_pixels_w", &label_pixels_w);
    }

    void StrangenessSelection::resetTTree(TTree* _tree) {
        m_reco_nu_vtx_x = std::numeric_limits<double>::min(); 
        m_reco_nu_vtx_y = std::numeric_limits<double>::min();
        m_reco_nu_vtx_z = std::numeric_limits<double>::min();
        m_in_fiducial = false;
        m_nhits_u = m_nhits_v = m_nhits_w = 0;
        m_charge_u = m_charge_v = m_charge_w = 0.0f;
        m_wirerange_u = m_wirerange_v = m_wirerange_w = 0.0f;
        m_timerange_u = m_timerange_v = m_timerange_w = 0.0f;
        calo_pixels_u.clear();
        calo_pixels_v.clear();
        calo_pixels_w.clear();
        reco_pixels_u.clear();
        reco_pixels_v.clear();
        reco_pixels_w.clear();
        label_pixels_u.clear();
        label_pixels_v.clear();
        label_pixels_w.clear();
    }

    DEFINE_ART_CLASS_TOOL(StrangenessSelection)
}

#endif