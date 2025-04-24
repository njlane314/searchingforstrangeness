#ifndef SELECTION_FIDUCIALVOLUME_CXX
#define SELECTION_FIDUCIALVOLUME_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "../CommonDefs/Corrections.h" 

namespace selection
{
    class FiducialVolumeSelection : public SelectionToolBase {
    public:
        FiducialVolumeSelection(const fhicl::ParameterSet& pset);
        ~FiducialVolumeSelection() {}

        void configure(fhicl::ParameterSet const & pset);
        bool selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image<float>>& calo_images, 
                                const std::vector<image::Image<int>>& reco_images, 
                                const std::vector<image::Image<int>>& label_images);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);
        
        bool isFiducial(const double x[3]) const;

    private:
        double m_fidvol_xstart, m_fidvol_xend;
        double m_fidvol_ystart, m_fidvol_yend;
        double m_fidvol_zstart, m_fidvol_zend;

        art::ServiceHandle<geo::Geometry> _geometry;
        art::ServiceHandle<detinfo::DetectorPropertiesService> _detprop;

        double m_reco_nu_vtx_x; 
        double m_reco_nu_vtx_y;
        double m_reco_nu_vtx_z;
    };

    FiducialVolumeSelection::FiducialVolumeSelection(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void FiducialVolumeSelection::configure(fhicl::ParameterSet const & pset) {
        m_fidvol_xstart = pset.get<double>("fidvolXstart");
        m_fidvol_xend = pset.get<double>("fidvolXend");
        m_fidvol_ystart = pset.get<double>("fidvolYstart");
        m_fidvol_yend = pset.get<double>("fidvolYend");
        m_fidvol_zstart = pset.get<double>("fidvolZstart");
        m_fidvol_zend = pset.get<double>("fidvolZend");
    }

    void FiducialVolumeSelection::setBranches(TTree* _tree) {
        _tree->Branch("reco_nu_vtx_x", &m_reco_nu_vtx_x, "reco_nu_vtx_x/D"); 
        _tree->Branch("reco_nu_vtx_y", &m_reco_nu_vtx_y, "reco_nu_vtx_y/D"); 
        _tree->Branch("reco_nu_vtx_z", &m_reco_nu_vtx_z, "reco_nu_vtx_z/D"); 
    }

    void FiducialVolumeSelection::resetTTree(TTree* _tree) {
        m_reco_nu_vtx_x = std::numeric_limits<double>::min(); 
        m_reco_nu_vtx_y = std::numeric_limits<double>::min();
        m_reco_nu_vtx_z = std::numeric_limits<double>::min();
    }

    bool FiducialVolumeSelection::isFiducial(const double x[3]) const {
        geo::TPCGeo const &thisTPC = _geometry->TPC();
        geo::BoxBoundedGeo theTpcGeo = thisTPC.ActiveBoundingBox();
        std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
        bool is_x = x[0] > (bnd[0] + m_fidvol_xstart) && x[0] < (bnd[1] - m_fidvol_xend);
        bool is_y = x[1] > (bnd[2] + m_fidvol_ystart) && x[1] < (bnd[3] - m_fidvol_yend);
        bool is_z = x[2] > (bnd[4] + m_fidvol_zstart) && x[2] < (bnd[5] - m_fidvol_zend);
        return is_x && is_y && is_z;
    }

    bool FiducialVolumeSelection::selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image<float>>& calo_images, 
                                const std::vector<image::Image<int>>& reco_images, 
                                const std::vector<image::Image<int>>& label_images) {
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

                    bool in_fiducial = isFiducial(xyz);

                    return in_fiducial;
                }
            }
        }

        m_reco_nu_vtx_x = std::numeric_limits<double>::min();
        m_reco_nu_vtx_y = std::numeric_limits<double>::min();
        m_reco_nu_vtx_z = std::numeric_limits<double>::min();
        return false;
    }

    DEFINE_ART_CLASS_TOOL(FiducialVolumeSelection)
}

#endif