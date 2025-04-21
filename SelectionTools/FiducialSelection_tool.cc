#include <iostream>
#include "SelectionToolBase.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace selection
{
    class FiducialSelection : public SelectionToolBase {
    public:
        FiducialSelection(const fhicl::ParameterSet& pset);
        ~FiducialSelection() {}

        void configure(fhicl::ParameterSet const & pset);
        bool selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& pfp_pxy_v);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);

    private:
        double _fiducial_x_min, _fiducial_x_max;
        double _fiducial_y_min, _fiducial_y_max;
        double _fiducial_z_min, _fiducial_z_max;

        art::ServiceHandle<geo::Geometry> _geometry;
        art::ServiceHandle<detinfo::DetectorPropertiesService> _detprop;

        double m_x_centroid; 
        double m_y_centroid;
        double m_z_centroid;
    };

    FiducialSelection::FiducialSelection(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void FiducialSelection::configure(fhicl::ParameterSet const & pset) {
        _fiducial_x_min = pset.get<double>("fiducial_x_min");
        _fiducial_x_max = pset.get<double>("fiducial_x_max");
        _fiducial_y_min = pset.get<double>("fiducial_y_min");
        _fiducial_y_max = pset.get<double>("fiducial_y_max");
        _fiducial_z_min = pset.get<double>("fiducial_z_min");
        _fiducial_z_max = pset.get<double>("fiducial_z_max");
    }

    void FiducialSelection::setBranches(TTree* _tree) {
        _tree->Branch("x_centroid", &m_x_centroid, "x_centroid/D");
        _tree->Branch("y_centroid", &m_y_centroid, "y_centroid/D");
        _tree->Branch("z_centroid", &m_z_centroid, "z_centroid/D");
    }

    void FiducialSelection::resetTTree(TTree* _tree) {
        m_x_centroid = std::numeric_limits<double>::min(); 
        m_y_centroid = std::numeric_limits<double>::min();
        m_z_centroid = std::numeric_limits<double>::min();
    }

    bool FiducialSelection::selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> slice_hits;
        for (const auto& pfp_pxy : pfp_pxy_v) {
            auto hits = pfp_pxy.get<recob::Hit>();
            slice_hits.insert(slice_hits.end(), hits.begin(), hits.end());
        }

        if (slice_hits.empty()) {
            m_x_centroid = std::numeric_limits<double>::min();
            m_y_centroid = std::numeric_limits<double>::min();
            m_z_centroid = std::numeric_limits<double>::min();
            return false;
        }

        double total_charge = 0.0;
        double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;

        for (const auto& hit : slice_hits) {
            double charge = hit->Integral();
            total_charge += charge;

            unsigned int wire = hit->WireID().Wire;
            unsigned int plane = hit->WireID().Plane;
            double time = hit->PeakTime();

            double x = _detprop->ConvertTicksToX(time, plane);

            geo::WireID wire_id(hit->WireID());
            auto wire_pos = _geometry->WireIDToWireGeo(wire_id).GetCenter();
            double y = wire_pos.Y();
            double z = wire_pos.Z();

            x_sum += charge * x;
            y_sum += charge * y;
            z_sum += charge * z;
        }

        if (total_charge <= 0) {
            m_x_centroid = std::numeric_limits<double>::min();
            m_y_centroid = std::numeric_limits<double>::min();
            m_z_centroid = std::numeric_limits<double>::min();
            return false;
        }

        double x_centroid = x_sum / total_charge;
        double y_centroid = y_sum / total_charge;
        double z_centroid = z_sum / total_charge;

        m_x_centroid = x_centroid;  
        m_y_centroid = y_centroid;
        m_z_centroid = z_centroid;

        bool in_fiducial = (x_centroid >= _fiducial_x_min && x_centroid <= _fiducial_x_max &&
                            y_centroid >= _fiducial_y_min && y_centroid <= _fiducial_y_max &&
                            z_centroid >= _fiducial_z_min && z_centroid <= _fiducial_z_max);

        return in_fiducial;
    }
}

DEFINE_ART_CLASS_TOOL(selection::FiducialSelection)