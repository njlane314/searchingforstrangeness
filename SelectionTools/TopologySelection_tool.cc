#ifndef SELECTION_TOPOLOGYSELECTION_H
#define SELECTION_TOPOLOGYSELECTION_H

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <limits>
#include <Eigen/Dense>
#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace selection
{
    class TopologySelection : public SelectionToolBase {
    public:
        TopologySelection(const fhicl::ParameterSet& pset);
        ~TopologySelection() {}

        void configure(fhicl::ParameterSet const & pset);
        bool selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& pfp_pxy_v);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);

    private:
        std::string _wire_tag;    // Tag for wire data product
        double _threshold;        // Threshold for lambda2 / lambda1 ratio

        art::ServiceHandle<geo::Geometry> _geometry;
        art::ServiceHandle<detinfo::DetectorPropertiesService> _detprop;

        double m_lambda1;         // First eigenvalue
        double m_lambda2;         // Second eigenvalue
        double m_lambda3;         // Third eigenvalue
        double m_ratio;           // Ratio of lambda2 / lambda1
    };

    TopologySelection::TopologySelection(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void TopologySelection::configure(fhicl::ParameterSet const & pset) {
        _wire_tag = pset.get<std::string>("wire_tag");
        _threshold = pset.get<double>("threshold");
    }

    void TopologySelection::setBranches(TTree* _tree) {
        _tree->Branch("lambda1", &m_lambda1, "lambda1/D");
        _tree->Branch("lambda2", &m_lambda2, "lambda2/D");
        _tree->Branch("lambda3", &m_lambda3, "lambda3/D");
        _tree->Branch("ratio", &m_ratio, "ratio/D");
    }

    void TopologySelection::resetTTree(TTree* _tree) {
        m_lambda1 = std::numeric_limits<double>::lowest();
        m_lambda2 = std::numeric_limits<double>::lowest();
        m_lambda3 = std::numeric_limits<double>::lowest();
        m_ratio = -1.0;
    }

    bool TopologySelection::selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& pfp_pxy_v) {
        // Retrieve wires
        auto wire_handle = e.getValidHandle<std::vector<recob::Wire>>(_wire_tag);
        std::vector<art::Ptr<recob::Wire>> wire_vec;
        art::fill_ptr_vector(wire_vec, wire_handle);

        // Create wire map: (Plane, Wire) -> wire pointer
        std::map<std::pair<unsigned int, unsigned int>, art::Ptr<recob::Wire>> wire_map;
        for (const auto& wire_ptr : wire_vec) {
            auto channel = wire_ptr->Channel();
            auto wids = _geometry->ChannelToWire(channel);
            if (!wids.empty()) {
                auto wid = wids[0]; // Assuming one wire per channel
                wire_map[{wid.Plane, wid.Wire}] = wire_ptr;
            }
        }

        // Collect hits from slice
        std::vector<art::Ptr<recob::Hit>> slice_hits;
        for (const auto& pfp_pxy : pfp_pxy_v) {
            auto hits = pfp_pxy.get<recob::Hit>();
            slice_hits.insert(slice_hits.end(), hits.begin(), hits.end());
        }

        // Collect points from ROIs
        std::vector<double> x_pos, y_pos, z_pos, weights;
        for (const auto& hit : slice_hits) {
            auto wid = hit->WireID();
            auto key = std::make_pair(wid.Plane, wid.Wire);
            auto it = wire_map.find(key);
            if (it == wire_map.end()) continue;
            const auto& wire = *it->second;

            // Get the signal ROI
            const auto& signal = wire.SignalROI();
            int hit_start_tick = hit->StartTick();
            int hit_end_tick = hit->EndTick();

            for (const auto& range : signal.get_ranges()) {
                int range_start = range.begin_index();
                int range_end = range_start + range.data().size();
                int overlap_start = std::max(range_start, hit_start_tick);
                int overlap_end = std::min(range_end, hit_end_tick);
                if (overlap_start >= overlap_end) continue;

                for (int tick = overlap_start; tick < overlap_end; ++tick) {
                    float adc = range.data()[tick - range_start];
                    if (adc <= 0) continue;

                    double x = _detprop->ConvertTicksToX(tick, wid.Plane);
                    auto wire_pos = _geometry->WireIDToWireGeo(wid).GetCenter();
                    double y = wire_pos.Y();
                    double z = wire_pos.Z();

                    x_pos.push_back(x);
                    y_pos.push_back(y);
                    z_pos.push_back(z);
                    weights.push_back(adc);
                }
            }
        }

        size_t n_points = x_pos.size();
        if (n_points == 0) {
            m_lambda1 = std::numeric_limits<double>::lowest();
            m_lambda2 = std::numeric_limits<double>::lowest();
            m_lambda3 = std::numeric_limits<double>::lowest();
            m_ratio = -1.0;
            return false;
        }

        // Compute weighted mean
        double total_weight = 0.0;
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        for (size_t i = 0; i < n_points; ++i) {
            sum_x += weights[i] * x_pos[i];
            sum_y += weights[i] * y_pos[i];
            sum_z += weights[i] * z_pos[i];
            total_weight += weights[i];
        }

        if (total_weight <= 0) {
            m_lambda1 = std::numeric_limits<double>::lowest();
            m_lambda2 = std::numeric_limits<double>::lowest();
            m_lambda3 = std::numeric_limits<double>::lowest();
            m_ratio = -1.0;
            return false;
        }

        double mean_x = sum_x / total_weight;
        double mean_y = sum_y / total_weight;
        double mean_z = sum_z / total_weight;

        // Compute centered data
        Eigen::MatrixXd centered(n_points, 3);
        for (size_t i = 0; i < n_points; ++i) {
            centered(i, 0) = x_pos[i] - mean_x;
            centered(i, 1) = y_pos[i] - mean_y;
            centered(i, 2) = z_pos[i] - mean_z;
        }

        // Compute weighted covariance matrix
        Eigen::MatrixXd cov(3,3);
        cov.setZero();
        for (size_t i = 0; i < n_points; ++i) {
            Eigen::Vector3d centered_vec = centered.row(i);
            cov += weights[i] * (centered_vec * centered_vec.transpose());
        }
        cov /= total_weight;

        // Compute eigenvalues
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
        if (eig.info() != Eigen::Success) {
            m_lambda1 = std::numeric_limits<double>::lowest();
            m_lambda2 = std::numeric_limits<double>::lowest();
            m_lambda3 = std::numeric_limits<double>::lowest();
            m_ratio = -1.0;
            return false;
        }

        Eigen::VectorXd eigenvalues = eig.eigenvalues().reverse();
        m_lambda1 = eigenvalues(0);
        m_lambda2 = eigenvalues(1);
        m_lambda3 = eigenvalues(2);

        if (m_lambda1 > 0) {
            m_ratio = m_lambda2 / m_lambda1;
        } else {
            m_ratio = -1.0;
        }

        return m_ratio > _threshold;
    }
}

DEFINE_ART_CLASS_TOOL(selection::TopologySelection)
#endif