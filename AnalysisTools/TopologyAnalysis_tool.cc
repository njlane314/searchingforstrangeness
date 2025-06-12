#ifndef ANALYSIS_TOPOLOGYANALYSIS_CXX
#define ANALYSIS_TOPOLOGYANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include <canvas/Persistency/Common/FindManyP.h>
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <vector>
#include <limits>

namespace analysis 
{
    class TopologyAnalysis : public AnalysisToolBase {
    public:
        explicit TopologyAnalysis(const fhicl::ParameterSet& pset) { configure(pset); }
        ~TopologyAnalysis() = default;

        void configure(fhicl::ParameterSet const& pset) override;
        void analyseEvent(art::Event const& e, bool fData) override;
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag fCLSproducer;

        int m_nhit_u, m_nhits_v, m_nhits_w;
        float m_charge_u, m_charge_v, m_charge_w;
        float m_wirerange_u, m_wirerange_v, m_wirerange_w;
        float m_timerange_u, m_timerange_v, m_timerange_w;

        void computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
    };

    void TopologyAnalysis::configure(fhicl::ParameterSet const& pset) {
        fCLSproducer = pset.get<art::InputTag>("CLSproducer", "pandora");
    }

    void TopologyAnalysis::analyseEvent(art::Event const& e, bool fData) {}

    void TopologyAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) {
        std::cout << "TopologyAnalysis::analyseSlice: " << slice_pfp_v.size() << " PFPs in slice." << std::endl;
        this->computeFeatures(e, slice_pfp_v);
    }

    void TopologyAnalysis::computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
        auto clusterHandle = e.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
        art::FindManyP<recob::Hit> clusterHits(clusterHandle, e, fCLSproducer);
        
        std::map<int, std::vector<art::Ptr<recob::Hit>>> sliceHitsPerPlane;
        for (const auto& pfp : slice_pfp_v) {
            auto clusters = pfp.get<recob::Cluster>();
            for (const auto& cluster : clusters) {
                const std::vector<art::Ptr<recob::Hit>>& hits = clusterHits.at(cluster.key());
                for (const auto& hit : hits) {
                    int plane = hit->WireID().planeID().Plane;
                    if (plane < 0 || plane > 2) continue;  
                    sliceHitsPerPlane[plane].push_back(hit);
                }
            }
        }

        for (int plane = 0; plane < 3; ++plane) {
            auto& hits = sliceHitsPerPlane[plane];
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
                m_nhit_u = nhits; m_charge_u = charge_sum; m_wirerange_u = wirerange; m_timerange_u = timerange;
            } else if (plane == 1) {
                m_nhits_v = nhits; m_charge_v = charge_sum; m_wirerange_v = wirerange; m_timerange_v = timerange;
            } else {
                m_nhits_w = nhits; m_charge_w = charge_sum; m_wirerange_w = wirerange; m_timerange_w = timerange;
            }
        }
    }

    void TopologyAnalysis::setBranches(TTree* tree) {
        tree->Branch("nhits_u", &m_nhit_u, "nhits_u/I");
        tree->Branch("nhit_v", &m_nhits_v, "nhit_v/I");
        tree->Branch("nhit_w", &m_nhits_w, "nhit_w/I");
        tree->Branch("charge_u", &m_charge_u, "charge_u/F");
        tree->Branch("charge_v", &m_charge_v, "charge_v/F");
        tree->Branch("charge_w", &m_charge_w, "charge_w/F");
        tree->Branch("wirerange_u", &m_wirerange_u, "wirerange_u/F");
        tree->Branch("wirerange_v", &m_wirerange_v, "wirerange_v/F");
        tree->Branch("wirerange_w", &m_wirerange_w, "wirerange_w/F");
        tree->Branch("timerange_u", &m_timerange_u, "timerange_u/F");
        tree->Branch("timerange_v", &m_timerange_v, "timerange_v/F");
        tree->Branch("timerange_w", &m_timerange_w, "timerange_w/F");
    }

    void TopologyAnalysis::resetTTree(TTree* tree) {
        m_nhit_u = m_nhits_v = m_nhits_w = 0;
        m_charge_u = m_charge_v = m_charge_w = 0.0f;
        m_wirerange_u = m_wirerange_v = m_wirerange_w = 0.0f;
        m_timerange_u = m_timerange_v = m_timerange_w = 0.0f;
    }

    DEFINE_ART_CLASS_TOOL(TopologyAnalysis)
}

#endif