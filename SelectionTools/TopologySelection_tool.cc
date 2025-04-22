#ifndef SELECTION_TOPOLOGY_CXX
#define SELECTION_TOPOLOGY_CXX

#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include <canvas/Persistency/Common/FindManyP.h>
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ubana/XGBoost/xgboost/c_api.h"
#include "cetlib/search_path.h"
#include <vector>
#include <limits>

namespace selection 
{
    class TopologySelection : public SelectionToolBase {
    public:
        explicit TopologySelection(const fhicl::ParameterSet& pset) { configure(pset); }
        ~TopologySelection() { if (m_bdt) XGBoosterFree(m_bdt); }

        void configure(fhicl::ParameterSet const& pset);
        bool selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image>& calo_images, 
                                const std::vector<image::Image>& reco_images, 
                                const std::vector<image::Image>& label_images);
        void analyseSlice(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
        void setBranches(TTree* tree);
        void resetTTree(TTree* tree);

    private:
        art::InputTag fCLSproducer;
        float m_bdtThreshold;
        bool m_inferenceMode;
        BoosterHandle m_bdt = nullptr;
        float m_bdt_score = -1.0f;

        int m_nhit_u, m_nhits_v, m_nhits_w;
        float m_charge_u, m_charge_v, m_charge_w;
        float m_wirerange_u, m_wirerange_v, m_wirerange_w;
        float m_timerange_u, m_timerange_v, m_timerange_w;

        void computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
        float predictBDT(const std::vector<float>& features);
    };

    void TopologySelection::configure(fhicl::ParameterSet const& pset) {
        m_bdtThreshold = pset.get<float>("BDTThreshold", 0.5);
        m_inferenceMode = pset.get<bool>("InferenceMode", false);
        fCLSproducer = pset.get<art::InputTag>("CLSproducer");

        if (m_inferenceMode) {
            cet::search_path sp("FW_SEARCH_PATH");
            std::string filename;
            if (!sp.find_file(pset.get<std::string>("BDTModel"), filename))
                throw cet::exception("TopologySelection") << "BDT model file not found.";
            if (XGBoosterCreate(NULL, 0, &m_bdt) != 0 || XGBoosterLoadModel(m_bdt, filename.c_str()) != 0)
                throw cet::exception("TopologySelection") << "Failed to initialize XGBoost booster.";
        }
    }

    void TopologySelection::analyseSlice(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
        computeFeatures(e, slice_pfp_v);
    }

    bool TopologySelection::selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image>& calo_images, 
                                const std::vector<image::Image>& reco_images, 
                                const std::vector<image::Image>& label_images) {
        if (!m_inferenceMode) {
            return false;
        }

        std::vector<float> features = {
            static_cast<float>(m_nhit_u), static_cast<float>(m_nhits_v), static_cast<float>(m_nhits_w),
            m_charge_u, m_charge_v, m_charge_w,
            m_wirerange_u, m_wirerange_v, m_wirerange_w,
            m_timerange_u, m_timerange_v, m_timerange_w
        };
        m_bdt_score = predictBDT(features);
        return m_bdt_score > m_bdtThreshold; 
    }

    void TopologySelection::computeFeatures(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
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

    float TopologySelection::predictBDT(const std::vector<float>& features) {
        DMatrixHandle dmat;
        if (XGDMatrixCreateFromMat(const_cast<float*>(features.data()), 1, features.size(), 0, &dmat) != 0)
            throw cet::exception("TopologySelection") << "Failed to create DMatrix.";
        bst_ulong out_len;
        const float* out_result;
        if (XGBoosterPredict(m_bdt, dmat, 0, 0, &out_len, &out_result) != 0) {
            XGDMatrixFree(dmat);
            throw cet::exception("TopologySelection") << "XGBoost prediction failed.";
        }
        float result = *out_result;
        XGDMatrixFree(dmat);
        return result;
    }

    void TopologySelection::setBranches(TTree* tree) {
        tree->Branch("nhits_u", &m_nhit_u, "nhits_u/I");
        tree->Branch("nhit_v", &m_nhits_v, "nhit_v/I");
        tree->Branch("nhit_w", &m_nhits_w, "nhit_w/I");
        tree->Branch("charge_u", &m_charge_u, "charge_u/F");
        tree->Branch("charge_v", &m_charge_v, "charge_v/F");
        tree->Branch("chare_w", &m_charge_w, "chare_w/F");
        tree->Branch("wirerange_u", &m_wirerange_u, "wirerange_u/F");
        tree->Branch("wirerange_v", &m_wirerange_v, "wirerange_v/F");
        tree->Branch("wirerange_w", &m_wirerange_w, "wirerange_w/F");
        tree->Branch("timerange_u", &m_timerange_u, "timerange_u/F");
        tree->Branch("timerange_v", &m_timerange_v, "timerange_v/F");
        tree->Branch("timerange_w", &m_timerange_w, "timerange_w/F");
        if (m_inferenceMode) {
            tree->Branch("bdt_score", &m_bdt_score, "bdt_score/F");
        }
    }

    void TopologySelection::resetTTree(TTree* tree) {
        m_nhit_u = m_nhits_v = m_nhits_w = 0;
        m_charge_u = m_charge_v = m_charge_w = 0.0f;
        m_wirerange_u = m_wirerange_v = m_wirerange_w = 0.0f;
        m_timerange_u = m_timerange_v = m_timerange_w = 0.0f;
        m_bdt_score = -1.0f;
    }

    DEFINE_ART_CLASS_TOOL(TopologySelection)
}

#endif