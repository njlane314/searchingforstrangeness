#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ubana/XGBoost/xgboost/c_api.h"
#include "cetlib/search_path.h"
#include <vector>
#include <limits>

namespace selection {
    class TopologySelection : public SelectionToolBase {
    public:
        explicit TopologySelection(const fhicl::ParameterSet& pset) { configure(pset); }
        ~TopologySelection() { if (m_bdt) XGBoosterFree(m_bdt); }

        void configure(fhicl::ParameterSet const& pset) override;
        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) override;
        void analyzeSlice(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        float m_bdtThreshold;
        bool m_inferenceMode;
        BoosterHandle m_bdt = nullptr;
        float m_bdtScore = -1.0f;

        int m_nHitsU, m_nHitsV, m_nHitsW;
        float m_totalChargeU, m_totalChargeV, m_totalChargeW;
        float m_wireRangeU, m_wireRangeV, m_wireRangeW;
        float m_timeRangeU, m_timeRangeV, m_timeRangeW;

        void computeFeatures(const std::vector<common::ProxyPfpElem_t>& slice_pfp_v);
        float predictBDT(const std::vector<float>& features);
    };

    void TopologySelection::configure(fhicl::ParameterSet const& pset) {
        m_bdtThreshold = pfp.get<float>("BDTThreshold", 0.5);
        m_inferenceMode = pset.get<bool>("InferenceMode", false);

        if (m_inferenceMode) {
            cet::search_path sp("FW_SEARCH_PATH");
            std::string filename;
            if (!sp.find_file(pset.get<std::string>("BDTModel"), filename))
                throw cet::exception("TopologySelection") << "BDT model file not found.";
            if (XGBoosterCreate(NULL, 0, &m_bdt) != 0 || XGBoosterLoadModel(m_bdt, filename.c_str()) != 0)
                throw cet::exception("TopologySelection") << "Failed to initialize XGBoost booster.";
        }
    }

    void TopologySelection::analyzeSlice(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
        computeFeatures(slice_pfp_v);
    }

    bool TopologySelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        if (!m_inferenceMode) {
            return false;
        }

        std::vector<float> features = {
            static_cast<float>(m_nHitsU), static_cast<float>(m_nHitsV), static_cast<float>(m_nHitsW),
            m_totalChargeU, m_totalChargeV, m_totalChargeW,
            m_wireRangeU, m_wireRangeV, m_wireRangeW,
            m_timeRangeU, m_timeRangeV, m_timeRangeW
        };
        m_bdtScore = predictBDT(features);
        return m_bdtScore > m_bdtThreshold; 
    }

    void TopologySelection::computeFeatures(const std::vector<common::ProxyPfpElem_t>& slice_pfp_v) {
        std::map<int, std::vector<art::Ptr<recob::Hit>>> sliceHitsPerPlane;
        for (const auto& pfp : slice_pfp_v) {
            auto hits = pfp.get<recob::Hit>();
            for (const auto& hit : hits) {
                int plane = hit->WireID().planeID().Plane;
                sliceHitsPerPlane[plane].push_back(hit);
            }
        }

        for (int plane = 0; plane < 3; ++plane) {
            auto& hits = sliceHitsPerPlane[plane];
            float minWire = std::numeric_limits<float>::max(), maxWire = std::numeric_limits<float>::lowest();
            float minTime = std::numeric_limits<float>::max(), maxTime = std::numeric_limits<float>::lowest();
            float chargeSum = 0.0f;

            for (const auto& hit : hits) {
                float wire = hit->WireID().Wire, time = hit->PeakTime();
                minWire = std::min(minWire, wire); maxWire = std::max(maxWire, wire);
                minTime = std::min(minTime, time); maxTime = std::max(maxTime, time);
                chargeSum += hit->Integral();
            }

            int nHits = hits.size();
            float wireRange = nHits ? maxWire - minWire : 0.0f;
            float timeRange = nHits ? maxTime - minTime : 0.0f;

            if (plane == 0) {
                m_nHitsU = nHits; m_totalChargeU = chargeSum; m_wireRangeU = wireRange; m_timeRangeU = timeRange;
            } else if (plane == 1) {
                m_nHitsV = nHits; m_totalChargeV = chargeSum; m_wireRangeV = wireRange; m_timeRangeV = timeRange;
            } else {
                m_nHitsW = nHits; m_totalChargeW = chargeSum; m_wireRangeW = wireRange; m_timeRangeW = timeRange;
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
        tree->Branch("nHitsU", &m_nHitsU, "nHitsU/I");
        tree->Branch("nHitsV", &m_nHitsV, "nHitsV/I");
        tree->Branch("nHitsW", &m_nHitsW, "nHitsW/I");
        tree->Branch("totalChargeU", &m_totalChargeU, "totalChargeU/F");
        tree->Branch("totalChargeV", &m_totalChargeV, "totalChargeV/F");
        tree->Branch("totalChargeW", &m_totalChargeW, "totalChargeW/F");
        tree->Branch("wireRangeU", &m_wireRangeU, "wireRangeU/F");
        tree->Branch("wireRangeV", &m_wireRangeV, "wireRangeV/F");
        tree->Branch("wireRangeW", &m_wireRangeW, "wireRangeW/F");
        tree->Branch("timeRangeU", &m_timeRangeU, "timeRangeU/F");
        tree->Branch("timeRangeV", &m_timeRangeV, "timeRangeV/F");
        tree->Branch("timeRangeW", &m_timeRangeW, "timeRangeW/F");
        if (m_inferenceMode) {
            tree->Branch("bdtScore", &m_bdtScore, "bdtScore/F");
        }
    }

    void TopologySelection::resetTTree(TTree* tree) {
        m_nHitsU = m_nHitsV = m_nHitsW = 0;
        m_totalChargeU = m_totalChargeV = m_totalChargeW = 0.0f;
        m_wireRangeU = m_wireRangeV = m_wireRangeW = 0.0f;
        m_timeRangeU = m_timeRangeV = m_timeRangeW = 0.0f;
        m_bdtScore = -1.0f;
    }
}

DEFINE_ART_CLASS_TOOL(selection::TopologySelection)