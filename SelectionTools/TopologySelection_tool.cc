#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ubana/XGBoost/xgboost/c_api.h"
#include "cetlib/search_path.h"

// containment, charge threshold and then this
// probably only use hits and extent here and place a loose cut
namespace selection
{
    class TopologySelection : public SelectionToolBase {
    public:
        explicit TopologySelection(const fhicl::ParameterSet& pset);
        ~TopologySelection();

        void configure(fhicl::ParameterSet const& pset) override;
        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        std::string hitLabel_;
        int minHitsPerPlane_;
        float bdtThreshold_;
        bool trainingMode_;
        std::string trainingDataFile_;
        BoosterHandle bdt_;
        TTree* myTree_;

        // Features
        int nHitsU_, nHitsV_, nHitsW_;
        float totalChargeU_, totalChargeV_, totalChargeW_;
        float wireRangeU_, wireRangeV_, wireRangeW_;
        float timeRangeU_, timeRangeV_, timeRangeW_;
        float bdtScore_;
        int isSignal_; // For training data

        // Helper functions
        std::vector<float> computeFeatures(const art::Event& e);
        float predictBDT(const std::vector<float>& features, BoosterHandle bdt, int nTrees);
        void writeTrainingData(const std::vector<float>& features, int label);
    };

    TopologySelection::TopologySelection(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    TopologySelection::~TopologySelection() {
        if (bdt_) XGBoosterFree(bdt_);
    }

    void TopologySelection::configure(fhicl::ParameterSet const& pset) {
        hitLabel_ = pset.get<std::string>("HitLabel", "gaushit");
        minHitsPerPlane_ = pset.get<int>("MinHitsPerPlane", 5);
        bdtThreshold_ = pset.get<float>("BDTThreshold", 0.5);
        trainingMode_ = pset.get<bool>("TrainingMode", false);
        trainingDataFile_ = pset.get<std::string>("TrainingDataFile", "training_data.txt");

        if (!trainingMode_) {
            int xgtest = -1;
            cet::search_path sp("FW_SEARCH_PATH");
            std::string filename;
            xgtest = XGBoosterCreate(NULL, 0, &bdt_);
            sp.find_file(pset.get<std::string>("BDTModel"), filename);
            xgtest = XGBoosterLoadModel(bdt_, filename.c_str());
            if (xgtest != 0) throw cet::exception("TopologySelection") << "XGBoost error: " << xgtest;
        }
    }

    bool TopologySelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<float> features = computeFeatures(e);

        if (trainingMode_) {
            int label = 0; // Placeholder: 0 = background, 1 = signal (requires simulation truth info)
            writeTrainingData(features, label);
            return false; // Training mode: no selection
        } else {
            bdtScore_ = predictBDT(features, bdt_, 100);
            return bdtScore_ > bdtThreshold_;
        }
    }

    std::vector<float> TopologySelection::computeFeatures(const art::Event& e) {
        art::Handle<std::vector<recob::Hit>> hitHandle;
        e.getByLabel(hitLabel_, hitHandle);
        if (!hitHandle.isValid()) {
            mf::LogWarning("TopologySelection") << "No hits found with label: " << hitLabel_;
            return std::vector<float>(12, 0.0f);
        }

        std::map<int, std::vector<art::Ptr<recob::Hit>>> hitsPerPlane;
        for (size_t i = 0; i < hitHandle->size(); ++i) {
            art::Ptr<recob::Hit> hit(hitHandle, i);
            int plane = hit->WireID().planeID().Plane;
            hitsPerPlane[plane].push_back(hit);
        }

        nHitsU_ = nHitsV_ = nHitsW_ = 0;
        totalChargeU_ = totalChargeV_ = totalChargeW_ = 0.0f;
        wireRangeU_ = wireRangeV_ = wireRangeW_ = 0.0f;
        timeRangeU_ = timeRangeV_ = timeRangeW_ = 0.0f;

        for (int plane = 0; plane < 3; ++plane) {
            auto& hits = hitsPerPlane[plane];
            if (hits.size() < static_cast<size_t>(minHitsPerPlane_)) continue;

            float minWire = std::numeric_limits<float>::max(), maxWire = std::numeric_limits<float>::lowest();
            float minTime = std::numeric_limits<float>::max(), maxTime = std::numeric_limits<float>::lowest();
            float chargeSum = 0.0f;

            for (const auto& hit : hits) {
                float wire = hit->WireID().Wire;
                float time = hit->PeakTime();
                minWire = std::min(minWire, wire);
                maxWire = std::max(maxWire, wire);
                minTime = std::min(minTime, time);
                maxTime = std::max(maxTime, time);
                chargeSum += hit->Integral();
            }

            int nHits = hits.size();
            float wireRange = (nHits > 0) ? maxWire - minWire : 0.0f;
            float timeRange = (nHits > 0) ? maxTime - minTime : 0.0f;

            if (plane == 0) {
                nHitsU_ = nHits;
                totalChargeU_ = chargeSum;
                wireRangeU_ = wireRange;
                timeRangeU_ = timeRange;
            } else if (plane == 1) {
                nHitsV_ = nHits;
                totalChargeV_ = chargeSum;
                wireRangeV_ = wireRange;
                timeRangeV_ = timeRange;
            } else {
                nHitsW_ = nHits;
                totalChargeW_ = chargeSum;
                wireRangeW_ = wireRange;
                timeRangeW_ = timeRange;
            }
        }

        return {
            static_cast<float>(nHitsU_), static_cast<float>(nHitsV_), static_cast<float>(nHitsW_),
            totalChargeU_, totalChargeV_, totalChargeW_,
            wireRangeU_, wireRangeV_, wireRangeW_,
            timeRangeU_, timeRangeV_, timeRangeW_
        };
    }

    float TopologySelection::predictBDT(const std::vector<float>& features, BoosterHandle bdt, int nTrees) {
        int xgtest = -1;
        DMatrixHandle dmat;
        xgtest = XGDMatrixCreateFromMat(const_cast<float*>(features.data()), 1, features.size(), 0, &dmat);
        bst_ulong out_len = 0;
        const float* out_result = NULL;
        xgtest = XGBoosterPredict(bdt, dmat, 0, nTrees, &out_len, &out_result);
        if (xgtest != 0) throw cet::exception("TopologySelection") << "XGBoost predict error: " << xgtest;
        float result = *out_result;
        XGDMatrixFree(dmat);
        return result;
    }

    void TopologySelection::writeTrainingData(const std::vector<float>& features, int label) {
        std::ofstream outfile(trainingDataFile_, std::ios_base::app);
        for (float f : features) {
            outfile << f << " ";
        }
        outfile << label << std::endl;
    }

    void TopologySelection::setBranches(TTree* _tree) {
        _tree->Branch("nHitsU", &nHitsU_, "nHitsU/I");
        _tree->Branch("nHitsV", &nHitsV_, "nHitsV/I");
        _tree->Branch("nHitsW", &nHitsW_, "nHitsW/I");
        _tree->Branch("totalChargeU", &totalChargeU_, "totalChargeU/F");
        _tree->Branch("totalChargeV", &totalChargeV_, "totalChargeV/F");
        _tree->Branch("totalChargeW", &totalChargeW_, "totalChargeW/F");
        _tree->Branch("wireRangeU", &wireRangeU_, "wireRangeU/F");
        _tree->Branch("wireRangeV", &wireRangeV_, "wireRangeV/F");
        _tree->Branch("wireRangeW", &wireRangeW_, "wireRangeW/F");
        _tree->Branch("timeRangeU", &timeRangeU_, "timeRangeU/F");
        _tree->Branch("timeRangeV", &timeRangeV_, "timeRangeV/F");
        _tree->Branch("timeRangeW", &timeRangeW_, "timeRangeW/F");
        _tree->Branch("bdtScore", &bdtScore_, "bdtScore/F");
        _tree->Branch("isSignal", &isSignal_, "isSignal/I");
    }

    void TopologySelection::resetTTree(TTree* _tree) {
        nHitsU_ = nHitsV_ = nHitsW_ = 0;
        totalChargeU_ = totalChargeV_ = totalChargeW_ = 0.0f;
        wireRangeU_ = wireRangeV_ = wireRangeW_ = 0.0f;
        timeRangeU_ = timeRangeV_ = timeRangeW_ = 0.0f;
        bdtScore_ = -1.0f;
        isSignal_ = -1;
    }

    DEFINE_ART_CLASS_TOOL(TopologySelection)
}