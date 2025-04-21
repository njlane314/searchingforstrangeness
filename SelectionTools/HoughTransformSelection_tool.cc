#ifndef SELECTION_HOUGHBDTSELECTION_CXX
#define SELECTION_HOUGHBDTSELECTION_CXX

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include "SelectionToolBase.h"
#include "larreco/RecoAlg/HoughBaseAlg.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Random/RanluxEngine.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubana/XGBoost/xgboost/c_api.h"
#include "cetlib/search_path.h"

namespace selection
{
    class HoughBDTSelection : public SelectionToolBase {
    public:
        explicit HoughBDTSelection(const fhicl::ParameterSet& pset);
        ~HoughBDTSelection();

        void configure(fhicl::ParameterSet const& pset) override;
        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        struct TrackFeatures {
            double length;
            double rmsd;
        };

        std::string hitLabel_;
        std::unique_ptr<HoughBaseAlg> houghAlg_;
        BoosterHandle mainBDT_;
        int minHitsPerLine_;
        float mainBDTThreshold_;
        bool outputTrainingData_;
        std::string trainingDataFile_;
        CLHEP::HepRandomEngine engine_;
        TTree* myTree_;

        // Features
        int nLinesU_, nLinesV_, nLinesW_;
        int totalHitsLinesU_, totalHitsLinesV_, totalHitsLinesW_;
        float wireRangeU_, wireRangeV_, wireRangeW_;
        float timeRangeU_, timeRangeV_, timeRangeW_;
        float pcaRatioU_, pcaRatioV_, pcaRatioW_;
        float avgTrackLengthU_, avgTrackLengthV_, avgTrackLengthW_;
        float avgRMSDU_, avgRMSDV_, avgRMSDW_;
        float mainBDTScore_;
        int isSignal_; // For training data

        // Helper functions
        std::vector<float> computeFeatures(const art::Event& e);
        void extractLineFeatures(const std::vector<protoTrack>& lines, int& nLines, int& totalHitsInLines);
        void computeExtentFeatures(const std::vector<art::Ptr<recob::Hit>>& hits, float& wireRange, float& timeRange, float& pcaRatio);
        float computePCARatio(const std::vector<art::Ptr<recob::Hit>>& hits);
        float predictBDT(const std::vector<float>& features, BoosterHandle bdt, int nTrees);
        void writeTrainingData(const std::vector<float>& features, int label);
        TrackFeatures computeTrackFeatures(const protoTrack& track);
    };

    HoughBDTSelection::HoughBDTSelection(const fhicl::ParameterSet& pset) : engine_(createEngine(0, "Ranlux64")) {
        configure(pset);
    }

    HoughBDTSelection::~HoughBDTSelection() {
        if (mainBDT_) XGBoosterFree(mainBDT_);
    }

    void HoughBDTSelection::configure(fhicl::ParameterSet const& pset) {
        hitLabel_ = pset.get<std::string>("HitLabel", "gaushit");
        minHitsPerLine_ = pset.get<int>("MinHitsPerLine", 10);
        mainBDTThreshold_ = pset.get<float>("MainBDTThreshold", 0.5);
        outputTrainingData_ = pset.get<bool>("OutputTrainingData", false);
        trainingDataFile_ = pset.get<std::string>("TrainingDataFile", "training_data.txt");
        houghAlg_ = std::make_unique<HoughBaseAlg>(pset.get<fhicl::ParameterSet>("HoughAlg"));

        if (!outputTrainingData_) {
            int xgtest = -1;
            cet::search_path sp("FW_SEARCH_PATH");
            std::string filename;
            xgtest = XGBoosterCreate(NULL, 0, &mainBDT_);
            sp.find_file(pset.get<std::string>("MainBDTModel"), filename);
            xgtest = XGBoosterLoadModel(mainBDT_, filename.c_str());
            if (xgtest != 0) throw cet::exception("HoughBDTSelection") << "XGBoost error: " << xgtest;
        }
    }

    bool HoughBDTSelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<float> features = computeFeatures(e);

        if (outputTrainingData_) {
            int label = 0; // Placeholder: 0 = background, 1 = signal (requires simulation truth info)
            writeTrainingData(features, label);
            return false; // Training mode: no selection
        } else {
            mainBDTScore_ = predictBDT(features, mainBDT_, 100);
            return mainBDTScore_ > mainBDTThreshold_;
        }
    }

    std::vector<float> HoughBDTSelection::computeFeatures(const art::Event& e) {
        art::Handle<std::vector<recob::Hit>> hitHandle;
        e.getByLabel(hitLabel_, hitHandle);
        if (!hitHandle.isValid()) {
            mf::LogWarning("HoughBDTSelection") << "No hits found with label: " << hitLabel_;
            return std::vector<float>(21, 0.0f);
        }

        std::map<int, std::vector<art::Ptr<recob::Hit>>> hitsPerPlane;
        for (size_t i = 0; i < hitHandle->size(); ++i) {
            art::Ptr<recob::Hit> hit(hitHandle, i);
            int plane = hit->WireID().planeID().Plane;
            hitsPerPlane[plane].push_back(hit);
        }

        nLinesU_ = nLinesV_ = nLinesW_ = 0;
        totalHitsLinesU_ = totalHitsLinesV_ = totalHitsLinesW_ = 0;
        wireRangeU_ = wireRangeV_ = wireRangeW_ = 0.0f;
        timeRangeU_ = timeRangeV_ = timeRangeW_ = 0.0f;
        pcaRatioU_ = pcaRatioV_ = pcaRatioW_ = 0.0f;
        avgTrackLengthU_ = avgTrackLengthV_ = avgTrackLengthW_ = 0.0f;
        avgRMSDU_ = avgRMSDV_ = avgRMSDW_ = 0.0f;

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);

        for (int plane = 0; plane < 3; ++plane) {
            auto& hits = hitsPerPlane[plane];
            if (hits.empty()) continue;

            std::vector<unsigned int> pointIdToClusterId(hits.size(), 0);
            unsigned int clusterId = 0, nClusters = 0;
            std::vector<protoTrack> linesFound;
            houghAlg_->Transform(clockData, detProp, hits, &engine_, &pointIdToClusterId, clusterId, &nClusters, &linesFound);

            int nLines = 0, totalHitsInLines = 0;
            extractLineFeatures(linesFound, nLines, totalHitsInLines);

            std::vector<TrackFeatures> sigTrackFeatures;
            for (const auto& line : linesFound) {
                if (line.hits.size() >= static_cast<size_t>(minHitsPerLine_)) {
                    TrackFeatures tf = computeTrackFeatures(line);
                    sigTrackFeatures.push_back(tf);
                }
            }
            if (!sigTrackFeatures.empty()) {
                double sumLength = 0.0, sumRMSD = 0.0;
                for (const auto& tf : sigTrackFeatures) {
                    sumLength += tf.length;
                    sumRMSD += tf.rmsd;
                }
                float avgLength = static_cast<float>(sumLength / sigTrackFeatures.size());
                float avgRMSD = static_cast<float>(sumRMSD / sigTrackFeatures.size());
                if (plane == 0) {
                    avgTrackLengthU_ = avgLength;
                    avgRMSDU_ = avgRMSD;
                } else if (plane == 1) {
                    avgTrackLengthV_ = avgLength;
                    avgRMSDV_ = avgRMSD;
                } else {
                    avgTrackLengthW_ = avgLength;
                    avgRMSDW_ = avgRMSD;
                }
            }

            if (plane == 0) { nLinesU_ = nLines; totalHitsLinesU_ = totalHitsInLines; }
            else if (plane == 1) { nLinesV_ = nLines; totalHitsLinesV_ = totalHitsInLines; }
            else { nLinesW_ = nLines; totalHitsLinesW_ = totalHitsInLines; }

            float wireRange, timeRange, pcaRatio;
            computeExtentFeatures(hits, wireRange, timeRange, pcaRatio);
            if (plane == 0) { wireRangeU_ = wireRange; timeRangeU_ = timeRange; pcaRatioU_ = pcaRatio; }
            else if (plane == 1) { wireRangeV_ = wireRange; timeRangeV_ = timeRange; pcaRatioV_ = pcaRatio; }
            else { wireRangeW_ = wireRange; timeRangeW_ = timeRange; pcaRatioW_ = pcaRatio; }
        }

        return {
            static_cast<float>(nLinesU_), static_cast<float>(nLinesV_), static_cast<float>(nLinesW_),
            static_cast<float>(totalHitsLinesU_), static_cast<float>(totalHitsLinesV_), static_cast<float>(totalHitsLinesW_),
            wireRangeU_, wireRangeV_, wireRangeW_,
            timeRangeU_, timeRangeV_, timeRangeW_,
            pcaRatioU_, pcaRatioV_, pcaRatioW_,
            avgTrackLengthU_, avgTrackLengthV_, avgTrackLengthW_,
            avgRMSDU_, avgRMSDV_, avgRMSDW_
        };
    }

    void HoughBDTSelection::extractLineFeatures(const std::vector<protoTrack>& lines, int& nLines, int& totalHitsInLines) {
        nLines = 0;
        totalHitsInLines = 0;
        for (const auto& line : lines) {
            if (line.hits.size() >= static_cast<size_t>(minHitsPerLine_)) {
                nLines++;
                totalHitsInLines += line.hits.size();
            }
        }
    }

    void HoughBDTSelection::computeExtentFeatures(const std::vector<art::Ptr<recob::Hit>>& hits, float& wireRange, float& timeRange, float& pcaRatio) {
        if (hits.empty()) {
            wireRange = timeRange = pcaRatio = 0.0f;
            return;
        }

        float minWire = std::numeric_limits<float>::max(), maxWire = std::numeric_limits<float>::lowest();
        float minTime = std::numeric_limits<float>::max(), maxTime = std::numeric_limits<float>::lowest();
        for (const auto& hit : hits) {
            float wire = hit->WireID().Wire;
            float time = hit->PeakTime();
            minWire = std::min(minWire, wire);
            maxWire = std::max(maxWire, wire);
            minTime = std::min(minTime, time);
            maxTime = std::max(maxTime, time);
        }
        wireRange = maxWire - minWire;
        timeRange = maxTime - minTime;
        pcaRatio = computePCARatio(hits);
    }

    float HoughBDTSelection::computePCARatio(const std::vector<art::Ptr<recob::Hit>>& hits) {
        if (hits.size() < 2) return 0.0f;

        double meanWire = 0.0, meanTime = 0.0;
        for (const auto& hit : hits) {
            meanWire += hit->WireID().Wire;
            meanTime += hit->PeakTime();
        }
        meanWire /= hits.size();
        meanTime /= hits.size();

        double covWW = 0.0, covWT = 0.0, covTT = 0.0;
        for (const auto& hit : hits) {
            double w = hit->WireID().Wire - meanWire;
            double t = hit->PeakTime() - meanTime;
            covWW += w * w;
            covWT += w * t;
            covTT += t * t;
        }
        covWW /= hits.size();
        covWT /= hits.size();
        covTT /= hits.size();

        double trace = covWW + covTT;
        double det = covWW * covTT - covWT * covWT;
        double discriminant = trace * trace - 4 * det;
        if (discriminant < 0) return 0.0f;
        double sqrtDiscriminant = std::sqrt(discriminant);
        double lambda1 = (trace + sqrtDiscriminant) / 2;
        double lambda2 = (trace - sqrtDiscriminant) / 2;
        if (lambda1 <= 0 || lambda2 <= 0) return 0.0f;
        return lambda2 / lambda1;
    }

    float HoughBDTSelection::predictBDT(const std::vector<float>& features, BoosterHandle bdt, int nTrees) {
        int xgtest = -1;
        DMatrixHandle dmat;
        xgtest = XGDMatrixCreateFromMat(const_cast<float*>(features.data()), 1, features.size(), 0, &dmat);
        bst_ulong out_len = 0;
        const float* out_result = NULL;
        xgtest = XGBoosterPredict(bdt, dmat, 0, nTrees, &out_len, &out_result);
        if (xgtest != 0) throw cet::exception("HoughBDTSelection") << "XGBoost predict error: " << xgtest;
        float result = *out_result;
        XGDMatrixFree(dmat);
        return result;
    }

    void HoughBDTSelection::writeTrainingData(const std::vector<float>& features, int label) {
        std::ofstream outfile(trainingDataFile_, std::ios_base::app);
        for (float f : features) {
            outfile << f << " ";
        }
        outfile << label << std::endl;
    }

    HoughBDTSelection::TrackFeatures HoughBDTSelection::computeTrackFeatures(const protoTrack& track) {
        TrackFeatures features;

        // Compute length in physical space (cm)
        double dw = track.pMax0 - track.pMin0; // Wire direction distance
        double dt = track.pMax1 - track.pMin1; // Drift direction distance
        features.length = std::sqrt(dw * dw + dt * dt);

        // Compute RMSD in (wire, time) space
        if (!track.hits.empty()) {
            double sumSqDev = 0.0;
            for (const auto& hit : track.hits) {
                double w = hit->WireID().Wire;
                double t = hit->PeakTime();
                double t_exp = track.clusterSlope * w + track.clusterIntercept;
                double dev = t - t_exp;
                sumSqDev += dev * dev;
            }
            features.rmsd = std::sqrt(sumSqDev / track.hits.size());
        } else {
            features.rmsd = 0.0;
        }

        return features;
    }

    void HoughBDTSelection::setBranches(TTree* _tree) {
        _tree->Branch("nLinesU", &nLinesU_, "nLinesU/I");
        _tree->Branch("nLinesV", &nLinesV_, "nLinesV/I");
        _tree->Branch("nLinesW", &nLinesW_, "nLinesW/I");
        _tree->Branch("totalHitsLinesU", &totalHitsLinesU_, "totalHitsLinesU/I");
        _tree->Branch("totalHitsLinesV", &totalHitsLinesV_, "totalHitsLinesV/I");
        _tree->Branch("totalHitsLinesW", &totalHitsLinesW_, "totalHitsLinesW/I");
        _tree->Branch("wireRangeU", &wireRangeU_, "wireRangeU/F");
        _tree->Branch("wireRangeV", &wireRangeV_, "wireRangeV/F");
        _tree->Branch("wireRangeW", &wireRangeW_, "wireRangeW/F");
        _tree->Branch("timeRangeU", &timeRangeU_, "timeRangeU/F");
        _tree->Branch("timeRangeV", &timeRangeV_, "timeRangeV/F");
        _tree->Branch("timeRangeW", &timeRangeW_, "timeRangeW/F");
        _tree->Branch("pcaRatioU", &pcaRatioU_, "pcaRatioU/F");
        _tree->Branch("pcaRatioV", &pcaRatioV_, "pcaRatioV/F");
        _tree->Branch("pcaRatioW", &pcaRatioW_, "pcaRatioW/F");
        _tree->Branch("avgTrackLengthU", &avgTrackLengthU_, "avgTrackLengthU/F");
        _tree->Branch("avgTrackLengthV", &avgTrackLengthV_, "avgTrackLengthV/F");
        _tree->Branch("avgTrackLengthW", &avgTrackLengthW_, "avgTrackLengthW/F");
        _tree->Branch("avgRMSDU", &avgRMSDU_, "avgRMSDU/F");
        _tree->Branch("avgRMSDV", &avgRMSDV_, "avgRMSDV/F");
        _tree->Branch("avgRMSDW", &avgRMSDW_, "avgRMSDW/F");
        _tree->Branch("mainBDTScore", &mainBDTScore_, "mainBDTScore/F");
        _tree->Branch("isSignal", &isSignal_, "isSignal/I");
    }

    void HoughBDTSelection::resetTTree(TTree* _tree) {
        nLinesU_ = nLinesV_ = nLinesW_ = 0;
        totalHitsLinesU_ = totalHitsLinesV_ = totalHitsLinesW_ = 0;
        wireRangeU_ = wireRangeV_ = wireRangeW_ = 0.0f;
        timeRangeU_ = timeRangeV_ = timeRangeW_ = 0.0f;
        pcaRatioU_ = pcaRatioV_ = pcaRatioW_ = 0.0f;
        avgTrackLengthU_ = avgTrackLengthV_ = avgTrackLengthW_ = 0.0f;
        avgRMSDU_ = avgRMSDV_ = avgRMSDW_ = 0.0f;
        mainBDTScore_ = -1.0f;
        isSignal_ = -1;
    }

    DEFINE_ART_CLASS_TOOL(HoughBDTSelection)
}

#endif


/* 
import xgboost as xgb
import numpy as np
from sklearn.model_selection import train_test_split

# Load training data from the tool's output
data = np.loadtxt('training_data.txt')
X = data[:, :-1]  # Features
y = data[:, -1]   # Labels (0 = background, 1 = signal)

# Split into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train the BDT
dtrain = xgb.DMatrix(X_train, label=y_train)
params = {'max_depth': 3, 'eta': 0.1, 'objective': 'binary:logistic'}
bst = xgb.train(params, dtrain, num_boost_round=100)
bst.save_model('main_bdt.model')
*/