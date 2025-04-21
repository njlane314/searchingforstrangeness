#ifndef SELECTION_RESNETLOGISTICSELECTION_CXX
#define SELECTION_RESNETLOGISTICSELECTION_CXX

#include <iostream>
#include <vector>
#include <fstream>
#include <torch/script.h>
#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ImageProcessor.h"
#include "cetlib/search_path.h"

namespace selection
{
    class ContrastiveSelection : public SelectionToolBase {
    public:
        explicit ContrastiveSelection(const fhicl::ParameterSet& pset);
        ~ContrastiveSelection() override = default;

        void configure(const fhicl::ParameterSet& pset) override;
        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        std::string image_label_;
        bool training_mode_;
        std::string training_data_file_;
        float selection_threshold_;
        torch::jit::script::Module resnet_model_;
        torch::jit::script::Module logistic_model_;
        TTree* my_tree_;

        // Features and outputs
        std::vector<float> plane_scores_; // Scores for U, V, W planes
        float final_score_;              // Final logistic regression score
        int is_signal_;                  // Truth label for training

        // Helper functions
        std::vector<float> getPlaneScores(const art::Event& e);
        float getFinalScore(const std::vector<float>& plane_scores);
        void writeTrainingData(const std::vector<float>& plane_scores, int label);
    };

    ContrastiveSelection::ContrastiveSelection(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    void ContrastiveSelection::configure(const fhicl::ParameterSet& pset) {
        image_label_ = pset.get<std::string>("ImageLabel", "inputImages");
        training_mode_ = pset.get<bool>("TrainingMode", false);
        training_data_file_ = pset.get<std::string>("TrainingDataFile", "training_data.txt");
        selection_threshold_ = pset.get<float>("SelectionThreshold", 0.5);

        // Load pre-trained models
        cet::search_path sp("FW_SEARCH_PATH");
        std::string resnet_path, logistic_path;
        sp.find_file(pset.get<std::string>("ResNetModel"), resnet_path);
        sp.find_file(pset.get<std::string>("LogisticModel"), logistic_path);

        try {
            resnet_model_ = torch::jit::load(resnet_path);
            resnet_model_.eval();
            if (!training_mode_) {
                logistic_model_ = torch::jit::load(logistic_path);
                logistic_model_.eval();
            }
        } catch (const c10::Error& e) {
            throw cet::exception("ContrastiveSelection") << "Failed to load models: " << e.what();
        }
    }

    bool ContrastiveSelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        auto images = e.getValidHandle<std::vector<image::Image>>(image_label_);
        if (images->size() != 3) {
            mf::LogWarning("ContrastiveSelection") << "Expected 3 images, got " << images->size();
            return false;
        }

        // Extract scores for each plane
        plane_scores_ = getPlaneScores(e);

        if (training_mode_) {
            // Placeholder: Replace with actual truth label from simulation
            int label = 0; // 0 = background, 1 = signal (requires MCTruth integration)
            writeTrainingData(plane_scores_, label);
            return false; // No selection in training mode
        } else {
            // Compute final score using logistic regression
            final_score_ = getFinalScore(plane_scores_);
            return final_score_ > selection_threshold_;
        }
    }

    std::vector<float> ContrastiveSelection::getPlaneScores(const art::Event& e) {
        auto images = e.getValidHandle<std::vector<image::Image>>(image_label_);
        std::vector<float> scores;
        scores.reserve(3);

        for (const auto& img : *images) {
            // Convert image to tensor
            auto data = img.data();
            torch::Tensor tensor = torch::from_blob(data.data(), {1, 1, static_cast<long>(img.height()), static_cast<long>(img.width())})
                                      .to(torch::kFloat32);
            // Normalize tensor (assuming ResNet expects normalized input)
            tensor = (tensor - tensor.mean()) / tensor.std();

            // Pass through ResNet
            torch::NoGradGuard no_grad;
            torch::Tensor output = resnet_model_.forward({tensor}).toTensor();
            // Assume output is [1, num_classes], take signal class probability (class 1)
            float score = torch::softmax(output, 1)[0][1].item<float>();
            scores.push_back(score);
        }

        return scores;
    }

    float ContrastiveSelection::getFinalScore(const std::vector<float>& plane_scores) {
        // Convert scores to tensor
        torch::Tensor input_tensor = torch::from_blob(const_cast<float*>(plane_scores.data()), {1, 3}).to(torch::kFloat32);
        // Pass through logistic regression
        torch::NoGradGuard no_grad;
        torch::Tensor output = logistic_model_.forward({input_tensor}).toTensor();
        // Assume output is [1, 1] (sigmoid output for binary classification)
        return output[0][0].item<float>();
    }

    void ContrastiveSelection::writeTrainingData(const std::vector<float>& plane_scores, int label) {
        std::ofstream outfile(training_data_file_, std::ios_base::app);
        if (!outfile.is_open()) {
            mf::LogError("ContrastiveSelection") << "Failed to open training data file: " << training_data_file_;
            return;
        }
        for (float score : plane_scores) {
            outfile << score << " ";
        }
        outfile << label << std::endl;
    }

    void ContrastiveSelection::setBranches(TTree* _tree) {
        my_tree_ = _tree;
        my_tree_->Branch("plane_scores", &plane_scores_);
        my_tree_->Branch("final_score", &final_score_, "final_score/F");
        my_tree_->Branch("is_signal", &is_signal_, "is_signal/I");
    }

    void ContrastiveSelection::resetTTree(TTree* _tree) {
        plane_scores_.clear();
        final_score_ = -1.0f;
        is_signal_ = -1;
    }

    DEFINE_ART_CLASS_TOOL(ContrastiveSelection)
}

#endif