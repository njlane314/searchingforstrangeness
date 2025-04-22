#ifndef SELECTION_CONTRASTIVE_CXX
#define SELECTION_CONTRASTIVE_CXX

#include <iostream>
#include <vector>
#include <torch/script.h>
#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "../CommonDefs/Image.h"
#include "cetlib/search_path.h"

namespace selection
{
    class ContrastiveSelection : public SelectionToolBase {
    public:
        explicit ContrastiveSelection(const fhicl::ParameterSet& pset);
        ~ContrastiveSelection() override = default;
        void configure(const fhicl::ParameterSet& pset);
        bool selectEvent(art::Event const& e, 
                        const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                        const std::vector<image::Image>& calo_images, 
                        const std::vector<image::Image>& reco_images, 
                        const std::vector<image::Image>& label_images);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);
    private:
        std::string image_label_;
        float selection_threshold_;
        torch::jit::script::Module resnet_model_;
        std::vector<float> plane_scores_;
        float final_score_;
        std::vector<float> getPlaneScores(const std::vector<image::Image>& images);
    };

    ContrastiveSelection::ContrastiveSelection(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    void ContrastiveSelection::configure(const fhicl::ParameterSet& pset) {
        image_label_ = pset.get<std::string>("ImageLabel", "inputImages");
        selection_threshold_ = pset.get<float>("SelectionThreshold", 0.5);
        cet::search_path sp("FW_SEARCH_PATH");
        std::string resnet_path;
        sp.find_file(pset.get<std::string>("ResNetModel"), resnet_path);
        try {
            resnet_model_ = torch::jit::load(resnet_path);
            resnet_model_.eval();
        } catch (const c10::Error& e) {
            throw cet::exception("ContrastiveSelection") << "Failed to load ResNet model: " << e.what();
        }
    }

    bool ContrastiveSelection::selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image>& calo_images, 
                                const std::vector<image::Image>& reco_images, 
                                const std::vector<image::Image>& label_images) {
        auto images = e.getValidHandle<std::vector<image::Image>>(image_label_);
        if (images->size() != 3) {
            mf::LogWarning("ContrastiveSelection") << "Expected 3 images, got " << images->size();
            return false;
        }
        plane_scores_ = getPlaneScores(*images);
        if (plane_scores_.size() != 3) {
            return false;
        }
        final_score_ = (plane_scores_[0] + plane_scores_[1] + plane_scores_[2]) / 3.0;
        return final_score_ > selection_threshold_;
    }

    std::vector<float> ContrastiveSelection::getPlaneScores(const std::vector<image::Image>& images) {
        std::vector<float> scores;
        scores.reserve(3);
        for (const auto& img : images) {
            auto data = img.data();
            torch::Tensor tensor = torch::from_blob(const_cast<float*>(data.data()), {1, 1, static_cast<long>(img.height()), static_cast<long>(img.width())})
                                      .to(torch::kFloat32);
            tensor = (tensor - tensor.mean()) / tensor.std();
            torch::NoGradGuard no_grad;
            torch::Tensor output = resnet_model_.forward({tensor}).toTensor();
            float score = torch::softmax(output, 1)[0][1].item<float>();
            scores.push_back(score);
        }
        return scores;
    }

    void ContrastiveSelection::setBranches(TTree* _tree) {
        _tree->Branch("plane_scores", &plane_scores_);
        _tree->Branch("final_score", &final_score_, "final_score/F");
    }

    void ContrastiveSelection::resetTTree(TTree* _tree) {
        plane_scores_.clear();
        final_score_ = -1.0f;
    }

    DEFINE_ART_CLASS_TOOL(ContrastiveSelection)
}

#endif