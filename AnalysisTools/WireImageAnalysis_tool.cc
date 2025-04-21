#ifndef WIREIMAGESAVER_H
#define WIREIMAGESAVER_H

#include "AnalysisToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "ImageProcessor.h"
#include <vector>
#include <string>

namespace analysis
{
    class WireImageAnalysis : public AnalysisToolBase {
    public:
        explicit WireImageAnalysis(const fhicl::ParameterSet& pset);
        void configure(fhicl::ParameterSet const& pset) override;
        void analyzeEvent(art::Event const& e, bool is_data) override;
        void analyzeSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override {}
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        std::string input_images_label_;
        std::string truth_images_label_;
        std::string label_images_label_;
        std::vector<std::vector<float>> input_images_;
        std::vector<std::vector<float>> truth_images_;
        std::vector<std::vector<float>> label_images_;
    };

    WireImageAnalysis::WireImageAnalysis(const fhicl::ParameterSet& pset) {
        configure(pset);
    }

    void WireImageAnalysis::configure(fhicl::ParameterSet const& pset) {
        input_images_label_ = pset.get<std::string>("InputImagesLabel", "inputImages");
        truth_images_label_ = pset.get<std::string>("TruthImagesLabel", "truthImages");
        label_images_label_ = pset.get<std::string>("LabelImagesLabel", "labelImages");
    }

    void WireImageAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("input_images", &input_images_);
        _tree->Branch("truth_images", &truth_images_);
        _tree->Branch("label_images", &label_images_);
    }

    void WireImageAnalysis::resetTTree(TTree* _tree) {
        input_images_.clear();
        truth_images_.clear();
        label_images_.clear();
    }

    void WireImageAnalysis::analyzeEvent(art::Event const& e, bool is_data) {
        art::Handle<std::vector<image::Image>> input_handle;
        e.getByLabel(input_images_label_, input_handle);
        if (input_handle.isValid()) {
            for (const auto& img : *input_handle) {
                input_images_.push_back(img.data());
            }
        } else {
            mf::LogWarning("WireImageAnalysis") << "No input images found with label: " << input_images_label_;
        }

        if (!is_data) {
            art::Handle<std::vector<image::Image>> truth_handle;
            e.getByLabel(truth_images_label_, truth_handle);
            if (truth_handle.isValid()) {
                for (const auto& img : *truth_handle) {
                    truth_images_.push_back(img.data());
                }
            } else {
                mf::LogWarning("WireImageAnalysis") << "No truth images found with label: " << truth_images_label_;
            }

            art::Handle<std::vector<image::Image>> label_handle;
            e.getByLabel(label_images_label_, label_handle);
            if (label_handle.isValid()) {
                for (const auto& img : *label_handle) {
                    label_images_.push_back(img.data());
                }
            } else {
                mf::LogWarning("WireImageAnalysis") << "No label images found with label: " << label_images_label_;
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(WireImageAnalysis)
}

#endif