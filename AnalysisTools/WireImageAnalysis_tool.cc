#ifndef WIREIMAGEGENERATOR_H
#define WIREIMAGEGENERATOR_H

#include "ImageGeneratorBase.h"
#include "../ImageProcessor.h"
#include "../LabelClassifier.h"
#include <vector>
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace analysis 
{
    class WireImageAnalysis : public ImageGeneratorBase {
    public:
        WireImageAnalysis(const fhicl::ParameterSet& pset) 
            : ImageGeneratorBase(pset), _pset(pset), _event_classifier(pset.get<fhicl::ParameterSet>("EventClassifier")), _label_classifier(pset.get<fhicl::ParameterSet>("LabelClassifier")) {}
        void analyseEvent(const art::Event& e, bool is_data) override;
        void analyseSlice(const art::Event& e, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        fhicl::ParameterSet _pset;
        signature::EventClassifier _event_classifier;
        signature::LabelClassifier _label_classifier;
        std::vector<std::vector<float>> event_wire_images_;
        std::vector<std::vector<float>> slice_wire_images_;
        std::vector<std::vector<float>> event_truth_wire_images_;
        std::vector<std::vector<float>> slice_truth_wire_images_;
        std::vector<std::vector<float>> event_label_wire_images_; 
        std::vector<std::vector<float>> slice_label_wire_images_; 
        bool contained_;
    };

    void WireImageAnalysis::analyseEvent(const art::Event& e, bool is_data) {}

    void WireImageAnalysis::analyseSlice(const art::Event& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
        if (is_data) return;

        std::vector<art::Ptr<recob::Wire>> wire_vec;
        if (auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(_WREproducer)) 
            art::fill_ptr_vector(wire_vec, wireHandle);

        art::FindManyP<recob::Hit> wire_hit_assoc(wire_vec, e, _HITproducer);

        auto neutrino_hits = this->collectNeutrinoHits(e, slice_pfp_v);
        std::vector<double> centroid_wires(planes_.size()), centroid_drifts(planes_.size());
        for (size_t i = 0; i < planes_.size(); ++i) {
            auto view = static_cast<common::PandoraView>(planes_[i]);
            auto [wire, drift] = this->calculateCentroid(e, view, neutrino_hits);
            centroid_wires[i] = wire;
            centroid_drifts[i] = drift;
        }
        auto properties = this->buildImageProperties(centroid_wires, centroid_drifts);

        std::vector<art::Ptr<recob::Hit>> hit_vec;
        if (auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(_HITproducer))
            art::fill_ptr_vector(hit_vec, hitHandle);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hit_vec, e, _BKTproducer);
        signature::Pattern pattern = _event_classifier.getPattern(e);
        std::vector<art::Ptr<simb::MCParticle>> mcp_vec;
        if (auto mcpHandle = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer))
            art::fill_ptr_vector(mcp_vec, mcpHandle);
        std::vector<signature::Label> particle_labels = _label_classifier.classifyParticles(e);

        float adc_threshold = _pset.get<float>("ADCThreshold", 0.0f);

        std::vector<image::Image> input_imgs, truth_imgs, label_imgs;
        image::constructAllWireImages(
            e, properties, wire_vec, wire_hit_assoc, mcp_bkth_assoc, pattern, particle_labels,
            _geo, mcp_vec, _bad_channel_mask, input_imgs, truth_imgs, label_imgs, adc_threshold
        );

        slice_wire_images_ = image::extractImages(input_imgs);
        slice_truth_wire_images_ = image::extractImages(truth_imgs);
        slice_label_wire_images_ = image::extractImages(label_imgs);
        contained_ = this->isNeutrinoContained(e, properties, centroid_drifts);
    }

    void WireImageAnalysis::setBranches(TTree* tree) {
        tree->Branch("event_wire_images", &event_wire_images_);
        tree->Branch("slice_wire_images", &slice_wire_images_);
        tree->Branch("event_truth_wire_images", &event_truth_wire_images_);
        tree->Branch("slice_truth_wire_images", &slice_truth_wire_images_);
        tree->Branch("event_label_wire_images", &event_label_wire_images_); 
        tree->Branch("slice_label_wire_images", &slice_label_wire_images_); 
        tree->Branch("slice_wire_contained", &contained_);
    }

    void WireImageAnalysis::resetTTree(TTree* tree) {
        event_wire_images_.clear();
        slice_wire_images_.clear();
        event_truth_wire_images_.clear();
        slice_truth_wire_images_.clear();
        event_label_wire_images_.clear(); 
        slice_label_wire_images_.clear();
        contained_ = false;
    }

    DEFINE_ART_CLASS_TOOL(WireImageAnalysis)
}

#endif