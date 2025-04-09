#ifndef HITIMAGEGENERATOR_H
#define HITIMAGEGENERATOR_H

#include "ImageGeneratorBase.h"
#include "ImageProcessor.h"
#include <vector>
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace analysis 
{
    class HitImageAnalysis : public ImageGeneratorBase {
    public:
        HitImageAnalysis(const fhicl::ParameterSet& pset) : ImageGeneratorBase(pset), _pset(pset) {}
        void analyseEvent(const art::Event& e, bool is_data) override {}
        void analyseSlice(const art::Event& e, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        fhicl::ParameterSet _pset;
        std::vector<std::vector<float>> slice_hit_images_;
        std::vector<std::vector<float>> slice_truth_hit_images_;
        bool contained_;
    };

    void HitImageAnalysis::analyseSlice(const art::Event& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
        if (is_data) return;
        auto neutrino_hits = this->collectNeutrinoHits(e, slice_pfp_v);
        std::vector<double> centroid_wires(planes_.size()), centroid_drifts(planes_.size());
        for (size_t i = 0; i < planes_.size(); ++i) {
            auto view = static_cast<common::PandoraView>(planes_[i]);
            auto [wire, drift] = this->calculateCentroid(e, view, neutrino_hits);
            centroid_wires[i] = wire;
            centroid_drifts[i] = drift;
        }
        auto properties = this->buildImageProperties(centroid_wires, centroid_drifts);
        slice_hit_images_ = image::extractImages(image::constructInputHitImages(properties, neutrino_hits, *_geo));
        contained_ = this->isNeutrinoContained(e, properties, centroid_drifts);

        std::vector<art::Ptr<recob::Hit>> hit_vec;
        if (auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(_HITproducer))
            art::fill_ptr_vector(hit_vec, hitHandle);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hit_vec, e, _BKTproducer);
        signature::EventClassifier classifier(_pset);
        signature::Pattern pattern = classifier.getPattern(e);
        std::vector<art::Ptr<simb::MCParticle>> mcp_vec;
        if (auto mcpHandle = e.getValidHandle<std::vector<simb::MCParticle>>(_MCTproducer))
            art::fill_ptr_vector(mcp_vec, mcpHandle);

        slice_truth_hit_images_ = image::extractImages(
            image::constructTruthHitImages(properties, hit_vec, mcp_bkth_assoc, pattern, classifier, *_geo, mcp_vec, _bad_channel_mask)
        );
    }

    void HitImageAnalysis::setBranches(TTree* tree) {
        tree->Branch("slice_hit_images", &slice_hit_images_);
        tree->Branch("slice_truth_hit_images", &slice_truth_hit_images_);
        tree->Branch("slice_hit_contained", &contained_);
    }

    void HitImageAnalysis::resetTTree(TTree* tree) {
        slice_hit_images_.clear();
        slice_truth_hit_images_.clear();
        contained_ = false;
    }

    DEFINE_ART_CLASS_TOOL(HitImageAnalysis)
}

#endif