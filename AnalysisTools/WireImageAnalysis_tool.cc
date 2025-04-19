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
            : ImageGeneratorBase(pset), _pset(pset), _event_classifier(pset.get<fhicl::ParameterSet>("EventClassifier")), _label_classifier(pset.get<fhicl::ParameterSet>("LabelClassifier")) {
            detector_min_wire_.resize(planes_.size());
            detector_max_wire_.resize(planes_.size());
            for (size_t i = 0; i < planes_.size(); ++i) {
                double min_z = std::numeric_limits<double>::max();
                double max_z = std::numeric_limits<double>::min();
                geo::View_t view = static_cast<geo::View_t>(planes_[i]);
                for (geo::TPCID::TPCID_t tpc = 0; tpc < _geo->NTPC(); ++tpc) {
                    for (geo::WireID::WireID_t wire = 0; wire < _geo->Nwires(geo::PlaneID(0, tpc, planes_[i])); ++wire) {
                        geo::WireID wire_id(0, tpc, planes_[i], wire);
                        TVector3 wire_center = _geo->WireIDToWireGeo(wire_id).GetCenter();
                        double transformed_Z = (view == geo::kW) ? wire_center.Z() :
                                               (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                                   (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));
                        min_z = std::min(min_z, transformed_Z);
                        max_z = std::max(max_z, transformed_Z);
                    }
                }
                detector_min_wire_[i] = min_z;
                detector_max_wire_[i] = max_z;
            }
            detector_max_drift_ = _detp->ConvertTicksToX(_detp->NumberTimeSamples(), geo::PlaneID(0, 0, 0));
        }
        void analyseEvent(const art::Event& e, bool is_data) override;
        void analyseSlice(const art::Event& e, std::vector<common::ProxyPfpElem_t>& slicePfpVector, bool is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        fhicl::ParameterSet _pset;
        signature::EventClassifier _event_classifier;
        signature::LabelClassifier _label_classifier;
        std::vector<double> detector_min_wire_;
        std::vector<double> detector_max_wire_;
        double detector_max_drift_;
        std::vector<std::vector<float>> event_wire_images_;
        std::vector<std::vector<float>> slice_wire_images_;
        std::vector<std::vector<float>> event_truth_wire_images_;
        std::vector<std::vector<float>> slice_truth_wire_images_;
        std::vector<std::vector<float>> event_label_wire_images_; 
        std::vector<std::vector<float>> slice_label_wire_images_; 
        bool contained_;
    };

    void WireImageAnalysis::analyseEvent(const art::Event& e, bool is_data) {
        if (is_data) return;
        this->loadBadChannelMask(e);

        std::vector<art::Ptr<recob::Wire>> wire_vec;
        if (auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(_WREproducer)) 
            art::fill_ptr_vector(wire_vec, wireHandle);

        std::vector<image::ImageProperties> properties;
        for (size_t i = 0; i < planes_.size(); ++i) {
            double pixel_w = (planes_[i] == 0) ? _wire_pitch_u : (planes_[i] == 1) ? _wire_pitch_v : _wire_pitch_w;
            double pixel_h = _drift_step;
            size_t w = static_cast<size_t>((detector_max_wire_[i] - detector_min_wire_[i]) / pixel_w);
            size_t h = static_cast<size_t>(detector_max_drift_ / pixel_h);
            properties.emplace_back(detector_min_wire_[i], 0.0, h, w, pixel_w, pixel_h, static_cast<geo::View_t>(planes_[i]), 0);
        }

        event_wire_images_ = image::extractImages(
            image::constructInputWireImages(properties, wire_vec, *_geo, _bad_channel_mask)
        );
        
        std::vector<art::Ptr<recob::Hit>> hit_vec;
        if (auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(_HITproducer))
            art::fill_ptr_vector(hit_vec, hitHandle);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hit_vec, e, _BKTproducer);

        art::FindManyP<recob::Hit> wire_hit_assoc(wire_vec, e, _HITproducer); 
 
        signature::Pattern pattern = _event_classifier.getPattern(e);
        std::vector<art::Ptr<simb::MCParticle>> mcp_vec;
        if (auto mcpHandle = e.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer))
            art::fill_ptr_vector(mcp_vec, mcpHandle);

        event_truth_wire_images_ = image::extractImages(
            image::constructTruthWireImages(properties, wire_vec, wire_hit_assoc, mcp_bkth_assoc, pattern, _event_classifier, *_geo, mcp_vec, _bad_channel_mask)
        );

        std::vector<signature::Label> particle_labels = _label_classifier.classifyParticles(e);
        event_label_wire_images_ = image::extractImages(
            image::constructLabelWireImages(properties, wire_vec, wire_hit_assoc, mcp_bkth_assoc, particle_labels, mcp_vec, *_geo, _bad_channel_mask)
        );
    }

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
            properties, wire_vec, wire_hit_assoc, mcp_bkth_assoc, pattern, particle_labels,
            *_geo, mcp_vec, _bad_channel_mask, input_imgs, truth_imgs, label_imgs, adc_threshold
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