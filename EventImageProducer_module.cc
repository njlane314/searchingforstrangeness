#ifndef EVENTIMAGEPRODUCER_H
#define EVENTIMAGEPRODUCER_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "detinfo/DetectorPropertiesService.h"
#include "TVector3.h"

#include "CommonDefs/Image.h"
#include "CommonDefs/Centroid.h"
#include "CommonDefs/ReconstructionLabeler.h"
#include "CommonDefs/PrimaryLabeler.h"

class EventImageProducer : public art::EDProducer {
public:
    explicit EventImageProducer(fhicl::ParameterSet const& pset);
    void produce(art::Event& e) override;

private:
    void collectSlicePFPs(
        const art::Ptr<recob::PFParticle>& pfp,
        const std::vector<art::Ptr<recob::PFParticle>>& pfp_vec,
        const std::map<unsigned int, size_t>& pfpmap,
        std::vector<art::Ptr<recob::PFParticle>>& slice_pfps
    );

    std::string selection_filter_label_;
    std::string _PFPproducer;
    std::string _WIREproducer;
    std::string hit_tag_;
    std::string _HITproducer;
    float q_threshold_;
    art::InputTag particle_label_;
    double gamma_threshold_;
    double hadron_threshold_;
    std::vector<image::ImageProperties> image_properties_;
    art::InputTag mctruth_label_;
};

EventImageProducer::EventImageProducer(fhicl::ParameterSet const& pset)
    : EDProducer{pset},
      selection_filter_label_(pset.get<std::string>("SelectionFilterLabel")),
      _PFPproducer(pset.get<std::string>("PfpTag")),
      _WIREproducer(pset.get<std::string>("WireTag")),
      hit_tag_(pset.get<std::string>("HitTag")),
      _HITproducer(pset.get<std::string>("MCParticleTag")),
      q_threshold_(pset.get<float>("ADCThreshold")),
      particle_label_(pset.get<art::InputTag>("MCParticleTag", "largeant")),
      gamma_threshold_(pset.get<double>("gamma_threshold", 0.02)),
      hadron_threshold_(pset.get<double>("hadron_threshold", 0.2)),
      mctruth_label_(pset.get<art::InputTag>("mctruth_label", "generator"))
{
    auto props_pset = pset.get<std::vector<fhicl::ParameterSet>>("ImageProperties");
    for (const auto& prop_pset : props_pset) {
        image_properties_.emplace_back(
            prop_pset.get<double>("center_x"),
            prop_pset.get<double>("center_y"),
            prop_pset.get<size_t>("width"),
            prop_pset.get<size_t>("height"),
            prop_pset.get<double>("pixel_h"),
            prop_pset.get<double>("pixel_w"),
            static_cast<geo::View_t>(prop_pset.get<int>("view")),
            prop_pset.get<int>("kernel")
        );
    }

    produces<std::vector<image::Image>>("neutrinoInputImages");
    produces<std::vector<image::Image>>("neutrinoTruthImages");
    produces<std::vector<image::Image>>("neutrinoLabelImages");
    produces<art::Assns<recob::PFParticle, image::Image>>("neutrinoInputAssns");
    produces<art::Assns<recob::PFParticle, image::Image>>("neutrinoTruthAssns");
    produces<art::Assns<recob::PFParticle, image::Image>>("neutrinoLabelAssns");
}

void EventImageProducer::collectSlicePFPs(
    const art::Ptr<recob::PFParticle>& pfp,
    const std::vector<art::Ptr<recob::PFParticle>>& pfp_vec,
    const std::map<unsigned int, size_t>& pfpmap,
    std::vector<art::Ptr<recob::PFParticle>>& slice_pfps
) {
    slice_pfps.push_back(pfp);
    for (unsigned int d : pfp->Daughters()) {
        auto it = pfpmap.find(d);
        if (it != pfpmap.end()) {
            collectSlicePFPs(pfp_vec[it->second], pfp_vec, pfpmap, slice_pfps);
        }
    }
}

void EventImageProducer::produce(art::Event& e) {
    auto filter_handle = e.getHandle<bool>(selection_filter_label_);
    if (!filter_handle.isValid() || !*filter_handle) {
        return;
    }

    auto selectedPFPs = e.getValidHandle<std::vector<art::Ptr<recob::PFParticle>>>(selection_filter_label_ + ".selectedNeutrinoPFPs");
    auto pfp_handle = e.getValidHandle<std::vector<recob::PFParticle>>(_PFPproducer);
    std::vector<art::Ptr<recob::PFParticle>> pfp_vec;
    art::fill_ptr_vector(pfp_vec, pfp_handle);
    std::map<unsigned int, size_t> pfpmap;
    for (size_t i = 0; i < pfp_vec.size(); ++i) {
        pfpmap[pfp_vec[i]->Self()] = i;
    }

    auto wires = e.getValidHandle<std::vector<recob::Wire>>(_WIREproducer);
    auto hits = e.getValidHandle<std::vector<recob::Hit>>(hit_tag_);
    auto mcparticles = e.getValidHandle<std::vector<simb::MCParticle>>(_HITproducer);

    art::FindManyP<recob::Hit> wire_hit_assoc(wires, e, hit_tag_);
    art::FindManyP<recob::Hit> pfp_hit_assoc(pfp_vec, e, hit_tag_);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hits, e, _HITproducer);

    auto geom = art::ServiceHandle<geo::Geometry>()->provider();
    std::vector<reco_labels::ReconstructionLabel> particle_labels = 
        reco_labels::classifyParticles(e, particle_label_, gamma_threshold_, hadron_threshold_);
    std::vector<bool> bad_channel_mask(geom->Nchannels(), false);

    std::vector<truth_labels::PrimaryLabel> primary_labels = 
        truth_labels::classifyPrimaryParticles(e, mctruth_label_, particle_label_);

    std::unordered_map<int, truth_labels::PrimaryLabel> trackid_to_primary_label;
    for (size_t i = 0; i < mcparticles->size() && i < primary_labels.size(); ++i) {
        trackid_to_primary_label[(*mcparticles)[i].TrackId()] = primary_labels[i];
    }

    std::unique_ptr<std::vector<image::Image>> input_images(new std::vector<image::Image>);
    std::unique_ptr<std::vector<image::Image>> truth_images(new std::vector<image::Image>);
    std::unique_ptr<std::vector<image::Image>> label_images(new std::vector<image::Image>);
    std::unique_ptr<art::Assns<recob::PFParticle, image::Image>> input_assns(new art::Assns<recob::PFParticle, image::Image>);
    std::unique_ptr<art::Assns<recob::PFParticle, image::Image>> truth_assns(new art::Assns<recob::PFParticle, image::Image>);
    std::unique_ptr<art::Assns<recob::PFParticle, image::Image>> label_assns(new art::Assns<recob::PFParticle, image::Image>);

    for (const auto& pfp : *selectedPFPs) {
        std::vector<art::Ptr<recob::PFParticle>> slice_pfps;
        collectSlicePFPs(pfp, pfp_vec, pfpmap, slice_pfps);

        std::vector<art::Ptr<recob::Hit>> slice_hits;
        for (const auto& slice_pfp : slice_pfps) {
            size_t index = slice_pfp.key();
            auto pfp_hits = pfp_hit_assoc.at(index);
            slice_hits.insert(slice_hits.end(), pfp_hits.begin(), pfp_hits.end());
        }

        if (slice_hits.empty()) continue;

        TVector3 centroid = common::calculateChargeCentroid(e, slice_hits);
        if (centroid.X() == -999) continue;

        std::vector<image::ImageProperties> properties;
        for (const auto& prop : image_properties_) {
            properties.emplace_back(
                centroid.X(), centroid.Z(),
                prop.width(), prop.height(),
                prop.pixel_h(), prop.pixel_w(),
                prop.view(), prop.kernel()
            );
        }

        std::vector<image::Image> slice_input_imgs, slice_truth_imgs, slice_label_imgs;
        for (const auto& prop : properties) {
            slice_input_imgs.emplace_back(prop);
            image::Image truth_img(prop);
            truth_img.clear(static_cast<float>(truth_labels::PrimaryLabel::undefined));
            slice_truth_imgs.push_back(std::move(truth_img));
            image::Image label_img(prop);
            label_img.clear(static_cast<float>(reco_labels::ReconstructionLabel::undefined));
            slice_label_imgs.push_back(std::move(label_img));
        }

        auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

        for (size_t wire_idx = 0; wire_idx < wires->size(); ++wire_idx) {
            const auto& wire = (*wires)[wire_idx];
            auto ch_id = wire.Channel();
            if (ch_id < bad_channel_mask.size() && bad_channel_mask[ch_id]) continue;

            std::vector<geo::WireID> wire_ids = wireReadout.ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = wireReadout.Plane(wire_ids.front().planeID()).View();
            size_t view_idx = view - geo::kU;

            geo::Point_t center = wireReadout.Wire(wire_ids.front()).GetCenter();
            TVector3 wire_center(center.X(), center.Y(), center.Z());
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                    (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

            const auto& hits = wire_hit_assoc.at(wire_idx);

            for (const auto& range : wire.SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t idx = 0; idx < adcs.size(); ++idx) {
                    int tick = start_tick + idx;
                    double x = det_props.ConvertTicksToX(tick, wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(x);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    truth_labels::PrimaryLabel truth_label = truth_labels::PrimaryLabel::undefined;
                    reco_labels::ReconstructionLabel particle_label = reco_labels::ReconstructionLabel::undefined;

                    for (const auto& hit : hits) {
                        if (tick >= hit->StartTick() && tick < hit->EndTick()) {
                            auto data = mcp_bkth_assoc.data(hit.key());
                            if (data.empty()) {
                                truth_label = truth_labels::PrimaryLabel::other;
                                particle_label = reco_labels::ReconstructionLabel::cosmic;
                            } else {
                                for (size_t i = 0; i < data.size(); ++i) {
                                    if (data[i]->isMaxIDE == 1) {
                                        int track_id = mcp_bkth_assoc.at(hit.key())[i]->TrackId();
                                        if (auto it = trackid_to_primary_label.find(track_id); it != trackid_to_primary_label.end()) {
                                            truth_label = it->second;
                                        }
                                        break;
                                    }
                                }
                            }
                            break;
                        }
                    }

                    if (adcs[idx] > q_threshold_) {
                        slice_input_imgs[view_idx].set(row, col, adcs[idx]);
                        slice_truth_imgs[view_idx].set(row, col, static_cast<float>(truth_label), false);
                        slice_label_imgs[view_idx].set(row, col, static_cast<float>(particle_label), false);
                    }
                }
            }
        }

        size_t start_idx = input_images->size();
        input_images->insert(input_images->end(), slice_input_imgs.begin(), slice_input_imgs.end());
        truth_images->insert(truth_images->end(), slice_truth_imgs.begin(), slice_truth_imgs.end());
        label_images->insert(label_images->end(), slice_label_imgs.begin(), slice_label_imgs.end());

        for (size_t j = 0; j < slice_input_imgs.size(); ++j) {
            input_assns->addSingle(pfp, art::Ptr<image::Image>(input_images, start_idx + j));
            truth_assns->addSingle(pfp, art::Ptr<image::Image>(truth_images, start_idx + j));
            label_assns->addSingle(pfp, art::Ptr<image::Image>(label_images, start_idx + j));
        }
    }

    e.put(std::move(input_images), "neutrinoInputImages");
    e.put(std::move(truth_images), "neutrinoTruthImages");
    e.put(std::move(label_images), "neutrinoLabelImages");
    e.put(std::move(input_assns), "neutrinoInputAssns");
    e.put(std::move(truth_assns), "neutrinoTruthAssns");
    e.put(std::move(label_assns), "neutrinoLabelAssns");
}

DEFINE_ART_MODULE(EventImageProducer)
#endif