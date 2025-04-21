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
#include "CommonDefs/PFPUtils.h"
#include "ReconstructionLabeler.h"

class EventImageProducer;

class EventImageProducer : public art::EDProducer {
public:
    explicit EventImageProducer(fhicl::ParameterSet const& pset);
    void produce(art::Event& e) override;

private:
    std::string selection_filter_label_; ///< Label for event selection filter
    std::string _PFPproducer;            ///< Producer label for PFParticles
    std::string _WIREproducer;           ///< Producer label for wires
    std::string hit_tag_;                ///< Tag for hits
    std::string _HITproducer;            ///< Producer label for MC particles
    float q_threshold_;                  ///< ADC threshold for hit inclusion
    art::InputTag particle_label_;       ///< Input tag for particle labels
    double gamma_threshold_;             ///< Threshold for gamma classification
    double hadron_threshold_;            ///< Threshold for hadron classification
    std::vector<image::ImageProperties> image_properties_; ///< Predefined image properties

    void construct_images(
        const art::Event& e,
        const std::vector<ImageProperties>& properties,
        const std::vector<art::Ptr<recob::Wire>>& wires,
        const art::FindManyP<recob::Hit>& wire_hit_assoc,
        const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
        const std::vector<reco_labels::ReconstructionLabel>& particle_labels,
        const geo::GeometryCore* geom,
        const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
        const std::vector<bool>& bad_channel_mask,
        float q_threshold,
        std::vector<Image>& input_images,
        std::vector<Image>& truth_images,
        std::vector<Image>& label_images
    );
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
      hadron_threshold_(pset.get<double>("hadron_threshold", 0.2))
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

    // Initialize output collections
    auto input_images = std::make_unique<std::vector<image::Image>>();
    auto truth_images = std::make_unique<std::vector<image::Image>>();
    auto label_images = std::make_unique<std::vector<image::Image>>();
    auto input_assns = std::make_unique<art::Assns<recob::PFParticle, image::Image>>();
    auto truth_assns = std::make_unique<art::Assns<recob::PFParticle, image::Image>>();
    auto label_assns = std::make_unique<art::Assns<recob::PFParticle, image::Image>>();

    // Access geometry and classify particles
    auto geom = art::ServiceHandle<geo::Geometry>()->provider();
    std::vector<reco_labels::ReconstructionLabel> particle_labels = 
        reco_labels::classifyParticles(e, particle_label_, gamma_threshold_, hadron_threshold_);
    std::vector<bool> bad_channel_mask(geom->Nchannels(), false);

    // Process each selected neutrino slice
    for (const auto& pfp : *selectedPFPs) {
        // Collect all PFParticles in the slice
        std::vector<art::Ptr<recob::PFParticle>> slice_pfps;
        common::collectSlicePFPs(pfp, pfp_vec, pfpmap, slice_pfps);

        // Collect hits for the slice
        std::vector<art::Ptr<recob::Hit>> slice_hits;
        for (const auto& slice_pfp : slice_pfps) {
            size_t index = slice_pfp.key();
            auto pfp_hits = pfp_hit_assoc.at(index);
            slice_hits.insert(slice_hits.end(), pfp_hits.begin(), pfp_hits.end());
        }

        if (slice_hits.empty()) continue;

        // Calculate charge centroid
        TVector3 centroid = common::calculateChargeCentroid(e, slice_hits);
        if (centroid.X() == -999) continue; // Skip invalid centroids

        // Adjust image properties to center on centroid
        std::vector<image::ImageProperties> properties;
        for (const auto& prop : image_properties_) {
            properties.emplace_back(
                centroid.X(), centroid.Z(),
                prop.width(), prop.height(),
                prop.pixel_h(), prop.pixel_w(),
                prop.view(), prop.kernel()
            );
        }

        // Generate images for the slice
        std::vector<image::Image> slice_input_imgs, slice_truth_imgs, slice_label_imgs;
        construct_images(
            e, properties, *wires, wire_hit_assoc, mcp_bkth_assoc,
            particle_labels, geom, *mcparticles, bad_channel_mask, q_threshold_,
            slice_input_imgs, slice_truth_imgs, slice_label_imgs
        );

        // Store images and create associations
        size_t start_idx = input_images->size();
        input_images->insert(input_images->end(), slice_input_imgs.begin(), slice_input_imgs.end());
        truth_images->insert(truth_images->end(), slice_truth_imgs.begin(), slice_truth_imgs.end());
        label_images->insert(label_images->end(), slice_label_imgs.begin(), slice_label_imgs.end());

        for (size_t j = 0; j < slice_input_imgs.size(); ++j) {
            input_assns->addSingle(pfp, art::Ptr<image::Image>(input_images.get(), start_idx + j));
            truth_assns->addSingle(pfp, art::Ptr<image::Image>(truth_images.get(), start_idx + j));
            label_assns->addSingle(pfp, art::Ptr<image::Image>(label_images.get(), start_idx + j));
        }
    }

    // Put data products into the event
    e.put(std::move(input_images), "neutrinoInputImages");
    e.put(std::move(truth_images), "neutrinoTruthImages");
    e.put(std::move(label_images), "neutrinoLabelImages");
    e.put(std::move(input_assns), "neutrinoInputAssns");
    e.put(std::move(truth_assns), "neutrinoTruthAssns");
    e.put(std::move(label_assns), "neutrinoLabelAssns");
}

/**
 * @brief Constructs the input, truth, and label images.
 */
void EventImageProducer::construct_images(
    const art::Event& e,
    const std::vector<ImageProperties>& properties,
    const std::vector<art::Ptr<recob::Wire>>& wires,
    const art::FindManyP<recob::Hit>& wire_hit_assoc,
    const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
    const std::vector<reco_labels::ReconstructionLabel>& particle_labels,
    const geo::GeometryCore* geom,
    const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
    const std::vector<bool>& bad_channel_mask,
    float q_threshold,
    std::vector<Image>& input_images,
    std::vector<Image>& truth_images,
    std::vector<Image>& label_images
) {
    input_images.clear();
    truth_images.clear();
    label_images.clear();

    // Initialize images for each view
    for (const auto& prop : properties) {
        input_images.emplace_back(prop);
        Image truth_img(prop);
        truth_img.clear(reco_labels::kEmptySignature);
        truth_images.push_back(std::move(truth_img));
        Image label_img(prop);
        label_img.clear(static_cast<float>(reco_labels::ReconstructionLabel::undefined));
        label_images.push_back(std::move(label_img));
    }

    // Set all particles as neutrino signatures for truth maps
    std::array<std::unordered_map<int, reco_labels::SignatureType>, 3> truth_maps;
    for (const auto& particle : mcp_vec) {
        for (auto& map : truth_maps) {
            map[particle->TrackId()] = reco_labels::kNeutrinoSignature;
        }
    }

    // Map track IDs to reconstruction labels
    std::unordered_map<int, reco_labels::ReconstructionLabel> trackid_to_label;
    for (size_t i = 0; i < mcp_vec.size() && i < particle_labels.size(); ++i) {
        trackid_to_label[mcp_vec[i]->TrackId()] = particle_labels[i];
    }

    auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

    // Process each wire
    for (size_t wire_idx = 0; wire_idx < wires.size(); ++wire_idx) {
        const auto& wire = wires[wire_idx];
        auto ch_id = wire->Channel();
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

        // Process signal ROIs
        for (const auto& range : wire->SignalROI().get_ranges()) {
            const auto& adcs = range.data();
            int start_tick = range.begin_index();
            for (size_t idx = 0; idx < adcs.size(); ++idx) {
                int tick = start_tick + idx;
                double x = det_props.ConvertTicksToX(tick, wire_ids.front().planeID());
                size_t row = properties[view_idx].row(x);
                size_t col = properties[view_idx].col(wire_coord);
                if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                reco_labels::SignatureType truth_label = reco_labels::kEmptySignature;
                reco_labels::ReconstructionLabel particle_label = reco_labels::ReconstructionLabel::undefined;

                // Associate hits with truth and labels
                for (const auto& hit : hits) {
                    if (tick >= hit->StartTick() && tick < hit->EndTick()) {
                        auto data = mcp_bkth_assoc.data(hit.key());
                        if (data.empty()) {
                            truth_label = reco_labels::kCosmicSignature;
                            particle_label = reco_labels::ReconstructionLabel::cosmic;
                        } else {
                            for (size_t i = 0; i < data.size(); ++i) {
                                if (data[i]->isMaxIDE == 1) {
                                    int track_id = mcp_bkth_assoc.at(hit.key())[i]->TrackId();
                                    if (auto it = truth_maps[view_idx].find(track_id); it != truth_maps[view_idx].end()) {
                                        truth_label = it->second;
                                    }
                                    if (auto it = trackid_to_label.find(track_id); it != trackid_to_label.end()) {
                                        particle_label = it->second;
                                    }
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }

                // Set pixel values if above threshold
                if (adcs[idx] > q_threshold) {
                    input_images[view_idx].set(row, col, adcs[idx]);
                    truth_images[view_idx].set(row, col, static_cast<float>(truth_label), false);
                    label_images[view_idx].set(row, col, static_cast<float>(particle_label), false);
                }
            }
        }
    }
}

} // namespace image

DEFINE_ART_MODULE(image::EventImageProducer)

#endif
