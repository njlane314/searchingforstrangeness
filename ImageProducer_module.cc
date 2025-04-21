#ifndef WIREIMAGEPRODUCER_H
#define WIREIMAGEPRODUCER_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "ImageProcessor.h"
#include <vector>
#include <string>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "EventClassifier.h"
#include "LabelClassifier.h"

namespace image {
    class ImageProducer : public art::EDProducer {
    public:
        explicit ImageProducer(fhicl::ParameterSet const& pset);
        void produce(art::Event& e) override;

    private:
        std::string selection_filter_label_;
        std::string pfp_tag_;
        std::string wire_tag_;
        std::string hit_tag_;
        std::string mcparticle_tag_;
        float adc_threshold_;
        std::vector<image::ImageProperties> image_properties_;

        void construct_images(const art::Event& e,
                            const std::vector<ImageProperties>& properties,
                            const std::vector<art::Ptr<recob::Wire>>& wires,
                            const art::FindManyP<recob::Hit>& wire_hit_assoc,
                            const art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>& mcp_bkth_assoc,
                            const signature::Pattern& pattern,
                            const std::vector<signature::Label>& particle_labels,
                            const geo::GeometryCore* geom,
                            const std::vector<art::Ptr<simb::MCParticle>>& mcp_vec,
                            const std::vector<bool>& bad_channel_mask,
                            std::vector<Image>& input_images,
                            std::vector<Image>& truth_images,
                            std::vector<Image>& label_images,
                            float adc_threshold = 10.0f) {
            input_images.clear();
            truth_images.clear();
            label_images.clear();

            for (const auto& prop : properties) {
                input_images.emplace_back(prop); 
                Image truth_img(prop);
                truth_img.clear(signature::kEmptySignature); 
                truth_images.push_back(std::move(truth_img));
                Image label_img(prop);
                label_img.clear(static_cast<float>(signature::Label::undefined)); 
                label_images.push_back(std::move(label_img));
            }

            std::array<std::unordered_map<int, signature::SignatureType>, 3> truth_maps;
            for (const auto& particle : mcp_vec) {
                for (auto& map : truth_maps) {
                    map[particle->TrackId()] = signature::kNeutrinoSignature; 
                }
            }
            for (const auto& [type, sig] : pattern) {
                for (const auto& particle : sig) {
                    for (auto& map : truth_maps) {
                        map[particle->TrackId()] = type; 
                    }
                }
            }

            std::unordered_map<int, signature::Label> trackid_to_label;
            for (size_t i = 0; i < mcp_vec.size() && i < particle_labels.size(); ++i) {
                trackid_to_label[mcp_vec[i]->TrackId()] = particle_labels[i];
            }

            auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
            auto const& wireReadout = art::ServiceHandle<geo::WireReadout const>()->Get();

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

                for (const auto& range : wire->SignalROI().get_ranges()) {
                    const auto& adcs = range.data();
                    int start_tick = range.begin_index();
                    for (size_t idx = 0; idx < adcs.size(); ++idx) {
                        int tick = start_tick + idx;
                        double x = det_props.ConvertTicksToX(tick, wire_ids.front().planeID());
                        size_t row = properties[view_idx].row(x);
                        size_t col = properties[view_idx].col(wire_coord);
                        if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue; 

                        signature::SignatureType truth_label = signature::kEmptySignature;
                        signature::Label particle_label = signature::Label::undefined;

                        for (const auto& hit : hits) {
                            if (tick >= hit->StartTick() && tick < hit->EndTick()) {
                                auto data = mcp_bkth_assoc.data(hit.key());
                                if (data.empty()) {
                                    truth_label = signature::kCosmicSignature;
                                    particle_label = signature::Label::cosmic;
                                }
                                else {
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

                        if (adcs[idx] > adc_threshold) {
                            input_images[view_idx].set(row, col, adcs[idx]);
                            truth_images[view_idx].set(row, col, static_cast<float>(truth_label), false);
                            label_images[view_idx].set(row, col, static_cast<float>(particle_label), false);
                        }
                    }
                }
            }
        }

        std::vector<std::vector<float>> extractImages(const std::vector<image::Image>& images) {
            std::vector<std::vector<float>> data;
            data.reserve(images.size());
            for (const auto& img : images) {
                data.push_back(img.data());
            }
            return data;
        }
    };

    ImageProducer::ImageProducer(fhicl::ParameterSet const& pset)
        : EDProducer{pset},
          selection_filter_label_(pset.get<std::string>("SelectionFilterLabel")),
          pfp_tag_(pset.get<std::string>("PfpTag")),
          wire_tag_(pset.get<std::string>("WireTag")),
          hit_tag_(pset.get<std::string>("HitTag")),
          mcparticle_tag_(pset.get<std::string>("MCParticleTag")),
          adc_threshold_(pset.get<float>("ADCThreshold"))
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

        produces<std::vector<image::Image>>("inputImages");
        produces<std::vector<image::Image>>("truthImages");
        produces<std::vector<image::Image>>("labelImages");
        produces<art::Assns<recob::PFParticle, image::Image>>("inputAssns");
        produces<art::Assns<recob::PFParticle, image::Image>>("truthAssns");
        produces<art::Assns<recob::PFParticle, image::Image>>("labelAssns");
    }

    void ImageProducer::produce(art::Event& e) {
        // Check if event passed the initial selection filter
        auto filter_handle = e.getHandle<bool>(selection_filter_label_);
        if (!filter_handle.isValid() || !*filter_handle) {
            return; // Skip if event didn't pass initial selection
        }

        // Retrieve PFParticles
        auto pfp_handle = e.getValidHandle<std::vector<recob::PFParticle>>(pfp_tag_);
        std::vector<art::Ptr<recob::PFParticle>> pfp_vec;
        art::fill_ptr_vector(pfp_vec, pfp_handle);

        // Retrieve wires, hits, and MC particles
        auto wires = e.getValidHandle<std::vector<recob::Wire>>(wire_tag_);
        auto hits = e.getValidHandle<std::vector<recob::Hit>>(hit_tag_);
        auto mcparticles = e.getValidHandle<std::vector<simb::MCParticle>>(mcparticle_tag_);

        art::FindManyP<recob::Hit> wire_hit_assoc(wires, e, hit_tag_);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hits, e, mcparticle_tag_);

        // Initialize output collections
        auto input_images = std::make_unique<std::vector<image::Image>>();
        auto truth_images = std::make_unique<std::vector<image::Image>>();
        auto label_images = std::make_unique<std::vector<image::Image>>();
        auto input_assns = std::make_unique<art::Assns<recob::PFParticle, image::Image>>();
        auto truth_assns = std::make_unique<art::Assns<recob::PFParticle, image::Image>>();
        auto label_assns = std::make_unique<art::Assns<recob::PFParticle, image::Image>>();

        // Geometry and classifier
        auto geom = art::ServiceHandle<geo::Geometry>()->provider();
        signature::EventClassifier event_classifier(fhicl::ParameterSet());
        signature::LabelClassifier label_classifier(fhicl::ParameterSet());
        signature::Pattern pattern = event_classifier.getPattern(e);
        std::vector<signature::Label> particle_labels = label_classifier.classifyParticles(e);
        std::vector<bool> bad_channel_mask(geom->Nchannels(), false);

        // Identify slices (simplified: assume primary PFPs define slices)
        for (size_t i = 0; i < pfp_vec.size(); ++i) {
            if (!pfp_vec[i]->IsPrimary()) continue; // Process only primary PFPs

            // Collect hits for this slice (simplified)
            std::vector<art::Ptr<recob::Hit>> slice_hits;
            art::FindManyP<recob::Hit> pfp_hit_assoc(pfp_vec, e, hit_tag_);
            auto hits = pfp_hit_assoc.at(i);
            slice_hits.insert(slice_hits.end(), hits.begin(), hits.end());

            if (slice_hits.empty()) continue;

            // Calculate charge-weighted centroid
            double total_charge = 0.0;
            double x_sum = 0.0, y_sum = 0.0, z_sum = 0.0;
            auto det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
            for (const auto& hit : slice_hits) {
                double charge = hit->Integral();
                total_charge += charge;
                double x = det_props.ConvertTicksToX(hit->PeakTime(), hit->WireID().Plane);
                geo::WireID wire_id(hit->WireID());
                auto wire_pos = geom->WireIDToWireGeo(wire_id).GetCenter();
                double y = wire_pos.Y();
                double z = wire_pos.Z();
                x_sum += charge * x;
                y_sum += charge * y;
                z_sum += charge * z;
            }

            if (total_charge <= 0) continue;

            double x_centroid = x_sum / total_charge;
            double y_centroid = y_sum / total_charge;
            double z_centroid = z_sum / total_charge;

            // Build image properties centered on centroid
            std::vector<image::ImageProperties> properties;
            for (const auto& prop : image_properties_) {
                properties.emplace_back(
                    x_centroid, z_centroid, // Adjust y_centroid if needed
                    prop.width(), prop.height(),
                    prop.pixel_h(), prop.pixel_w(),
                    prop.view(), prop.kernel()
                );
            }

            // Generate images
            std::vector<image::Image> slice_input_imgs, slice_truth_imgs, slice_label_imgs;
            image::construct_images(
                e, properties, *wires, wire_hit_assoc, mcp_bkth_assoc,
                pattern, particle_labels, geom, *mcparticles, bad_channel_mask,
                slice_input_imgs, slice_truth_imgs, slice_label_imgs, adc_threshold_
            );

            // Store images and create associations
            size_t start_idx = input_images->size();
            input_images->insert(input_images->end(), slice_input_imgs.begin(), slice_input_imgs.end());
            truth_images->insert(truth_images->end(), slice_truth_imgs.begin(), slice_truth_imgs.end());
            label_images->insert(label_images->end(), slice_label_imgs.begin(), slice_label_imgs.end());

            for (size_t j = 0; j < slice_input_imgs.size(); ++j) {
                input_assns->addSingle(pfp_vec[i], art::Ptr<image::Image>(input_images.get(), start_idx + j));
                truth_assns->addSingle(pfp_vec[i], art::Ptr<image::Image>(truth_images.get(), start_idx + j));
                label_assns->addSingle(pfp_vec[i], art::Ptr<image::Image>(label_images.get(), start_idx + j));
            }
        }

        // Put data products into the event
        e.put(std::move(input_images), "inputImages");
        e.put(std::move(truth_images), "truthImages");
        e.put(std::move(label_images), "labelImages");
        e.put(std::move(input_assns), "inputAssns");
        e.put(std::move(truth_assns), "truthAssns");
        e.put(std::move(label_assns), "labelAssns");
    }
}

DEFINE_ART_MODULE(image::ImageProducer)
