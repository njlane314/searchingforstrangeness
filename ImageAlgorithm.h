#ifndef IMAGE_ALGORITHM_H
#define IMAGE_ALGORITHM_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <array>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <limits>
#include <optional>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "Image.h"
#include "CommonDefs/HitGeometry.h"
#include "CommonDefs/SliceProperties.h"
#include "CommonDefs/SemanticLabels.h"

namespace analysis 
{
    class ImageAlgorithm {
    public:
        ImageAlgorithm(const fhicl::ParameterSet& pset) {
            this->configure(pset);
        }

        void generateDetectorPlaneImages(const art::Event& event,
                    const std::vector<art::Ptr<recob::Hit>>& slice_hits,
                    std::vector<Image<float>>& detector_plane_tensor,
                    std::vector<Image<int>>& semantic_plane_tensor) {
            auto geo = art::ServiceHandle<geo::Geometry>()->provider();
            auto detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
     
            std::pair<double, double> charge_centroid_u, charge_centroid_v, charge_centroid_w;
            common::CalculateHitChargeCentroid(event, common::PandoraView::TPC_VIEW_U, slice_hits, charge_centroid_u);
            common::CalculateHitChargeCentroid(event, common::PandoraView::TPC_VIEW_V, slice_hits, charge_centroid_v);
            common::CalculateHitChargeCentroid(event, common::PandoraView::TPC_VIEW_W, slice_hits, charge_centroid_w);

            std::vector<ImageProperties> plane_properties;
            plane_properties.reserve(geo->Nplanes());
            std::pair<double, double> charge_centroid;
            for (size_t i = 0; i < geo->Nplanes(); ++i) {
                geo::View_t view = geo->Plane(i).View();
                double wire_pitch = _wire_pitch_w;
                if (view == geo::kU) {
                    wire_pitch = _wire_pitch_u;
                    charge_centroid = charge_centroid_u;
                } else if (view == geo::kV) {
                    wire_pitch = _wire_pitch_v;
                    charge_centroid = charge_centroid_v;
                } else if (view == geo::kW) {
                    wire_pitch = _wire_pitch_w;
                    charge_centroid = charge_centroid_w;
                } else {
                    throw std::runtime_error("Unknown view type: " + std::to_string(view));
                }
                ImageProperties image_props(
                    charge_centroid.first,
                    charge_centroid.second,
                    _image_width,
                    _image_height,
                    wire_pitch,
                    _drift_step,
                    view
                ); 
                plane_properties.push_back(image_props);
            }

            detector_plane_tensor.clear();
            detector_plane_tensor.reserve(geo->Nplanes());
            semantic_plane_tensor.clear();
            semantic_plane_tensor.reserve(geo->Nplanes());
            for (size_t i = 0; i < geo->Nplanes(); ++i) {
                Image<float> detector_image(plane_properties[i]);
                Image<int> semantic_image(plane_properties[i]);
                detector_image.clear(0.0f);
                semantic_image.clear(0);
                detector_plane_tensor.push_back(detector_image);
                semantic_plane_tensor.push_back(semantic_image);
            }

            auto hits = event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
            if (!hits.isValid()) {
                throw std::runtime_error("Invalid hits handle: " + fHITproducer.encode());
            }
            auto wires = event.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
            if (!wires.isValid()) {
                throw std::runtime_error("Invalid wires handle: " + fWIREproducer.encode());
            }
            auto mcps = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
            if (!mcps.isValid()) {
                throw std::runtime_error("Invalid MCParticles handle: " + fMCPproducer.encode());
            }
            art::FindManyP<recob::Hit> wire_hit_assoc(wires, event, fHITproducer);
            if (!wire_hit_assoc.isValid()) {
                throw std::runtime_error("Invalid wire association handle for hits: " + fHITproducer.encode());
            }
            art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcp_bkth_assoc(hits, event, fBTproducer);
            if (!mcp_bkth_assoc.isValid()) {
                throw std::runtime_error("Invalid BackTrackerMatchingData association handle for hits: " + fBTproducer);
            }

            std::map<int, size_t> track_to_mcp_map;
            for (size_t mcp_index = 0; mcp_index < mcps->size(); ++mcp_index) {
                track_to_mcp_map[mcps->at(mcp_index).TrackId()] = mcp_index;
            }

            std::vector<common::SemanticLabel> semantic_labels = common::ClassifySemanticParticles(event, fMCPproducer);
            for (size_t wire_index = 0; wire_index < wires->size(); ++wire_index) {
                const auto& wire = wires->at(wire_index);
                raw::ChannelID_t channel = wire.Channel();
                std::vector<geo::WireID> wire_ids = geo->ChannelToWire(channel);
                if (wire_ids.empty()) continue;
                geo::View_t view = geo->View(wire_ids.front().planeID());
                size_t view_idx = view - geo::kU;

                const geo::WireGeo* wire_geo = geo->WirePtr(wire_ids.front());
                TVector3 center = wire_geo->GetCenter();
                TVector3 wire_center(center.X(), center.Y(), center.Z());
                double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                    (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                        (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

                const auto& hits = wire_hit_assoc.at(wire_index);

                for (const auto& range : wire.SignalROI().get_ranges()) {
                    const auto& adcs = range.data();
                    int start_tick = range.begin_index();
                    for (size_t idx = 0; idx < adcs.size(); ++idx) {
                        int tick = start_tick + idx;
                        double drift_coord = detp->ConvertTicksToX(tick, wire_ids.front().planeID());
                        size_t row = plane_properties[view_idx].row(drift_coord);
                        size_t col = plane_properties[view_idx].col(wire_coord);
                        if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                        common::SemanticLabel semantic_label = common::SemanticLabel::kCosmic;

                        for (const auto& hit : hits) {
                            if (tick >= hit->StartTick() && tick < hit->EndTick()) {
                                auto bkth_data = mcp_bkth_assoc.data(hit.key()); 
                                if (bkth_data.empty()) continue;
                                for (size_t bkth_index = 0; bkth_index < bkth_data.size(); ++bkth_index) {
                                    int track_id = mcp_bkth_assoc.at(hit.key())[bkth_index]->TrackId();
                                    auto it = track_to_mcp_map.find(track_id);
                                    if (it != track_to_mcp_map.end()) {
                                        size_t mcp_index = it->second;
                                        if (mcp_index < mcps->size()) {
                                            semantic_label = semantic_labels[mcp_index];
                                        }
                                    }
                                }
                            }
                        }

                        if (adcs[idx] > 0) {
                            detector_plane_tensor[view_idx].set(row, col, adcs[idx]);
                            semantic_plane_tensor[view_idx].set(row, col, static_cast<int>(semantic_label));
                        }
                    }
                }
            }
        }

    private: 
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        std::string fBTproducer;

        size_t _image_width;
        size_t _image_height;

        double _drift_step;
        double _wire_pitch_u;
        double _wire_pitch_v;
        double _wire_pitch_w;

        void configure(const fhicl::ParameterSet& pset) {
            fHITproducer = pset.get<art::InputTag>("HITproducer", "gaushit");
            fWIREproducer = pset.get<art::InputTag>("WIREproducer", "buthcer");
            fMCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
            fBTproducer = pset.get<std::string>("BTproducer", "gaushit");

            _image_width = pset.get<size_t>("image_width", 512);
            _image_height = pset.get<size_t>("image_height", 512);
            
            auto geo = art::ServiceHandle<geo::Geometry>()->provider();
            auto detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
            auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
            double time_tick = clock->TPCClock().TickPeriod();

            std::cout << "Time tick: " << time_tick << " seconds" << std::endl;
            double drift_velocity = detp->DriftVelocity();
            std::cout << "Drift velocity: " << drift_velocity << " cm/us" << std::endl;
            _drift_step = time_tick * drift_velocity * 1e1; 
            std::cout << "Drift step: " << _drift_step << " cm" << std::endl;
            _wire_pitch_u = geo->WirePitch(geo::kU);
            _wire_pitch_v = geo->WirePitch(geo::kV);
            _wire_pitch_w = geo->WirePitch(geo::kW);
            std::cout << "Wire pitch U: " << _wire_pitch_u << " cm" << std::endl;
            std::cout << "Wire pitch V: " << _wire_pitch_v << " cm" << std::endl;
            std::cout << "Wire pitch Y: " << _wire_pitch_w << " cm" << std::endl;
        }

    };
}

#endif // IMAGE_ALGORITHM_H

     