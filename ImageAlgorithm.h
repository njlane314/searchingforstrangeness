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
#include "CommonDefs/Pandora.h"

namespace analysis 
{
    enum class SemanticLabel {
        kEmpty = 0,
        kCosmic,
        kNeutrino,
        kMuon,
        kProton,
        kPion,
        kChargedKaon, 
        kNeutralKaon,
        kLambda,
        kChargedSigma
    };

    class ImageAlgorithm {
    public:
        ImageAlgorithm(const fhicl::ParameterSet& pset) {
            this->configure(pset);
        }

        void configure(const fhicl::ParameterSet& pset) {
            _image_height = pset.get<size_t>("ImageHeight", 512);
            _image_width = pset.get<size_t>("ImageWidth", 512);
            fWIREproducer = pset.get<art::InputTag>("WireProducer", "butcher");
            fHITproducer = pset.get<art::InputTag>("HitProducer", "gaushit");
            fBKTproducer = pset.get<art::InputTag>("BackTrackerProducer", "gaushitTruthMatch");
        }

        void generateViewImages(
            const art::Event& event,            
            const std::vector<ImageProperties>& image_props,     
            const std::vector<art::Ptr<recob::Wire>>& wires,           
            const std::vector<art::Ptr<simb::MCParticle>>& mcps) {
            auto geo = art::ServiceHandle<geo::Geometry>()->provider();
            auto detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
            auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
            double time_tick = clock->TPCClock().TickPeriod();
            double drift_velocity = detp->DriftVelocity();
            _drift_step = time_tick * drift_velocity * 1e1;
            _wire_pitch_u = geo->WirePitch(geo::kU);
            _wire_pitch_v = geo->WirePitch(geo::kV);
            _wire_pitch_w = geo->WirePitch(geo::kW);
            
            std::unordered_map<int, SemanticLabel> track_id_to_label = classifyParticles(mcps);
            
            art::FindManyP<recob::Hit> wire_to_hits(wires, event, fHITproducer);
            std::vector<std::vector<art::Ptr<recob::Hit>>> sorted_hits(wires.size());
            std::vector<art::Ptr<recob::Hit>> all_hits;
            for (size_t i = 0; i < wires.size(); ++i) {
                sorted_hits[i] = wire_to_hits.at(i);
                std::sort(sorted_hits[i].begin(), sorted_hits[i].end(),
                    [](const auto& a, const auto& b) { return a->StartTick() < b->StartTick(); });
                all_hits.insert(all_hits.end(), sorted_hits[i].begin(), sorted_hits[i].end());
            }
            art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> hit_to_mcparticles(all_hits, event, fBKTproducer);
            std::unordered_map<art::Ptr<recob::Hit>, SemanticLabel> hit_to_label;
            for (size_t i = 0; i < all_hits.size(); ++i) {
                const auto& mcparticles = hit_to_mcparticles.at(i);
                const auto& matching_data = hit_to_mcparticles.data(i);
                if (!mcparticles.empty()) {
                    auto max_it = std::max_element(matching_data.begin(), matching_data.end(),
                        [](const auto* a, const auto* b) { return a->ideFraction < b->ideFraction; });
                    size_t idx = std::distance(matching_data.begin(), max_it);
                    int track_id = mcparticles[idx]->TrackId();
                    auto it = track_id_to_label.find(track_id);
                    hit_to_label[all_hits[i]] = (it != track_id_to_label.end()) ? it->second : SemanticLabel::kEmpty;
                } else {
                    hit_to_label[all_hits[i]] = SemanticLabel::kEmpty;
                }
            }

            detector_images_.clear();
            semantic_images_.clear();
            for (const auto& prop : image_props) {
                geo::View_t view = prop.view();
                Image<float> detector_image(prop);
                Image<int> semantic_image(prop);

                for (size_t w = 0; w < wires.size(); ++w) {
                    const auto& wire = wires[w];
                    raw::ChannelID_t channel = wire->Channel();
                    std::vector<geo::WireID> wire_ids = geo->ChannelToWire(channel);
                    if (wire_ids.empty() || geo->View(wire_ids.front()) != view) continue;

                    TVector3 wire_center = geo->WireIDToWireGeo(wire_ids.front()).GetCenter();
                    double wire_coord = common::ProjectYZToWire(wire_center.Y(), wire_center.Z(), view);

                    for (const auto& range : wire->SignalROI().get_ranges()) {
                        const auto& adcs = range.data();
                        int start_tick = range.begin_index();
                        for (size_t idx = 0; idx < adcs.size(); ++idx) {
                            if (adcs[idx] <= 1) continue;
                            int tick = start_tick + idx;
                            double drift_coord = detp->ConvertTicksToX(tick, wire_ids.front());
                            size_t row = prop.row(drift_coord);
                            size_t col = prop.col(wire_coord);
                            if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;
                            
                            // std::log10(adcs[idx]) 
                            detector_image.set(row, col, adcs[idx], true);
                            SemanticLabel label = getSemanticLabelForWire(wire, tick, sorted_hits[w], hit_to_label);
                            semantic_image.set(row, col, static_cast<int>(label), false);
                        }
                    }
                }

                detector_images_.push_back(std::move(detector_image));
                semantic_images_.push_back(std::move(semantic_image));
            }
        }

        const std::vector<Image<float>>& getDetectorImages() const {
            return detector_images_;
        }

        const std::vector<Image<int>>& getSemanticImages() const {
            return semantic_images_;
        }

        double getDriftStep() const {
            return _drift_step;
        }

        double getWirePitch(geo::View_t view) const {
            switch (view) {
                case geo::kU: return _wire_pitch_u;
                case geo::kV: return _wire_pitch_v;
                case geo::kW: case geo::kY: return _wire_pitch_w;
                default: return 0.0;
            }
        }

        size_t getImageWidth() const { 
            return _image_width; 
        }
        
        size_t getImageHeight() const { 
            return _image_height; 
        }

    private:
        size_t _image_height;
        size_t _image_width;
        double _drift_step;
        double _wire_pitch_u;
        double _wire_pitch_v;
        double _wire_pitch_w;
        std::vector<Image<float>> detector_images_;
        std::vector<Image<int>> semantic_images_;
        art::InputTag fWIREproducer;
        art::InputTag fHITproducer;
        art::InputTag fBKTproducer;

        static SemanticLabel getSemanticLabel(int pdg) {
            pdg = std::abs(pdg);
            if (pdg == 13) return SemanticLabel::kMuon;
            if (pdg == 2212) return SemanticLabel::kProton;
            if (pdg == 211) return SemanticLabel::kPion;
            if (pdg == 321) return SemanticLabel::kChargedKaon;
            if (pdg == 311) return SemanticLabel::kNeutralKaon;
            if (pdg == 3122) return SemanticLabel::kLambda;
            if (pdg == 3222 || pdg == 3112) return SemanticLabel::kChargedSigma;
            return SemanticLabel::kNeutrino;
        }

        std::unordered_map<int, SemanticLabel> classifyParticles(
            const std::vector<art::Ptr<simb::MCParticle>>& mcps) {
            std::unordered_map<int, SemanticLabel> track_id_to_label;
            for (const auto& mcp : mcps) {
                if (mcp->Mother() == 0) {
                    SemanticLabel label = getSemanticLabel(mcp->PdgCode());
                    track_id_to_label[mcp->TrackId()] = label;
                    for (int d = 0; d < mcp->NumberDaughters(); ++d) {
                        int daughter_id = mcp->Daughter(d);
                        track_id_to_label[daughter_id] = label;
                    }
                }
            }
            return track_id_to_label;
        }

        SemanticLabel getSemanticLabelForWire(
            const art::Ptr<recob::Wire>& wire,
            int tick,
            const std::vector<art::Ptr<recob::Hit>>& hits,
            const std::unordered_map<art::Ptr<recob::Hit>, SemanticLabel>& hit_to_label) {
            auto it = std::lower_bound(hits.begin(), hits.end(), tick,
                [](const art::Ptr<recob::Hit>& h, int t) { return h->EndTick() <= t; });
            if (it != hits.end() && (*it)->StartTick() <= tick && tick < (*it)->EndTick()) {
                auto label_it = hit_to_label.find(*it);
                return (label_it != hit_to_label.end()) ? label_it->second : SemanticLabel::kEmpty;
            }
            return SemanticLabel::kEmpty;
        }
    };
}

#endif