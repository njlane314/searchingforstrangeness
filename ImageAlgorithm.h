#ifndef IMAGEALGO_H
#define IMAGEALGO_H

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

namespace analysis {

enum class ReconstructionLabel {
    Empty,
    Cosmic,
    MIP,
    HIP,
    Shower,
    Michel,
    Diffuse,
    Invisible
};

enum class TruthPrimaryLabel {
    Empty,
    Cosmic,
    Electron,
    Muon,
    ChargedPion,
    NeutralPion,
    Proton,
    ChargedKaon,
    NeutralKaon,
    Lambda,
    ChargedSigma,
    Other
};

class ImageAlgorithm {
public:
    ImageAlgorithm(const fhicl::ParameterSet& p) {
        configure(p);
    }

    void configure(const fhicl::ParameterSet& p) {
        fHitProducer = p.get<art::InputTag>("HITproducer", "gaushit");
        fWireProducer = p.get<art::InputTag>("WIREproducer", "butcher");
        fMCPproducer = p.get<art::InputTag>("MCPproducer", "largeant");
        fBackTrackerLabel = p.get<std::string>("BackTrackerLabel", "gaushit");
        fProcessMC = p.get<bool>("ProcessMC", true);
        _image_width = p.get<int>("ImageWidth", 512);
        _image_height = p.get<int>("ImageHeight", 512);
        _adc_image_threshold = p.get<float>("ADCthreshold", 1.0);
        fGammaThreshold = p.get<double>("GammaThreshold", 0.1);
        fHadronThreshold = p.get<double>("HadronThreshold", 0.1);
        fUseCheatRecoLabels = p.get<bool>("UseCheatRecoLabels", true);

        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
        double tick_period = clock->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1;
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);
    }

    void generateSliceImages(const art::Event& e, const std::vector<art::Ptr<recob::Hit>>& hits,
                             std::vector<Image<float>>& raw_images,
                             std::vector<Image<int>>& reco_images,
                             std::vector<Image<int>>& true_images) {
        auto centroid_u = this->calculateChargeCentroid(e, geo::kU, hits);
        auto centroid_v = this->calculateChargeCentroid(e, geo::kV, hits);
        auto centroid_w = this->calculateChargeCentroid(e, geo::kW, hits);

        std::vector<ImageProperties> properties;
        properties.emplace_back(centroid_u.first, centroid_u.second, _image_width, _image_height, _drift_step, _wire_pitch_u, geo::kU);
        properties.emplace_back(centroid_v.first, centroid_v.second, _image_width, _image_height, _drift_step, _wire_pitch_v, geo::kV);
        properties.emplace_back(centroid_w.first, centroid_w.second, _image_width, _image_height, _drift_step, _wire_pitch_w, geo::kW);

        this->constructPixelImages(e, properties, raw_images, reco_images, true_images);
    }

private:
    art::InputTag fHitProducer;
    art::InputTag fWireProducer;
    art::InputTag fMCPproducer;
    std::string fBackTrackerLabel;
    bool fProcessMC;
    int _image_width;
    int _image_height;
    float _adc_image_threshold;
    double fGammaThreshold;
    double fHadronThreshold;
    bool fUseCheatRecoLabels;

    const geo::GeometryCore* _geo;
    const detinfo::DetectorProperties* _detp;
    float _drift_step;
    float _wire_pitch_u;
    float _wire_pitch_v;
    float _wire_pitch_w;

    std::vector<TruthPrimaryLabel> classified_true_labels;
    std::vector<ReconstructionLabel> classified_reco_labels;
    std::map<int, size_t> trackid_to_index;

    std::pair<double, double> calculateChargeCentroid(const art::Event& e, geo::View_t view, const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0;
        double sum_wire_coord = 0.0;
        double sum_drift_coord = 0.0;
        double wire_pitch = (view == geo::kU) ? _wire_pitch_u : (view == geo::kV) ? _wire_pitch_v : _wire_pitch_w;
        for (const auto& hit : hits) {
            if (_geo->View(hit->Channel()) != view) continue;
            double charge = hit->Integral();
            double wire_number = hit->WireID().Wire;
            double wire_coord = wire_number * wire_pitch;
            double drift_coord = _detp->ConvertTicksToX(hit->PeakTime(), hit->WireID().planeID());
            sum_charge += charge;
            sum_wire_coord += wire_coord * charge;
            sum_drift_coord += drift_coord * charge;
        }
        if (sum_charge == 0.0) return {0.0, 0.0};
        return {sum_wire_coord / sum_charge, sum_drift_coord / sum_charge};
    }

    void constructPixelImages(const art::Event& e, const std::vector<ImageProperties>& properties,
                          std::vector<Image<float>>& raw_images,
                          std::vector<Image<int>>& reco_images,
                          std::vector<Image<int>>& true_images) {
        for (const auto& prop : properties) {
            Image<float> raw_image(prop);
            Image<int> reco_image(prop);
            Image<int> true_image(prop);
            raw_image.clear(0.0);
            reco_image.clear(static_cast<int>(ReconstructionLabel::Empty));
            true_image.clear(static_cast<int>(TruthPrimaryLabel::Empty));
            raw_images.push_back(std::move(raw_image));
            reco_images.push_back(std::move(reco_image));
            true_images.push_back(std::move(true_image));
        }

        bool hasMCInfo = false;
        art::Handle<std::vector<simb::MCParticle>> mcpHandle;
        if (e.getByLabel(fMCPproducer, mcpHandle)) {
            if (mcpHandle.isValid() && !mcpHandle->empty()) {
                hasMCInfo = true;
            }
        }

        if (fProcessMC && hasMCInfo) {
            classified_true_labels = this->classifyTruthParticles(e);
            classified_reco_labels = this->classifyRecoParticles(e);
            trackid_to_index.clear();
            const auto& particles = *mcpHandle;
            for (size_t i = 0; i < particles.size(); ++i) {
                trackid_to_index[particles[i].TrackId()] = i;
            }
        } else {
            classified_true_labels.clear();
            classified_reco_labels.clear();
            trackid_to_index.clear();
        }

        auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(fWireProducer);
        auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitProducer);
        art::FindManyP<recob::Hit> wire_hit_assoc(wireHandle, e, fHitProducer);

        for (size_t wire_idx = 0; wire_idx < wireHandle->size(); ++wire_idx) {
            const auto& wire = wireHandle->at(wire_idx);
            auto ch_id = wire.Channel();
            std::vector<geo::WireID> wire_ids = _geo->ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = _geo->View(wire_ids.front().planeID());
            size_t view_idx = static_cast<size_t>(view);
            double wire_pitch = (view == geo::kU) ? _wire_pitch_u : (view == geo::kV) ? _wire_pitch_v : _wire_pitch_w;
            double wire_coord = wire_ids.front().Wire * wire_pitch;

            auto hits_for_wire = wire_hit_assoc.at(wire_idx);
            std::vector<art::Ptr<recob::Hit>> sorted_hits = hits_for_wire;
            std::sort(sorted_hits.begin(), sorted_hits.end(), [](const art::Ptr<recob::Hit>& a, const art::Ptr<recob::Hit>& b) {
                return a->StartTick() < b->StartTick();
            });

            size_t current_hit_idx = 0;

            for (const auto& range : wire.SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t idx = 0; idx < adcs.size(); ++idx) {
                    int tick = start_tick + idx;
                    double drift_coord = _detp->ConvertTicksToX(static_cast<double>(tick), wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(drift_coord);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                    if (adcs[idx] > _adc_image_threshold) {
                        float adc = std::log10(adcs[idx]);
                        raw_images[view_idx].set(row, col, adc);
                        int reco_label = static_cast<int>(ReconstructionLabel::Cosmic);
                        int true_label = static_cast<int>(TruthPrimaryLabel::Cosmic);

                        if (hasMCInfo) {
                            auto it = std::upper_bound(
                                sorted_hits.begin() + current_hit_idx, sorted_hits.end(), tick,
                                [](int t, const art::Ptr<recob::Hit>& hit) { return t < hit->StartTick(); }
                            );

                            current_hit_idx = std::distance(sorted_hits.begin(), it);
                            if (current_hit_idx > 0) {
                                size_t prev_idx = current_hit_idx - 1;
                                const auto& prev_hit = sorted_hits[prev_idx];
                                if (prev_hit->StartTick() <= tick && tick < prev_hit->EndTick()) {
                                    int temp_reco = this->getRecoLabel(e, prev_hit);
                                    int temp_true = this->getTrueLabel(e, prev_hit);
                                    if (temp_reco != static_cast<int>(ReconstructionLabel::Empty)) {
                                        reco_label = temp_reco;
                                    }
                                    if (temp_true != static_cast<int>(TruthPrimaryLabel::Empty)) {
                                        true_label = temp_true;
                                    }
                                    return; 
                                }
                            }

                            if (current_hit_idx < sorted_hits.size()) {
                                const auto& hit = sorted_hits[current_hit_idx];
                                if (hit->StartTick() <= tick && tick < hit->EndTick()) {
                                    int temp_reco = this->getRecoLabel(e, hit);
                                    int temp_true = this->getTrueLabel(e, hit);
                                    if (temp_reco != static_cast<int>(ReconstructionLabel::Empty)) {
                                        reco_label = temp_reco;
                                    }
                                    if (temp_true != static_cast<int>(TruthPrimaryLabel::Empty)) {
                                        true_label = temp_true;
                                    }
                                }
                            }
                        }
                        reco_images[view_idx].set(row, col, reco_label);
                        true_images[view_idx].set(row, col, true_label);
                    } else {
                        reco_images[view_idx].set(row, col, static_cast<int>(ReconstructionLabel::Empty));
                        true_images[view_idx].set(row, col, static_cast<int>(TruthPrimaryLabel::Empty));
                    }
                }
            }
        }
    }

    TruthPrimaryLabel getTruthPrimaryLabelFromPDG(int pdg) const {
        pdg = std::abs(pdg);
        if (pdg == 11) return TruthPrimaryLabel::Electron;
        if (pdg == 13) return TruthPrimaryLabel::Muon;
        if (pdg == 211) return TruthPrimaryLabel::ChargedPion;
        if (pdg == 111) return TruthPrimaryLabel::NeutralPion;
        if (pdg == 2212) return TruthPrimaryLabel::Proton;
        if (pdg == 321) return TruthPrimaryLabel::ChargedKaon;
        if (pdg == 311 || pdg == 130 || pdg == 310) return TruthPrimaryLabel::NeutralKaon;
        if (pdg == 3122) return TruthPrimaryLabel::Lambda;
        if (pdg == 3222 || pdg == 3112 || pdg == 3212) return TruthPrimaryLabel::ChargedSigma;
        return TruthPrimaryLabel::Other;
    }

    void assignTruthLabelRecursively(
        size_t particle_index,
        const std::vector<simb::MCParticle>& particles,
        std::vector<TruthPrimaryLabel>& particle_labels,
        const std::unordered_map<int, size_t>& track_id_to_index,
        TruthPrimaryLabel label_to_assign
    ) const {
        std::cout << "truth loop" << std::endl;
        if (particle_index >= particles.size() || particle_index >= particle_labels.size()) return;
        particle_labels[particle_index] = label_to_assign;
        const auto& particle = particles[particle_index];
        for (int daughter_idx = 0; daughter_idx < particle.NumberDaughters(); ++daughter_idx) {
            int daughter_track_id = particle.Daughter(daughter_idx);
            auto it = track_id_to_index.find(daughter_track_id);
            if (it != track_id_to_index.end() && it->second < particles.size()) {
                this->assignTruthLabelRecursively(it->second, particles, particle_labels, track_id_to_index, label_to_assign);
            }
        }
    }

    std::vector<TruthPrimaryLabel> classifyTruthParticles(const art::Event& event) const {
        auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
        const auto& particles = *particle_collection_handle;
        std::unordered_map<int, size_t> track_id_to_index;
        for (size_t i = 0; i < particles.size(); ++i) {
            track_id_to_index[particles[i].TrackId()] = i;
        }
        std::vector<TruthPrimaryLabel> labels(particles.size(), TruthPrimaryLabel::Empty);
        for (size_t i = 0; i < particles.size(); ++i) {
            if (particles[i].Mother() == 0) {
                auto it = track_id_to_index.find(particles[i].TrackId());
                if (it != track_id_to_index.end()) {
                    TruthPrimaryLabel label = getTruthPrimaryLabelFromPDG(particles[i].PdgCode());
                    this->assignTruthLabelRecursively(it->second, particles, labels, track_id_to_index, label);
                }
            }
        }
        return labels;
    }

    std::pair<ReconstructionLabel, ReconstructionLabel> getReconstructionLabelAndPropagation(
        const std::vector<simb::MCParticle>& particles,
        const simb::MCParticle& particle,
        ReconstructionLabel parent_label,
        const std::map<int, size_t>& track_id_to_index
    ) const {
        if (parent_label != ReconstructionLabel::Empty) return {parent_label, parent_label};
        ReconstructionLabel label = ReconstructionLabel::Invisible;
        ReconstructionLabel child_label = ReconstructionLabel::Empty;
        int pdg = particle.PdgCode();
        double momentum = particle.P();
        std::string start_process = particle.Process();
        std::string end_process = particle.EndProcess();
        int parent_track_id = particle.Mother();
        int parent_pdg = 0;
        if (parent_track_id != 0) {
            auto it = track_id_to_index.find(parent_track_id);
            if (it != track_id_to_index.end() && it->second < particles.size()) {
                parent_pdg = particles[it->second].PdgCode();
            }
        }
        if (std::abs(pdg) == 211 || std::abs(pdg) == 13) {
            label = ReconstructionLabel::MIP;
        } else if (std::abs(pdg) == 321 || (std::abs(pdg) == 2212 && momentum >= fHadronThreshold)) {
            label = ReconstructionLabel::HIP;
        } else if (std::abs(pdg) == 11) {
            if (start_process == "primary") {
                label = ReconstructionLabel::Shower;
                child_label = ReconstructionLabel::Shower;
            } else if (std::abs(parent_pdg) == 13 && (start_process == "muMinusCaptureAtRest" || start_process == "muPlusCaptureAtRest" || start_process == "Decay")) {
                label = ReconstructionLabel::Michel;
                child_label = ReconstructionLabel::Michel;
            } else if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                if (momentum >= fGammaThreshold) {
                    label = ReconstructionLabel::Shower;
                    child_label = ReconstructionLabel::Shower;
                } else {
                    label = ReconstructionLabel::Diffuse;
                }
            } else {
                label = ReconstructionLabel::Diffuse;
            }
        } else if (pdg == 22) {
            if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt" || start_process == "primary") {
                if (momentum >= fGammaThreshold) {
                    label = ReconstructionLabel::Shower;
                    child_label = ReconstructionLabel::Shower;
                } else {
                    label = ReconstructionLabel::Diffuse;
                }
            } else {
                label = ReconstructionLabel::Diffuse;
            }
        } else if (std::abs(pdg) == 2212 && momentum < fHadronThreshold) {
            label = ReconstructionLabel::Diffuse;
        }
        return {label, child_label};
    }

    void assignRecoLabelRecursively(
        size_t particle_index,
        const std::vector<simb::MCParticle>& particles,
        ReconstructionLabel parent_label,
        std::vector<ReconstructionLabel>& labels,
        const std::map<int, size_t>& track_id_to_index
    ) const {
        std::cout << "reco loop" << std::endl;
        if (particle_index >= particles.size() || particle_index >= labels.size()) return;
        const auto& particle = particles[particle_index];
        auto [label, child_label] = this->getReconstructionLabelAndPropagation(particles, particle, parent_label, track_id_to_index);
        labels[particle_index] = label;
        for (int i = 0; i < particle.NumberDaughters(); ++i) {
            int daughter_track_id = particle.Daughter(i);
            auto it = track_id_to_index.find(daughter_track_id);
            if (it != track_id_to_index.end() && it->second < particles.size()) {
                this->assignRecoLabelRecursively(it->second, particles, child_label, labels, track_id_to_index);
            }
        }
    }

    std::vector<ReconstructionLabel> classifyRecoParticles(const art::Event& event) const {
        if (fUseCheatRecoLabels) {
            auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
            const auto& particles = *particle_collection_handle;
            std::map<int, size_t> track_id_to_index;
            for (size_t i = 0; i < particles.size(); ++i) {
                track_id_to_index[particles[i].TrackId()] = i;
            }
            std::vector<ReconstructionLabel> labels(particles.size(), ReconstructionLabel::Empty);
            for (size_t i = 0; i < particles.size(); ++i) {
                if (particles[i].Mother() == 0) {
                    this->assignRecoLabelRecursively(i, particles, ReconstructionLabel::Empty, labels, track_id_to_index);
                }
            }
            return labels;
        } else {
            return {};
        }
    }

    int getRecoLabel(const art::Event& e, const art::Ptr<recob::Hit>& hit) {
        if (!fProcessMC) return static_cast<int>(ReconstructionLabel::Empty);

        auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitProducer);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> fm(hitHandle, e, fBackTrackerLabel);
        
        const auto& mcParts = fm.at(hit.key());
        const auto& datas = fm.data(hit.key());

        if (!mcParts.empty()) {
            for (size_t i = 0; i < mcParts.size(); ++i) {
                if (datas[i] && datas[i]->isMaxIDE) {
                    int track_id = mcParts[i]->TrackId();
                    auto it_trackid = trackid_to_index.find(track_id);
                    if (it_trackid != trackid_to_index.end() && it_trackid->second < classified_reco_labels.size()) {
                        return static_cast<int>(classified_true_labels[it_trackid->second]);
                    }
                    break; 
                }
            }
        }
        return static_cast<int>(ReconstructionLabel::Empty);
    }

    int getTrueLabel(const art::Event& e, const art::Ptr<recob::Hit>& hit) {
        if (!fProcessMC) return static_cast<int>(TruthPrimaryLabel::Empty);

        auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitProducer);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> fm(hitHandle, e, fBackTrackerLabel);
        
        const auto& mcParts = fm.at(hit.key());
        const auto& datas = fm.data(hit.key());

        if (!mcParts.empty()) {
            for (size_t i = 0; i < mcParts.size(); ++i) {
                if (datas[i] && datas[i]->isMaxIDE) {
                    int track_id = mcParts[i]->TrackId();
                    auto it_trackid = trackid_to_index.find(track_id);
                    if (it_trackid != trackid_to_index.end() && it_trackid->second < classified_true_labels.size()) {
                        return static_cast<int>(classified_true_labels[it_trackid->second]);
                    }
                    break; 
                }
            }
        }
        return static_cast<int>(TruthPrimaryLabel::Empty);
    }
};

}

#endif // IMAGEALGO_H