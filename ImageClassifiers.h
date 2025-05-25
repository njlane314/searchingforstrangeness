#ifndef IMAGE_CLASSIFIERS_H
#define IMAGE_CLASSIFIERS_H

#include <vector>
#include <string>
#include <map>
#include <array>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <limits>
#include <optional>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include "lardataobj/AnalysisBase/MVAOutput.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "CommonDefs/Pandora.h"
#include "CommonDefs/Types.h"

namespace analysis 
{
    class TruthLabelClassifier {
    public:
        enum class TruthPrimaryLabel {
            Empty = 0,
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

        static inline const std::array<std::string, 12> truth_primary_label_names = {
            "Empty", "Cosmic", "Electron", "Muon", "ChargedPion", "NeutralPion", "Proton", "ChargedKaon",
            "NeutralKaon", "Lambda", "ChargedSigma", "Other"
        };

        explicit TruthLabelClassifier(const art::InputTag& MCPproducer)
        : fMCPproducer(MCPproducer) {}

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

        void assignLabelToProgenyRecursively(
            size_t particle_index,
            const std::vector<simb::MCParticle>& particles,
            std::vector<TruthPrimaryLabel>& particle_labels,
            const std::unordered_map<int, size_t>& track_id_to_index,
            TruthPrimaryLabel primary_label_to_assign
        ) const {
            if (particle_index >= particles.size() || particle_index >= particle_labels.size()) {
                return;
            }
            particle_labels[particle_index] = primary_label_to_assign;
            const auto& particle = particles[particle_index];

            for (int daughter_idx = 0; daughter_idx < particle.NumberDaughters(); ++daughter_idx) {
                int daughter_track_id = particle.Daughter(daughter_idx);
                auto it = track_id_to_index.find(daughter_track_id);
                if (it != track_id_to_index.end()) {
                    if (it->second < particles.size()) {
                        assignLabelToProgenyRecursively(it->second, particles, particle_labels, track_id_to_index, primary_label_to_assign);
                    }
                }
            }
        }

        std::vector<TruthPrimaryLabel> classifyParticles(
            const art::Event& event
        ) const {
            const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
            const auto& particles = *particle_collection_handle;

            std::unordered_map<int, size_t> track_id_to_vector_index;
            for (size_t i = 0; i < particles.size(); ++i) {
                track_id_to_vector_index[particles[i].TrackId()] = i;
            }

            std::vector<TruthPrimaryLabel> classified_particle_labels(particles.size(), TruthPrimaryLabel::Empty);

            for (size_t i = 0; i < particles.size(); ++i) {
                if (particles[i].Mother() == 0) {
                    if (auto it = track_id_to_vector_index.find(particles[i].TrackId()); it != track_id_to_vector_index.end()) {
                        TruthPrimaryLabel initial_label = getTruthPrimaryLabelFromPDG(particles[i].PdgCode());
                        assignLabelToProgenyRecursively(it->second, particles, classified_particle_labels, track_id_to_vector_index, initial_label);
                    }
                }
            }
            return classified_particle_labels;
        }

    private:
        art::InputTag fMCPproducer;
    };

    class RecoLabelClassifier {
    public:
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

        static inline const std::array<std::string, 8> reco_label_names = {
            "Empty", "Cosmic", "MIP", "HIP", "Shower", "Michel", "Diffuse", "Invisible"
        };

        explicit RecoLabelClassifier(
            const art::InputTag& mc_producer_tag,
            double gamma_threshold,
            double hadron_threshold,
            bool use_cheat)
        : fMCPproducer(mc_producer_tag),
        fGammaThreshold(gamma_threshold),
        fHadronThreshold(hadron_threshold),
        fUseCheat(use_cheat) {}

        std::pair<ReconstructionLabel, ReconstructionLabel> getReconstructionLabelAndPropagation(
            const std::vector<simb::MCParticle>& particles,
            const simb::MCParticle& particle_to_label,
            ReconstructionLabel label_from_parent,
            const std::map<int, size_t>& track_id_to_index
        ) const {
            if (label_from_parent != ReconstructionLabel::Empty) {
                return {label_from_parent, label_from_parent};
            }

            ReconstructionLabel determined_label = ReconstructionLabel::Invisible;
            ReconstructionLabel propagated_label_for_children = ReconstructionLabel::Empty;

            int pdg_code = particle_to_label.PdgCode();
            double momentum = particle_to_label.P();
            const std::string& start_process = particle_to_label.Process();
            const std::string& end_process = particle_to_label.EndProcess();
            int parent_track_id = particle_to_label.Mother();
            int parent_pdg_code = 0;

            if (parent_track_id != 0) {
                auto it = track_id_to_index.find(parent_track_id);
                if (it != track_id_to_index.end()) {
                    if (it->second < particles.size()) {
                        parent_pdg_code = particles[it->second].PdgCode();
                    }
                }
            }

            if (std::abs(pdg_code) == 211 || std::abs(pdg_code) == 13) {
                determined_label = ReconstructionLabel::MIP;
            } else if (std::abs(pdg_code) == 321 || (std::abs(pdg_code) == 2212 && momentum >= fHadronThreshold)) {
                determined_label = ReconstructionLabel::HIP;
            } else if (std::abs(pdg_code) == 11) {
                if (start_process == "primary") {
                    determined_label = ReconstructionLabel::Shower;
                    propagated_label_for_children = ReconstructionLabel::Shower;
                } else if (std::abs(parent_pdg_code) == 13 &&
                        (start_process == "muMinusCaptureAtRest" ||
                            start_process == "muPlusCaptureAtRest" ||
                            start_process == "Decay")) {
                    determined_label = ReconstructionLabel::Michel;
                    propagated_label_for_children = ReconstructionLabel::Michel;
                } else if (start_process == "conv" || end_process == "conv" ||
                        start_process == "compt" || end_process == "compt") {
                    if (momentum >= fGammaThreshold) {
                        determined_label = ReconstructionLabel::Shower;
                        propagated_label_for_children = ReconstructionLabel::Shower;
                    } else {
                        determined_label = ReconstructionLabel::Diffuse;
                    }
                } else {
                    determined_label = ReconstructionLabel::Diffuse;
                }
            } else if (pdg_code == 22) {
                if (start_process == "conv" || end_process == "conv" ||
                    start_process == "compt" || end_process == "compt" ||
                    start_process == "primary") {
                    if (momentum >= fGammaThreshold) {
                        determined_label = ReconstructionLabel::Shower;
                        propagated_label_for_children = ReconstructionLabel::Shower;
                    } else {
                        determined_label = ReconstructionLabel::Diffuse;
                    }
                } else {
                    determined_label = ReconstructionLabel::Diffuse;
                }
            } else if (std::abs(pdg_code) == 2212 && momentum < fHadronThreshold) {
                determined_label = ReconstructionLabel::Diffuse;
            }
            return {determined_label, propagated_label_for_children};
        }

        void assignLabelToProgenyRecursively(
            size_t particle_index,
            const std::vector<simb::MCParticle>& particles,
            ReconstructionLabel label_from_parent,
            std::vector<ReconstructionLabel>& particle_labels,
            const std::map<int, size_t>& track_id_to_index
        ) const {
            if (particle_index >= particles.size() || particle_index >= particle_labels.size()) {
                return;
            }
            const auto& current_particle = particles[particle_index];
            auto [label_for_current, label_for_children] = getReconstructionLabelAndPropagation(
                particles, current_particle, label_from_parent, track_id_to_index
            );
            particle_labels[particle_index] = label_for_current;

            for (int i = 0; i < current_particle.NumberDaughters(); ++i) {
                int daughter_track_id = current_particle.Daughter(i);
                auto it = track_id_to_index.find(daughter_track_id);
                if (it != track_id_to_index.end()) {
                    if (it->second < particles.size()) {
                        assignLabelToProgenyRecursively(it->second, particles, label_for_children, particle_labels, track_id_to_index);
                    }
                }
            }
        }

        std::vector<ReconstructionLabel> classifyParticlesFromReco(const art::Event& event) const {
            return {};
        }

        std::vector<ReconstructionLabel> classifyParticles(const art::Event& event) const {
            if (fUseCheat) {
                const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
                const auto& particles = *particle_collection_handle;

                std::map<int, size_t> track_id_to_vector_index;
                for (size_t i = 0; i < particles.size(); ++i) {
                    track_id_to_vector_index[particles[i].TrackId()] = i;
                }

                std::vector<ReconstructionLabel> classified_particle_labels(particles.size(), ReconstructionLabel::Empty);

                for (size_t i = 0; i < particles.size(); ++i) {
                    if (particles[i].Mother() == 0) {
                        assignLabelToProgenyRecursively(i, particles, ReconstructionLabel::Empty, classified_particle_labels, track_id_to_vector_index);
                    }
                }
                return classified_particle_labels;
            } else {
                return classifyParticlesFromReco(event);
            }
        }

        art::InputTag fMCPproducer;
        double fGammaThreshold;
        double fHadronThreshold;
        bool fUseCheat;
    };
}

#endif