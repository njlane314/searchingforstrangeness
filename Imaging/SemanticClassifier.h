#ifndef SEMANTICCLASSIFIER_H
#define SEMANTICCLASSIFIER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/FindManyInChainP.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

#include "Imaging/Image.h"
#include "Common/PandoraUtilities.h"
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace image {
class SemanticClassifier {
    public:
    enum class SemanticLabel {
        Empty = 0,
        Cosmic,
        Muon,
        Electron,
        Photon,
        ChargedPion,
        NeutralPion,
        Neutron,
        Proton,
        ChargedKaon,
        NeutralKaon,
        Lambda,
        ChargedSigma,
        NeutralSigma,
        Other,
        Ambiguous
    };

    static inline const std::array<std::string, 16> semantic_label_names = {
        "Empty",    "Cosmic",        "Muon",       "Electron",
        "Photon",   "ChargedPion",   "NeutralPion", "Neutron",
        "Proton",   "ChargedKaon",   "NeutralKaon", "Lambda",
        "ChargedSigma", "NeutralSigma", "Other", "Ambiguous"};

    explicit SemanticClassifier(const art::InputTag &MCPproducer)
        : fMCPproducer(MCPproducer) {}

    SemanticLabel getSemanticLabel(int pdg) const {
        pdg = std::abs(pdg);
        if (pdg == 13)
            return SemanticLabel::Muon;
        if (pdg == 11)
            return SemanticLabel::Electron;
        if (pdg == 22)
            return SemanticLabel::Photon;
        if (pdg == 2112)
            return SemanticLabel::Neutron;
        if (pdg == 211)
            return SemanticLabel::ChargedPion;
        if (pdg == 111)
            return SemanticLabel::NeutralPion;
        if (pdg == 321)
            return SemanticLabel::ChargedKaon;
        if (pdg == 311 || pdg == 130 || pdg == 310)
            return SemanticLabel::NeutralKaon;
        if (pdg == 2212)
            return SemanticLabel::Proton;
        if (pdg == 3122)
            return SemanticLabel::Lambda;
        if (pdg == 3222 || pdg == 3112)
            return SemanticLabel::ChargedSigma;
        if (pdg == 3212)
            return SemanticLabel::NeutralSigma;

        return SemanticLabel::Other;
    }

    void assignLabelToProgenyRecursively(
        size_t particle_index,
        const std::vector<simb::MCParticle> &particles,
        std::vector<SemanticLabel> &particle_labels,
        const std::unordered_map<int, size_t> &track_id_to_index,
        SemanticLabel primary_label_to_assign) const {
        if (particle_index >= particles.size() || particle_index >= particle_labels.size()) {
            return;
        }
        particle_labels[particle_index] = primary_label_to_assign;
        const auto &particle = particles[particle_index];

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

    std::vector<SemanticLabel> classifyParticles(
        const art::Event &event) const {
        const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
        const auto &particles = *particle_collection_handle;

        std::unordered_map<int, size_t> track_id_to_vector_index;
        for (size_t i = 0; i < particles.size(); ++i) {
            track_id_to_vector_index[particles[i].TrackId()] = i;
        }

        std::vector<SemanticLabel> classified_particle_labels(particles.size(), SemanticLabel::Empty);

        for (size_t i = 0; i < particles.size(); ++i) {
            if (particles[i].Mother() == 0) {
                if (auto it = track_id_to_vector_index.find(particles[i].TrackId()); it != track_id_to_vector_index.end()) {
                    SemanticLabel initial_label = getSemanticLabel(particles[i].PdgCode());
                    assignLabelToProgenyRecursively(it->second, particles, classified_particle_labels, track_id_to_vector_index, initial_label);
                }
            }
        }
        return classified_particle_labels;
    }

    private:
    art::InputTag fMCPproducer;
};
} // namespace image

#endif // SEMANTICCLASSIFIER_H

