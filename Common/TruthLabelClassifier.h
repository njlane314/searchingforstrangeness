#ifndef TRUTHLABELCLASSIFIER_H
#define TRUTHLABELCLASSIFIER_H

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

#include "ImageTypes.h"
#include "PandoraUtilities.h"
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

namespace analysis {
class TruthLabelClassifier {
    public:
    enum class TruthPrimaryLabel {
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
        Other
    };

    static inline const std::array<std::string, 15> truth_primary_label_names = {
        "Empty", "Cosmic", "Muon", "Electron", "Photon",
        "ChargedPion", "NeutralPion", "Neutron", "Proton",
        "ChargedKaon", "NeutralKaon", "Lambda", "ChargedSigma",
        "NeutralSigma", "Other"};

    explicit TruthLabelClassifier(const art::InputTag &MCPproducer)
        : fMCPproducer(MCPproducer) {}

    TruthPrimaryLabel getTruthPrimaryLabel(int pdg) const {
        pdg = std::abs(pdg);
        if (pdg == 13)
            return TruthPrimaryLabel::Muon;
        if (pdg == 11)
            return TruthPrimaryLabel::Electron;
        if (pdg == 22)
            return TruthPrimaryLabel::Photon;
        if (pdg == 2112)
            return TruthPrimaryLabel::Neutron;
        if (pdg == 211)
            return TruthPrimaryLabel::ChargedPion;
        if (pdg == 111)
            return TruthPrimaryLabel::NeutralPion;
        if (pdg == 321)
            return TruthPrimaryLabel::ChargedKaon;
        if (pdg == 311 || pdg == 130 || pdg == 310)
            return TruthPrimaryLabel::NeutralKaon;
        if (pdg == 2212)
            return TruthPrimaryLabel::Proton;
        if (pdg == 3122)
            return TruthPrimaryLabel::Lambda;
        if (pdg == 3222 || pdg == 3112)
            return TruthPrimaryLabel::ChargedSigma;
        if (pdg == 3212)
            return TruthPrimaryLabel::NeutralSigma;

        return TruthPrimaryLabel::Other;
    }

    void assignLabelToProgenyRecursively(
        size_t particle_index,
        const std::vector<simb::MCParticle> &particles,
        std::vector<TruthPrimaryLabel> &particle_labels,
        const std::unordered_map<int, size_t> &track_id_to_index,
        TruthPrimaryLabel primary_label_to_assign) const {
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

    std::vector<TruthPrimaryLabel> classifyParticles(
        const art::Event &event) const {
        const auto particle_collection_handle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
        const auto &particles = *particle_collection_handle;

        std::unordered_map<int, size_t> track_id_to_vector_index;
        for (size_t i = 0; i < particles.size(); ++i) {
            track_id_to_vector_index[particles[i].TrackId()] = i;
        }

        std::vector<TruthPrimaryLabel> classified_particle_labels(particles.size(), TruthPrimaryLabel::Empty);

        for (size_t i = 0; i < particles.size(); ++i) {
            if (particles[i].Mother() == 0) {
                if (auto it = track_id_to_vector_index.find(particles[i].TrackId()); it != track_id_to_vector_index.end()) {
                    TruthPrimaryLabel initial_label = getTruthPrimaryLabel(particles[i].PdgCode());
                    assignLabelToProgenyRecursively(it->second, particles, classified_particle_labels, track_id_to_vector_index, initial_label);
                }
            }
        }
        return classified_particle_labels;
    }

    private:
    art::InputTag fMCPproducer;
};
} // namespace analysis

#endif // TRUTHLABELCLASSIFIER_H

