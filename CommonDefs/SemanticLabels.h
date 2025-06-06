#ifndef PRIMARYLABELER_H
#define PRIMARYLABELER_H

#include <vector>
#include <map>
#include <string>
#include <array>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace common {

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

inline SemanticLabel GetSemanticLabel(int pdg) {
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

inline void ProcessSemanticParticle(
    size_t particle_index,
    const std::vector<simb::MCParticle>& particles,
    std::vector<SemanticLabel>& particle_labels,
    const std::unordered_map<int, size_t>& track_id_to_index,
    SemanticLabel primary_label
) {
    particle_labels[particle_index] = primary_label; 
    const auto& particle = particles[particle_index];

    for (int daughter_index = 0; daughter_index < particle.NumberDaughters(); ++daughter_index) {
        if (auto it = track_id_to_index.find(particle.Daughter(daughter_index)); it != track_id_to_index.end()) {
            ProcessSemanticParticle(it->second, particles, particle_labels, track_id_to_index, primary_label);
        }
    }
}

inline std::vector<SemanticLabel> ClassifySemanticParticles(
    const art::Event& event,
    const art::InputTag& particle_label
) {
    auto const& particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(particle_label);
    const auto& particles = *particle_handle;

    std::unordered_map<int, size_t> track_id_to_index;
    for (size_t i = 0; i < particles.size(); ++i) {
        track_id_to_index[particles[i].TrackId()] = i;
    }

    std::vector<SemanticLabel> particle_labels(particles.size(), SemanticLabel::kNeutrino);

    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i].Mother() == 0) {  
            if (auto it = track_id_to_index.find(particles[i].TrackId()); it != track_id_to_index.end()) {
                SemanticLabel label = GetSemanticLabel(particles[i].PdgCode());
                ProcessSemanticParticle(it->second, particles, particle_labels, track_id_to_index, label);
            }
        }
    }
    
    return particle_labels;
}


}

#endif