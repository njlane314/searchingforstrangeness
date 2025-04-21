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

namespace truth_labels {

enum class PrimaryLabel {
    empty = 0,
    cosmic,
    neutrino,
    proton,
    neutron,
    charged_pion,
    neutral_pion,
    muon,
    electron,
    photon,
    charged_kaon,
    neutral_kaon,
    lambda,
    charged_sigma
};

const std::array<std::string, 16> label_names = {
    "empty", "cosmic", "neutrino", "proton", "neutron", "charged_pion", "neutral_pion", "muon", "electron", "photon",
    "charged_kaon", "neutral_kaon", "lambda", "charged_sigma"
};

inline std::vector<PrimaryLabel> classifyParticles(
    const art::Event& event,
    const art::InputTag& particle_label,
    double gamma_threshold,
    double hadron_threshold
) {
    auto const& particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(particle_label);
    const auto& particles = *particle_handle;

    std::map<int, size_t> id_to_index;
    for (size_t i = 0; i < particles.size(); ++i) {
        id_to_index[particles[i].TrackId()] = i;
    }

    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i].Mother() == 0) {
             if (auto it = track_id_to_particle_index.find(mcp.TrackId()); it != track_id_to_particle_index.end()) {
                PrimaryLabel label = getPrimaryLabel(mcp.PdgCode());
                processPrimaryParticle(it->second, particles, particle_labels, track_id_to_particle_index, label);
            }
        }
    }
    
    return particle_labels;
}

inline void processPrimaryParticle(
    size_t particle_index,
    const std::vector<simb::MCParticle>& particles,
    std::vector<PrimaryLabel>& particle_labels,
    const std::unordered_map<int, size_t>& track_id_to_particle_index,
    PrimaryLabel primary_label
) {
    particle_labels[particle_index] = primary_label; 
    const auto& particle = particles[particle_index];

    for (int daughter_index = 0; daughter_index < particle.NumberDaughters(); ++daughter_index) {
        if (auto it = track_id_to_particle_index.find(particle.Daughter(daughter_index)); it != track_id_to_particle_index.end()) {
            processPrimaryParticle(it->second, particles, particle_labels, track_id_to_particle_index, primary_label);
        }
    }
}

inline PrimaryLabel getPrimaryLabel(int pdg) {
    pdg = std::abs(pdg);
    if (pdg == 2212) return PrimaryLabel::proton;
    if (pdg == 2112) return PrimaryLabel::neutron;
    if (pdg == 211) return PrimaryLabel::charged_pion;
    if (pdg == 111) return PrimaryLabel::neutral_pion;
    if (pdg == 13) return PrimaryLabel::muon;
    if (pdg == 11) return PrimaryLabel::electron;
    if (pdg == 22) return PrimaryLabel::photon;
    if (pdg == 321) return PrimaryLabel::charged_kaon;
    if (pdg == 311) return PrimaryLabel::neutral_kaon;
    if (pdg == 3122) return PrimaryLabel::lambda;
    if (pdg == 3222 || pdg == 3112) return PrimaryLabel::charged_sigma;
    if (pdg == 3212) return PrimaryLabel::neutral_sigma;
    if (pdg == 3322 || pdg == 3312 || pdg == 3334) return PrimaryLabel::other_hyperon;
    return PrimaryLabel::other;
}

}

#endif