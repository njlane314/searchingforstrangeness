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
    undefined,
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
    charged_sigma,
    neutral_sigma,
    other_hyperon,
    other
};

const std::array<std::string, 15> label_names = {
    "undefined", "proton", "neutron", "charged_pion", "neutral_pion", "muon", "electron", "photon",
    "charged_kaon", "neutral_kaon", "lambda", "charged_sigma", "neutral_sigma", "other_hyperon", "other"
};

inline std::string label_to_string(PrimaryLabel label) {
    size_t index = static_cast<size_t>(label);
    if (index < label_names.size()) {
        return label_names[index];
    }
    return "unknown";
}

inline std::vector<PrimaryLabel> classifyPrimaryParticles(
    const art::Event& event,
    const art::InputTag& mctruth_label,
    const art::InputTag& particle_label
) {
    auto const& particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(particle_label);
    const auto& particles = *particle_handle;

    std::map<int, size_t> id_to_index;
    for (size_t i = 0; i < particles.size(); ++i) {
        id_to_index[particles[i].TrackId()] = i;
    }

    std::vector<PrimaryLabel> particle_labels(particles.size(), PrimaryLabel::undefined);

    auto const& mctruth_handle = event.getValidHandle<std::vector<simb::MCTruth>>(mctruth_label);
    for (const auto& mctruth : *mctruth_handle) {
        if (mctruth.NeutrinoSet()) {
            for (int i = 0; i < mctruth.NParticles(); ++i) {
                const auto& mcp = mctruth.GetParticle(i);
                if (mcp.StatusCode() == 1) {
                    auto it = id_to_index.find(mcp.TrackId());
                    if (it != id_to_index.end()) {
                        process_primary_particle(it->second, particles, particle_labels, id_to_index);
                    }
                }
            }
        }
    }
    return particle_labels;
}

inline void process_primary_particle(
    size_t idx,
    const std::vector<simb::MCParticle>& particles,
    std::vector<PrimaryLabel>& particle_labels,
    const std::map<int, size_t>& id_to_index
) {
    const auto& part = particles[idx];
    particle_labels[idx] = get_primary_label(part.PdgCode());

    for (int i = 0; i < part.NumberDaughters(); ++i) {
        int daughter_id = part.Daughter(i);
        auto it = id_to_index.find(daughter_id);
        if (it != id_to_index.end()) {
            process_primary_particle(it->second, particles, particle_labels, id_to_index);
        }
    }
}

inline PrimaryLabel get_primary_label(int pdg) {
    if (pdg == 2212) return PrimaryLabel::proton;
    if (pdg == 2112) return PrimaryLabel::neutron;
    if (std::abs(pdg) == 211) return PrimaryLabel::charged_pion;
    if (pdg == 111) return PrimaryLabel::neutral_pion;
    if (std::abs(pdg) == 13) return PrimaryLabel::muon;
    if (std::abs(pdg) == 11) return PrimaryLabel::electron;
    if (pdg == 22) return PrimaryLabel::photon;
    if (std::abs(pdg) == 321) return PrimaryLabel::charged_kaon;
    if (std::abs(pdg) == 311) return PrimaryLabel::neutral_kaon;
    if (std::abs(pdg) == 3122) return PrimaryLabel::lambda;
    if (std::abs(pdg) == 3222 || std::abs(pdg) == 3112) return PrimaryLabel::charged_sigma;
    if (std::abs(pdg) == 3212) return PrimaryLabel::neutral_sigma;
    if (std::abs(pdg) == 3322 || std::abs(pdg) == 3312 || std::abs(pdg) == 3334) return PrimaryLabel::other_hyperon;
    return PrimaryLabel::other;
}

}

#endif