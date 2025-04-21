#ifndef NEUTRINOPRIMARYPARTICLECLASSIFIER_H
#define NEUTRINOPRIMARYPARTICLECLASSIFIER_H

#include <vector>
#include <map>
#include <string>
#include <array>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace signature {

enum class PrimaryLabel {
    undefined,
    proton,
    neutron,
    charged_pion,
    neutral_pion,
    muon,
    electron,
    photon,
    kaon,
    hyperon,
    other
};

const std::array<std::string, 11> label_names = {
    "undefined", "proton", "neutron", "charged_pion", "neutral_pion", "muon", "electron", "photon", "kaon", "hyperon", "other"
};

std::string label_to_string(PrimaryLabel label) {
    size_t index = static_cast<size_t>(label);
    if (index < label_names.size()) {
        return label_names[index];
    }
    return "unknown";
}

class NeutrinoPrimaryParticleClassifier {
public:
    explicit NeutrinoPrimaryParticleClassifier(const fhicl::ParameterSet& pset)
        : _mctruth_label{pset.get<art::InputTag>("mctruth_label", "generator")},
          _particle_label{pset.get<art::InputTag>("particle_label", "largeant")} {}

    std::vector<PrimaryLabel> classifyParticles(const art::Event& event) const;

private:
    void process_particle(size_t idx, const std::vector<simb::MCParticle>& particles,
                          std::vector<PrimaryLabel>& particle_labels,
                          const std::map<int, size_t>& id_to_index) const;

    PrimaryLabel get_label(int pdg) const;

    art::InputTag _mctruth_label;
    art::InputTag _particle_label;
};

std::vector<PrimaryLabel> NeutrinoPrimaryParticleClassifier::classifyParticles(const art::Event& event) const {
    auto const& particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(_particle_label);
    const auto& particles = *particle_handle;

    std::map<int, size_t> id_to_index;
    for (size_t i = 0; i < particles.size(); ++i) {
        id_to_index[particles[i].TrackId()] = i;
    }

    std::vector<PrimaryLabel> particle_labels(particles.size(), PrimaryLabel::undefined);

    auto const& mctruth_handle = event.getValidHandle<std::vector<simb::MCTruth>>(_mctruth_label);
    for (const auto& mctruth : *mctruth_handle) {
        if (mctruth.NeutrinoSet()) {
            for (int i = 0; i < mctruth.NParticles(); ++i) {
                const auto& mcp = mctruth.GetParticle(i);
                if (mcp.StatusCode() == 1) {
                    auto it = id_to_index.find(mcp.TrackId());
                    if (it != id_to_index.end()) {
                        process_particle(it->second, particles, particle_labels, id_to_index);
                    }
                }
            }
        }
    }
    return particle_labels;
}

void NeutrinoPrimaryParticleClassifier::process_particle(size_t idx,
    const std::vector<simb::MCParticle>& particles,
    std::vector<PrimaryLabel>& particle_labels,
    const std::map<int, size_t>& id_to_index) const {
    const auto& part = particles[idx];
    particle_labels[idx] = get_label(part.PdgCode());

    for (int i = 0; i < part.NumberDaughters(); ++i) {
        int daughter_id = part.Daughter(i);
        auto it = id_to_index.find(daughter_id);
        if (it != id_to_index.end()) {
            process_particle(it->second, particles, particle_labels, id_to_index);
        }
    }
}

PrimaryLabel NeutrinoPrimaryParticleClassifier::get_label(int pdg) const {
    switch (pdg) {
        case 2212: return PrimaryLabel::proton;
        case 2112: return PrimaryLabel::neutron;
        case 211:
        case -211: return PrimaryLabel::charged_pion;
        case 111: return PrimaryLabel::neutral_pion;
        case 13:
        case -13: return PrimaryLabel::muon;
        case 11:
        case -11: return PrimaryLabel::electron;
        case 22: return PrimaryLabel::photon;
        case 321:
        case -321:
        case 311:
        case -311: return PrimaryLabel::kaon;
        case 3122:
        case -3122:
        case 3222:
        case -3222:
        case 3212:
        case -3212:
        case 3112:
        case -3112:
        case 3322:
        case -3322:
        case 3312:
        case -3312:
        case 3334:
        case -3334: return PrimaryLabel::hyperon;
        default: return PrimaryLabel::other;
    }
}

}

#endif