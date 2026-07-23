#include "ImagePipeline/SemanticClassifier.h"

#include <cstdlib>

namespace image {

SemanticClassifier::SemanticClassifier(
    const art::InputTag &mcp_producer)
    : mcp_producer_{mcp_producer}
{}

SemanticClassifier::SemanticLabel
SemanticClassifier::getSemanticLabel(int pdg) {
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

void SemanticClassifier::propagateLabel(
    std::size_t particle_index,
    const std::vector<simb::MCParticle> &particles,
    std::vector<SemanticLabel> &particle_labels,
    const std::unordered_map<int, std::size_t> &
        track_id_to_index,
    SemanticLabel primary_label_to_assign) {
    std::vector<unsigned char> visited(
        particles.size(), 0U);
    propagateLabelImpl(
        particle_index, particles, particle_labels,
        track_id_to_index, primary_label_to_assign, visited);
}

void SemanticClassifier::propagateLabelImpl(
    std::size_t particle_index,
    const std::vector<simb::MCParticle> &particles,
    std::vector<SemanticLabel> &particle_labels,
    const std::unordered_map<int, std::size_t> &
        track_id_to_index,
    SemanticLabel primary_label_to_assign,
    std::vector<unsigned char> &visited) {
    if (particle_index >= particles.size() ||
        particle_index >= particle_labels.size() ||
        particle_index >= visited.size() ||
        visited[particle_index] != 0U) {
        return;
    }

    visited[particle_index] = 1U;
    particle_labels[particle_index] =
        primary_label_to_assign;
    auto const &particle = particles[particle_index];

    for (int daughter = 0;
         daughter < particle.NumberDaughters(); ++daughter) {
        int const track_id = particle.Daughter(daughter);
        auto const found = track_id_to_index.find(track_id);
        if (found != track_id_to_index.end()) {
            propagateLabelImpl(
                found->second, particles, particle_labels,
                track_id_to_index, primary_label_to_assign,
                visited);
        }
    }
}

std::vector<SemanticClassifier::SemanticLabel>
SemanticClassifier::classifyParticles(
    const art::Event &event) const {
    auto const particles_handle =
        event.getValidHandle<std::vector<simb::MCParticle>>(
            mcp_producer_);
    auto const &particles = *particles_handle;

    std::unordered_map<int, std::size_t>
        track_id_to_index;
    for (std::size_t i = 0U; i < particles.size(); ++i)
        track_id_to_index[particles[i].TrackId()] = i;

    std::vector<SemanticLabel> labels(
        particles.size(), SemanticLabel::Empty);
    for (std::size_t i = 0U; i < particles.size(); ++i) {
        if (particles[i].Mother() != 0)
            continue;

        auto const found =
            track_id_to_index.find(particles[i].TrackId());
        if (found == track_id_to_index.end())
            continue;

        propagateLabel(
            found->second, particles, labels,
            track_id_to_index,
            getSemanticLabel(particles[i].PdgCode()));
    }
    return labels;
}

} // namespace image
