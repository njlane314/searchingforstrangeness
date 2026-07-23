#ifndef SEMANTICCLASSIFIER_H
#define SEMANTICCLASSIFIER_H

#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
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

    static inline const std::array<std::string, 16>
        semantic_label_names = {
            "Empty",
            "Cosmic",
            "Muon",
            "Electron",
            "Photon",
            "ChargedPion",
            "NeutralPion",
            "Neutron",
            "Proton",
            "ChargedKaon",
            "NeutralKaon",
            "Lambda",
            "ChargedSigma",
            "NeutralSigma",
            "Other",
            "Ambiguous"};

    explicit SemanticClassifier(
        const art::InputTag &mcp_producer);

    static SemanticLabel getSemanticLabel(int pdg);
    static void propagateLabel(
        std::size_t particle_index,
        const std::vector<simb::MCParticle> &particles,
        std::vector<SemanticLabel> &particle_labels,
        const std::unordered_map<int, std::size_t> &
            track_id_to_index,
        SemanticLabel primary_label_to_assign);
    std::vector<SemanticLabel>
    classifyParticles(const art::Event &event) const;

  private:
    static void propagateLabelImpl(
        std::size_t particle_index,
        const std::vector<simb::MCParticle> &particles,
        std::vector<SemanticLabel> &particle_labels,
        const std::unordered_map<int, std::size_t> &
            track_id_to_index,
        SemanticLabel primary_label_to_assign,
        std::vector<unsigned char> &visited);

    art::InputTag mcp_producer_;
};

} // namespace image

#endif
