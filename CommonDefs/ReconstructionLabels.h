#ifndef SIMPLELABELCLASSIFIER_H
#define SIMPLELABELCLASSIFIER_H

#include <vector>
#include <map>
#include <string>
#include <array>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace reco_labels 
{
    enum class ReconstructionLabel {
        empty,
        cosmic,  
        MIP,       
        HIP,       
        shower,    
        michel,    
        diffuse,   
        invisible, 
    };

    const std::array<std::string, 8> label_names = {
        "empty", "cosmic", "MIP", "HIP", "shower", "michel", "diffuse", "invisible"
    };

    inline std::string label_to_string(ReconstructionLabel label) {
        size_t index = static_cast<size_t>(label);
        if (index < label_names.size()) {
            return label_names[index];
        }
        return "unknown";
    }

    inline std::pair<ReconstructionLabel, ReconstructionLabel> compute_label(
        const std::vector<simb::MCParticle>& particles,
        const simb::MCParticle& part,
        ReconstructionLabel sl_from_parent,
        const std::map<int, size_t>& id_to_index,
        double gamma_threshold,
        double hadron_threshold
    ) {
        if (sl_from_parent != ReconstructionLabel::empty) {
            return {sl_from_parent, sl_from_parent};
        }

        ReconstructionLabel sl = ReconstructionLabel::invisible;  
        ReconstructionLabel slc = ReconstructionLabel::empty; 

        int pdg = part.PdgCode();
        double momentum = part.P();
        std::string start_process = part.Process();
        std::string end_process = part.EndProcess();
        int parent_id = part.Mother();
        int parent_pdg = 0;
        if (parent_id != 0) {
            auto it = id_to_index.find(parent_id);
            if (it != id_to_index.end()) {
                parent_pdg = particles[it->second].PdgCode();
            }
        }

        if (pdg == 211 || pdg == -211 || pdg == 13 || pdg == -13) { 
            sl = ReconstructionLabel::MIP;
        }
        else if (pdg == 321 || pdg == -321 || (std::abs(pdg) == 2212 && momentum >= hadron_threshold)) { 
            sl = ReconstructionLabel::HIP;
        }
        else if (pdg == 11 || pdg == -11) { 
            if (start_process == "primary") {
                sl = ReconstructionLabel::shower;
                slc = ReconstructionLabel::shower;
            }
            else if (std::abs(parent_pdg) == 13 && (start_process == "muMinusCaptureAtRest" || 
                                                    start_process == "muPlusCaptureAtRest" || 
                                                    start_process == "Decay")) {
                sl = ReconstructionLabel::michel;
                slc = ReconstructionLabel::michel;
            }
            else if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                if (momentum >= gamma_threshold) {
                    sl = ReconstructionLabel::shower;
                    slc = ReconstructionLabel::shower;
                }
                else {
                    sl = ReconstructionLabel::diffuse;
                }
            }
            else {
                sl = ReconstructionLabel::diffuse;
            }
        }
        else if (pdg == 22) { 
            if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                if (momentum >= gamma_threshold) {
                    sl = ReconstructionLabel::shower;
                    slc = ReconstructionLabel::shower;
                }
                else {
                    sl = ReconstructionLabel::diffuse;
                }
            }
            else {
                sl = ReconstructionLabel::diffuse;
            }
        }
        else if (std::abs(pdg) == 2212 && momentum < hadron_threshold) { 
            sl = ReconstructionLabel::diffuse;
        }

        return {sl, slc};
    }

    inline void process_particle(
        size_t idx,
        const std::vector<simb::MCParticle>& particles,
        ReconstructionLabel sl_from_parent,
        std::vector<ReconstructionLabel>& particle_labels,
        const std::map<int, size_t>& id_to_index,
        double gamma_threshold,
        double hadron_threshold
    ) {
        const auto& part = particles[idx];
        auto [sl, slc] = compute_label(particles, part, sl_from_parent, id_to_index, gamma_threshold, hadron_threshold);
        particle_labels[idx] = sl;

        for (int i = 0; i < part.NumberDaughters(); ++i) {
            int daughter_id = part.Daughter(i);
            auto it = id_to_index.find(daughter_id);
            if (it != id_to_index.end()) {
                process_particle(it->second, particles, slc, particle_labels, id_to_index, gamma_threshold, hadron_threshold);
            }
        }
    }

    inline std::vector<ReconstructionLabel> classifyParticles(
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

        std::vector<ReconstructionLabel> particle_labels(particles.size(), ReconstructionLabel::empty);

        for (size_t i = 0; i < particles.size(); ++i) {
            if (particles[i].Mother() == 0) {
                process_particle(i, particles, ReconstructionLabel::empty, particle_labels, id_to_index, gamma_threshold, hadron_threshold);
            }
        }
        return particle_labels;
    }
}

#endif