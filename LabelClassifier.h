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

namespace signature 
{
    enum class Label {
        MIP,       
        HIP,       
        shower,    
        michel,    
        diffuse,   
        invisible, 
        undefined,
        cosmic  
    };

    const std::array<std::string, 8> label_names = {
        "undefined", "cosmic", "MIP", "HIP", "shower", "michel", "diffuse", "invisible"
    };

    std::string label_to_string(Label label) {
        size_t index = static_cast<size_t>(label);
        if (index < label_names.size()) {
            return label_names[index];
        }
        return "unknown";
    }

    class LabelClassifier {
    public:
        explicit LabelClassifier(const fhicl::ParameterSet& pset)
            : _particle_label{pset.get<art::InputTag>("particle_label", "largeant")},
              _gamma_threshold{pset.get<double>("gamma_threshold", 0.02)},
              _hadron_threshold{pset.get<double>("hadron_threshold", 0.2)} {}

        std::vector<signature::Label> classifyParticles(const art::Event& event) const {
            auto const& particle_handle = event.getValidHandle<std::vector<simb::MCParticle>>(_particle_label);
            const auto& particles = *particle_handle;

            std::map<int, size_t> id_to_index;
            for (size_t i = 0; i < particles.size(); ++i) {
                id_to_index[particles[i].TrackId()] = i;
            }

            std::vector<signature::Label> particle_labels(particles.size());
            
            for (size_t i = 0; i < particles.size(); ++i) {
                if (particles[i].Mother() == 0) {
                    process_particle(i, particles, Label::undefined, particle_labels, id_to_index);
                }
            }
            return particle_labels;
        }

    private:
        void process_particle(size_t idx, const std::vector<simb::MCParticle>& particles, Label sl_from_parent,
                              std::vector<signature::Label>& particle_labels, const std::map<int, size_t>& id_to_index) const {
            const auto& part = particles[idx];
            auto [sl, slc] = compute_label(particles, part, sl_from_parent, id_to_index);
            particle_labels[idx] = sl;

            for (int i = 0; i < part.NumberDaughters(); ++i) {
                int daughter_id = part.Daughter(i);
                auto it = id_to_index.find(daughter_id);
                if (it != id_to_index.end()) {
                    process_particle(it->second, particles, slc, particle_labels, id_to_index);
                }
            }
        }

        std::pair<Label, Label> compute_label(const std::vector<simb::MCParticle>& particles, const simb::MCParticle& part,
                                              Label sl_from_parent, const std::map<int, size_t>& id_to_index) const {
            if (sl_from_parent != Label::undefined) {
                return {sl_from_parent, sl_from_parent};
            }

            Label sl = Label::invisible;  
            Label slc = Label::undefined; 

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
                sl = Label::MIP;
            }
            else if (pdg == 321 || pdg == -321 || (std::abs(pdg) == 2212 && momentum >= _hadron_threshold)) { 
                sl = Label::HIP;
            }
            else if (pdg == 11 || pdg == -11) { 
                if (start_process == "primary") {
                    sl = Label::shower;
                    slc = Label::shower;
                }
                else if (std::abs(parent_pdg) == 13 && (start_process == "muMinusCaptureAtRest" || 
                                                        start_process == "muPlusCaptureAtRest" || 
                                                        start_process == "Decay")) {
                    sl = Label::michel;
                    slc = Label::michel;
                }
                else if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                    if (momentum >= _gamma_threshold) {
                        sl = Label::shower;
                        slc = Label::shower;
                    }
                    else {
                        sl = Label::diffuse;
                    }
                }
                else {
                    sl = Label::diffuse;
                }
            }
            else if (pdg == 22) { 
                if (start_process == "conv" || end_process == "conv" || start_process == "compt" || end_process == "compt") {
                    if (momentum >= _gamma_threshold) {
                        sl = Label::shower;
                        slc = Label::shower;
                    }
                    else {
                        sl = Label::diffuse;
                    }
                }
                else {
                    sl = Label::diffuse;
                }
            }
            else if (std::abs(pdg) == 2212 && momentum < _hadron_threshold) { 
                sl = Label::diffuse;
            }

            return {sl, slc};
        }

        art::InputTag _particle_label;  
        double _gamma_threshold;        
        double _hadron_threshold;      
    };
}

#endif