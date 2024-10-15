#ifndef SIGNATURETOOLBASE_H
#define SIGNATURETOOLBASE_H

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "CommonFunctions/Types.h"

#include "TTree.h"
#include <limits>

namespace signature {

struct Particle
{
    int pdg_code;
    int track_id; 
    std::array<float, 3> momentum;
};

class SignatureToolBase {

public:

    virtual ~SignatureToolBase() noexcept = default;
    
    void configure(const fhicl::ParameterSet&){
        _thresh_map[211] = pset.get<float>("PionThreshold", 0.1);
        _thresh_map[13] = pset.get<float>("MuonThreshold", 0.1);
        _thresh_map[2212] = pset.get<float>("ProtonThreshold", 0.1);
        _thresh_map[321] = pset.get<float>("KaonThreshold", 0.1);
    };

    bool identifySignalParticles(art::Event const& evt, std::vector<Particle>& particle_elements)
    {
        bool found_signature = false;
        particle_elements.clear();

        this->findSignature(evt, particle_elements, found_signature);

        if (!found_signature)  
            particle_elements.clear();

        return found_signature;
    }

    bool aboveThreshold(const simb::MCParticle& mc_particle) const {
        float mom_mag = mc_particle.Momentum().Vect().Mag()
        int abs_pdg = std::abs(mc_particle.PdgCode());

        if (it != _thresh_map.end()) {
            return mom_mag > it->second;
        }

        return false;
    }

protected:
    std::unordered_map<int, float> _thresh_map;

    void fillParticle(const art::Ptr<simb::MCParticle>& mc_particle, std::vector<Particle>& particle_elements) {
        Particle particle;
        particle.pdg_code = mc_particle->PdgCode();
        particle.track_id = mc_particle->TrackId();
        particle.momentum = {static_cast<float>(mc_particle->Px()), static_cast<float>(mc_particle->Py()), static_cast<float>(mc_particle->Pz())};
        particle_elements.push_back(particle);
    }

    virtual void findSignature(art::Event const& evt, std::vector<Particle>& particle_elements, bool& found_signature) = 0;
};

} 

#endif