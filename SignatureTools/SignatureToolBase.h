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
    
    void configure(const fhicl::ParameterSet&){};
    virtual void identifySignalParticles(art::Event const& evt, 
                                         std::vector<Particle>& particle_elements, 
                                         bool& found_signature, 
                                         std::array<float, 3>& nu_vtx) = 0;
};

} 

#endif