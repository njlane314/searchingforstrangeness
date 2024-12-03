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
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Containment.h"

namespace signature {

using Signature = std::vector<art::Ptr<simb::MCParticle>>;

class SignatureToolBase 
{
public:

    virtual ~SignatureToolBase() noexcept = default;
    
    virtual void configure(fhicl::ParameterSet const& pset)
    {
        const std::vector<int> pdg_vec = {211, 13, 2212, 321, 11, 3222, 3112, 14, 2112};
        const std::vector<std::string> part_vec = {"Pion", "Muon", "Proton", "Kaon", "Electron", "SigmaPlus", "SigmaMinus", "NuMu", "Neutron"};
        for (size_t i = 0; i < pdg_vec.size(); ++i) {
            _thresh_map[pdg_vec[i]] = pset.get<float>(part_vec[i] + "Threshold", 0.1);
        }

        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");

        _fv_x_start = pset.get<float>("fidvolXstart", 10.0);
        _fv_y_start = pset.get<float>("fidvolYstart", 15.0);
        _fv_z_start = pset.get<float>("fidvolZstart", 10.0);
        _fv_x_end = pset.get<float>("fidvolXend", 10.0);
        _fv_y_end = pset.get<float>("fidvolYend", 15.0);
        _fv_z_end = pset.get<float>("fidvolZend", 50.0);
    }

    bool identifySignalParticles(art::Event const& evt, Signature& signature)
    {
        auto const& truth_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
        if (truth_handle->size() != 1) 
        {
            signature.clear();
            return false;
        }

        const simb::MCTruth& truth = truth_handle->at(0);
        const TLorentzVector& nu_vertex = truth.GetNeutrino().Nu().Position();
        double vertex[3] = {nu_vertex.X(), nu_vertex.Y(), nu_vertex.Z()};

        if (!common::point_inside_fv(vertex, _fv_x_start, _fv_y_start, _fv_z_start, _fv_x_end, _fv_y_end, _fv_z_end)) 
        {
            signature.clear();
            return false;
        }

        bool found_signature = false;
        this->findSignature(evt, signature, found_signature);

        if (!found_signature)  
            signature.clear();

        return found_signature;
    }

protected:
    std::unordered_map<int, float> _thresh_map;

    art::InputTag _MCPproducer, _MCTproducer;

    float _fv_x_start, _fv_y_start, _fv_z_start;
    float _fv_x_end, _fv_y_end, _fv_z_end;

    bool aboveThreshold(const simb::MCParticle& mcp) const 
    {
        float mom_mag = mcp.Momentum().Vect().Mag();
        int abs_pdg = std::abs(mcp.PdgCode());

        auto it = _thresh_map.find(abs_pdg);
        if (it != _thresh_map.end()) 
            return mom_mag > it->second;

        return false;
    }

    void fillSignature(const art::Ptr<simb::MCParticle>& mcp, Signature& signature) 
    {
        signature.push_back(mcp);
    }

    virtual void findSignature(art::Event const& evt, Signature& signature, bool& found_signature) = 0;
};

} 

#endif
