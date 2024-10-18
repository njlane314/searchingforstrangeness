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

namespace signature {

struct Trace
{
    int pdg;
    int trckid; 
    std::array<float, 3> mom;
};

using TraceCollection = std::vector<Trace>&;

class SignatureToolBase 
{
    
public:

    virtual ~SignatureToolBase() noexcept = default;
    
    virtual void configure(fhicl::ParameterSet const& pset)
    {
        _thresh_map[211] = pset.get<float>("PionThreshold", 0.1);
        _thresh_map[13] = pset.get<float>("MuonThreshold", 0.1);
        _thresh_map[2212] = pset.get<float>("ProtonThreshold", 0.1);
        _thresh_map[321] = pset.get<float>("KaonThreshold", 0.1);
        _thresh_map[11] = pset.get<float>("ElectronThreshold", 0.1);

        _thresh_map[3222] = pset.get<float>("SigmaPlusThreshold", 0.1);     
        _thresh_map[3112] = pset.get<float>("SigmaMinusThreshold", 0.1);   
        _thresh_map[3312] = pset.get<float>("XiThreshold", 0.1);       
        _thresh_map[3334] = pset.get<float>("OmegaThreshold", 0.1); 

        _include_mesons = pset.get<bool>("IncludeMesons", true);
        _include_hyperons = pset.get<bool>("IncludeHyperons", true);

        _MCPproducer = pset.get<art::InputTag>("MCProducerTag", "largeant");
    };

    bool identifySignalParticles(art::Event const& evt, TraceCollection& trace_coll)
    {
        bool found_signature = false;
        trace_coll.clear();

        // add only single mctruth condition
        // and requirement on fiducial volume

        this->findSignature(evt, trace_coll, found_signature);

        if (!found_signature)  
        {
            trace_coll.clear();
            return false;
        }   

        if (!_include_mesons && this->hasAdditionalParticles(evt, trace_coll, [&](int pdg_code) {
                return isChargedMeson(pdg_code);
            })) 
        {
            trace_coll.clear();
            return false;
        }

        if (!_include_hyperons && this->hasAdditionalParticles(evt, trace_coll, [&](int pdg_code) {
                return isChargedHyperon(pdg_code);
            })) 
        {
            trace_coll.clear();
            return false;
        }

        return found_signature;
    }

    bool aboveThreshold(const simb::MCParticle& mc_particle) const 
    {
        std::cout << "Checking threshold..." << std::endl;
        float mom_mag = mc_particle.Momentum().Vect().Mag();
        int abs_pdg = std::abs(mc_particle.PdgCode());

        auto it = _thresh_map.find(abs_pdg);
        if (it != _thresh_map.end()) {
            std::cout << "Momentum " << mom_mag << " and threshold " << it->second << std::endl;
            return mom_mag > it->second;
        }

        return mom_mag > 0.1; 
    }

protected:
    std::unordered_map<int, float> _thresh_map;
    bool _include_mesons;
    bool _include_hyperons;

    art::InputTag _MCPproducer;

    bool isChargedMeson(int pdg_code) const 
    {
        int abs_pdg = std::abs(pdg_code);
        int thousands_digit = (abs_pdg / 1000) % 10;
        int hundreds_digit = (abs_pdg / 100) % 10;

        if (thousands_digit != 0) return false;
        if (hundreds_digit == 0) return false;

        TParticlePDG* particle_info = TDatabasePDG::Instance()->GetParticle(pdg_code);
        if (!particle_info) 
            return false;

        return (particle_info->Charge() != 0);
    }

    bool isChargedHyperon(int pdg_code) 
    {   
        TParticlePDG* particle_info = TDatabasePDG::Instance()->GetParticle(pdg_code);
        if (!particle_info) 
            return false;

        int abs_pdg = std::abs(pdg_code);

        int thousands_digit = (abs_pdg / 1000) % 10;
        int hundreds_digit = (abs_pdg / 100) % 10;

        bool is_hyperon = 
            (thousands_digit == 3 && (hundreds_digit == 1 || hundreds_digit == 2)) ||  // Strange hyperons
            (thousands_digit == 4 && (hundreds_digit >= 1 && hundreds_digit <= 3));    // Charmed hyperons

        return is_hyperon && (particle_info->Charge() != 0);
    }

    template<typename TraceFilter>
    bool hasAdditionalParticles(art::Event const& evt, const TraceCollection& trace_coll, TraceFilter trace_filter) const 
    {
        std::cout << "Running additional particles function..." << std::endl;
        auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);  
        std::set<int> primary_track_ids;

        for (const auto& particle : trace_coll) 
            primary_track_ids.insert(particle.trckid);

        for (const auto& mc_particle : *mcp_h) 
        {
            if (trace_filter(mc_particle.PdgCode()) && primary_track_ids.find(mc_particle.TrackId()) == primary_track_ids.end() && aboveThreshold(mc_particle))
                return true;
        }

        return false;  
    }

    void fillTrace(const art::Ptr<simb::MCParticle>& mc_particle, TraceCollection& trace_coll) 
    {
        Trace trace;
        trace.pdg = mc_particle->PdgCode();
        trace.trckid = mc_particle->TrackId();
        trace.mom = {static_cast<float>(mc_particle->Px()), static_cast<float>(mc_particle->Py()), static_cast<float>(mc_particle->Pz())};
        trace_coll.push_back(trace);
    }

    virtual void findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature) = 0;
};

} 

#endif