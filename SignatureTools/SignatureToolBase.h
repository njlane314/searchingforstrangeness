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
        const std::vector<int> particle_ids = {211, 13, 2212, 321, 11, 3222, 3112, 3312, 3334};
        const std::vector<std::string> particle_names = {"Pion", "Muon", "Proton", "Kaon", "Electron", 
                                                         "SigmaPlus", "SigmaMinus", "Xi", "Omega"};
        for (size_t i = 0; i < particle_ids.size(); ++i) {
            _thresh_map[particle_ids[i]] = pset.get<float>(particle_names[i] + "Threshold", 0.1);
        }

        _include_mesons = pset.get<bool>("IncludeMesons", true);
        _include_hyperons = pset.get<bool>("IncludeHyperons", true);

        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");

        _fv_x_start = pset.get<float>("FvXStart", 10.0);
        _fv_y_start = pset.get<float>("FvYStart", 10.0);
        _fv_z_start = pset.get<float>("FvZStart", 10.0);
        _fv_x_end = pset.get<float>("FvXEnd", 10.0);
        _fv_y_end = pset.get<float>("FvYEnd", 10.0);
        _fv_z_end = pset.get<float>("FvZEnd", 10.0);
    }

    bool identifySignalParticles(art::Event const& evt, TraceCollection& trace_coll)
    {
        auto const& truth_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
        if (truth_handle->size() != 1) 
        {
            mf::LogWarning("Skipping event with more than one truth object.");
            trace_coll.clear();
            return false;
        }

        const simb::MCTruth& truth = truth_handle->at(0);
        const TLorentzVector& nu_vertex = truth.GetNeutrino().Nu().Position();
        double vertex[3] = {nu_vertex.X(), nu_vertex.Y(), nu_vertex.Z()};

        if (!common::point_inside_fv(vertex, _fv_x_start, _fv_y_start, _fv_z_start, _fv_x_end, _fv_y_end, _fv_z_end)) 
        {
            mf::LogWarning("Neutrino interaction vertex is outside the fiducial volume.");
            trace_coll.clear();
            return false;
        }

        bool found_signature = false;
        this->findSignature(evt, trace_coll, found_signature);

        if (!found_signature)  
        {
            trace_coll.clear();
            return false;
        }   

        if ((!_include_mesons && hasAdditionalParticles(evt, trace_coll, [&](int pdgCode) { return isChargedMeson(pdgCode); })) ||
            (!_include_hyperons && hasAdditionalParticles(evt, trace_coll, [&](int pdgCode) { return isChargedHyperon(pdgCode); }))) 
        {
            trace_coll.clear();
            return false;
        }

        return found_signature;
    }

protected:
    std::unordered_map<int, float> _thresh_map;
    bool _include_mesons;
    bool _include_hyperons;

    art::InputTag _MCPproducer, _MCTproducer;

    float _fv_x_start, _fv_y_start, _fv_z_start;
    float _fv_x_end, _fv_y_end, _fv_z_end;

    bool aboveThreshold(const simb::MCParticle& mc_particle) const 
    {
        float mom_mag = mc_particle.Momentum().Vect().Mag();
        int abs_pdg = std::abs(mc_particle.PdgCode());

        auto it = _thresh_map.find(abs_pdg);
        if (it != _thresh_map.end()) {
            return mom_mag > it->second;
        }

        return mom_mag > 0.1; 
    }

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
        std::cout << trace_coll.size() << std::endl;
    }

    virtual void findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature) = 0;
};

} 

#endif