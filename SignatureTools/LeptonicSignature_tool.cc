#ifndef LEPTONICSIGNATURE_CXX
#define LEPTONICSIGNATURE_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature {

class LeptonicSignature : public SignatureToolBase 
{

public:
    LeptonicSignature(const fhicl::ParameterSet& pset);
    ~LeptonicSignature() {}

    void configure(fhicl::ParameterSet const& pset) override;

protected:
    void findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature) override;

private:
    std::string _lepton_type;  
    bool _is_charged_current; 
    art::InputTag _MCTproducer;  
};

LeptonicSignature::LeptonicSignature(const fhicl::ParameterSet& pset)
    : _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _lepton_type{pset.get<std::string>("Lepton", "muon")}
    , _is_charged_current{pset.get<bool>("ChargedCurrent", true)}
{
    configure(pset); 
}

void LeptonicSignature::configure(fhicl::ParameterSet const& pset)
{
    SignatureToolBase::configure(pset);
}

void LeptonicSignature::findSignature(art::Event const& evt, TraceCollection& trace_coll, bool& found_signature)
{
    found_signature = false;

    auto const& mct_h = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    const simb::MCTruth& mct = mct_h->at(0);
    if (!mct.NeutrinoSet()) return;

    auto const& neutrino = mct.GetNeutrino();

    if (_is_charged_current && neutrino.CCNC() == simb::kNC) return;
    if (!_is_charged_current && neutrino.CCNC() == simb::kCC) return;

    for (const auto& mc_particle : *mcp_h) 
    {
        int pdg_code = std::abs(mc_particle.PdgCode());

        if ((_lepton_type == "muon" && pdg_code == 13) ||
            (_lepton_type == "electron" && pdg_code == 11)) 
        {
            if (this->aboveThreshold(mc_particle) && mc_particle.Process() == "primary") 
            {
                art::Ptr<simb::MCParticle> mc_particle_ptr(mcp_h, &mc_particle - &(*mcp_h)[0]);
                this->fillTrace(mc_particle_ptr, trace_coll);  
                found_signature = true;

                break;  
            }
        }

        if (_lepton_type == "none" && !_is_charged_current) 
        {
            found_signature = true;  

            break;
        }
    }
}

DEFINE_ART_CLASS_TOOL(LeptonicSignature)

} 

#endif