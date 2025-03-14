#ifndef SIGNATURE_TOOLBASE_H
#define SIGNATURE_TOOLBASE_H

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

enum SignatureType {
    SignatureInvalid = -1,
    SignatureEmpty = 0,
    SignatureNoise,
    SignaturePrimaryMuon,
    SignatureChargedKaon,
    SignatureKaonShort,
    SignatureLambda,
    SignatureChargedSigma
};

using Signature = std::pair<int,std::vector<art::Ptr<simb::MCParticle>>>;
using Pattern = std::vector<Signature>;

std::string GetSignatureName(const Signature& sig)
{

  switch (sig.first){
        case SignatureInvalid: return "Invalid";
        case SignatureEmpty: return "Empty";
        case SignatureNoise: return "Noise";
        case SignaturePrimaryMuon: return "PrimaryMuon";
        case SignatureChargedKaon: return "ChargedKaon";
        case SignatureKaonShort: return "KaonShort";
        case SignatureLambda: return "Lambda";
        case SignatureChargedSigma: return "ChargedSigma";
        default: return "Invalid"; 
  }

  return "Invalid";

}

class SignatureToolBase 
{
public:

    virtual ~SignatureToolBase() noexcept = default;
    
    virtual void configure(fhicl::ParameterSet const& pset)
    {
        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
    }

    bool constructSignature(art::Event const& evt, Signature& signature)
    {
        signature.second.clear();
        auto const& truth_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
        if (truth_handle->size() != 1) 
            return false;

        bool signature_found = false;
        this->findSignature(evt, signature, signature_found);

        if (!signature_found){
          //signature.clear();
          signature.second.clear();
          signature.first = SignatureInvalid;
        }

        return signature_found;
    }

protected:
    art::InputTag _MCPproducer, _MCTproducer;

    bool assessParticle(const simb::MCParticle& mcp) const 
    {
        float mom_mag = mcp.Momentum().Vect().Mag();
        int abs_pdg = std::abs(mcp.PdgCode());

        const TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(mcp.PdgCode());
        if (particle->Charge() == 0.0) 
            return true;

        std::unordered_map<int, float> thresh_map = {
            {211, 0.1},    // pi
            {13, 0.1},     // mu
            {2212, 0.3},   // p
            {321, 0.2},    // K
            {11, 0.1},     // e
            {3222, 0.1},   // sigma+
            {3112, 0.1},   // sigma-
        };
       
        // Check start is inside the TPC
        double pos[3] = {(double)mcp.Vx(),(double)mcp.Vy(),(double)mcp.Vz()};
        if(!common::point_inside_fv(pos)) return false;

        auto it = thresh_map.find(abs_pdg);
        if (it != thresh_map.end())
            return mom_mag > it->second;

        return false;
    }

    void fillSignature(const art::Ptr<simb::MCParticle>& mcp, Signature& signature) 
    {
        signature.second.push_back(mcp);
    }

    virtual void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) = 0;
};

} 

#endif
