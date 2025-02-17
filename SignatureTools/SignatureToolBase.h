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

enum Type {
    SignatureEmpty = 0,
    SignatureNoise,
    SignaturePrimaryMuon,
    SignatureChargedKaon,
    SignatureKaonShort,
    SignatureLambda,
    SignatureChargedSigma
};

using Signature = std::vector<art::Ptr<simb::MCParticle>>;
using Pattern = std::vector<std::pair<Type, Signature>>;

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
        signature.clear();
        auto const& truth_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
        if (truth_handle->size() != 1) 
            return false;

        bool signature_found = false;
        this->findSignature(evt, signature, signature_found);

        if (!signature_found)  
            signature.clear();

        return signature_found;
    }

    virtual Type getSignatureType() const {
        return SignatureEmpty;
    }

protected:
    art::InputTag _MCPproducer, _MCTproducer;

    void fillSignature(const art::Ptr<simb::MCParticle>& mcp, Signature& signature) 
    {
        signature.push_back(mcp);
    }

    virtual void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) = 0;
};

} 

#endif
