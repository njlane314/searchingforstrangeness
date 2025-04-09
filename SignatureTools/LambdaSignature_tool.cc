#ifndef SIGNATURE_LAMBDA_CXX
#define SIGNATURE_LAMBDA_CXX

#include <iostream>
#include "SignatureToolBase.h"

namespace signature
{
    class LambdaSignature : public SignatureToolBase, public VertexToolBase {
    public:
        explicit LambdaSignature(const fhicl::ParameterSet& pset) {
            configure(pset);
        }
        
        ~LambdaSignature() override = default;
        
        void configure(fhicl::ParameterSet const& pset) override {
            SignatureToolBase::configure(pset);
        }
        
        SignatureType getSignatureType() const override {
            return kLambdaSignature;
        }
        
        TVector3 findVertex(art::Event const& evt) const override {
            auto vertex = findDecayVertex(evt, 3122, {-211, 2212});
            return vertex.value_or(TVector3());
        }
        
    protected:
        void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override;
    };

    void LambdaSignature::findSignature(art::Event const& evt, Signature& signature, bool& signature_found) {
        auto mcp_map = buildMCParticleMap(evt);
        auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
        
        for (const auto& mcp : *mcp_h) {
            if (abs(mcp.PdgCode()) == 3122 && mcp.Process() == "primary" && 
                mcp.EndProcess() == "Decay" && !signature_found) {
                
                if (matchesDecayProducts(mcp_map.at(mcp.TrackId()), mcp_map, {-211, 2212})) {
                    auto decay_products = common::GetDaughters(mcp_map.at(mcp.TrackId()), mcp_map);
                    signature_found = true;
                    
                    for (const auto& particle : decay_products) {
                        const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(particle->PdgCode());
                        if (info->Charge() == 0.0)
                            continue;
                            
                        fillSignature(particle, signature);
                        
                        traverseDecayChain(particle, mcp_map, [this, &signature, particle_pdg = particle->PdgCode()]
                                          (const art::Ptr<simb::MCParticle>& p) {
                            if (p->PdgCode() == particle_pdg) {
                                const TParticlePDG* info = TDatabasePDG::Instance()->GetParticle(p->PdgCode());
                                if (info && info->Charge() != 0.0) {
                                    fillSignature(p, signature);
                                }
                            }
                        });
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(LambdaSignature)
}

#endif