#ifndef SIGNATURE_KAONSHORT_CXX
#define SIGNATURE_KAONSHORT_CXX

#include <iostream>
#include <optional>
#include "SignatureToolBase.h"
#include "CommonFunctions/Types.h" 

namespace signature
{
    class KaonShortSignature : public SignatureToolBase, public VertexToolBase {
    public:
        explicit KaonShortSignature(const fhicl::ParameterSet& pset)
            : _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
            , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")} {
            configure(pset);
        }
        
        ~KaonShortSignature() override = default;

        void configure(fhicl::ParameterSet const& pset) override {
            SignatureToolBase::configure(pset);
        }

        SignatureType getSignatureType() const override {
            return kKaonShortSignature;
        }

        TVector3 findVertex(art::Event const& evt) const override {
            auto vertex = findDecayVertex(evt, 310, {211, -211});
            return vertex.value_or(TVector3());
        }

    protected:
        void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) override {
            auto mcp_map = buildParticleMap(evt); 
            auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
            
            for (const auto& mcp : *mcp_h) {
                if (abs(mcp.PdgCode()) == 311 && mcp.Process() == "primary" && 
                    mcp.EndProcess() == "Decay" && !signature_found) {
                    
                    auto neutral_kaon_daughters = common::GetDaughters(mcp_map.at(mcp.TrackId()), mcp_map);
                    for (const auto& kaon_shorts : neutral_kaon_daughters) {
                        if (kaon_shorts->PdgCode() == 310 && kaon_shorts->EndProcess() == "Decay") {
                            if (matchesDecayProducts(kaon_shorts, mcp_map, {211, -211})) {
                                auto decay_products = common::GetDaughters(kaon_shorts, mcp_map);
                                signature_found = true;
                                for (const auto& pion : decay_products) {
                                    fillSignature(pion, signature);
                                    auto addDaughterInteractions = [this, &signature, &mcp_map](
                                        const art::Ptr<simb::MCParticle>& particle, auto& self) -> void {
                                        auto daughters = common::GetDaughters(particle, mcp_map);
                                        for (const auto& daugh : daughters) {
                                            if (daugh->PdgCode() == particle->PdgCode()) {
                                                fillSignature(daugh, signature);
                                                self(daugh, self);
                                            }
                                        }
                                    };
                                    
                                    addDaughterInteractions(pion, addDaughterInteractions);
                                }
                                
                                if (signature_found) return;
                            }
                        }
                    }
                }
            }
        }

    private:
        art::InputTag _HitProducer;
        art::InputTag _BacktrackTag;
    };

    DEFINE_ART_CLASS_TOOL(KaonShortSignature)
} 

#endif