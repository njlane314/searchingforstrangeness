#ifndef SIGNATURE_TOOLBASE_H
#define SIGNATURE_TOOLBASE_H

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "../CommonFunctions/Types.h"
#include "TTree.h"
#include <limits>
#include <vector>
#include <map>
#include <functional>
#include <optional>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "../CommonFunctions/Pandora.h"
#include "../CommonFunctions/Scatters.h"
#include "../CommonFunctions/Corrections.h"
#include "../CommonFunctions/Containment.h"
#include "TVector3.h"
#include "VertexToolBase.h"
#include "lardataobj/MCBase/MCShower.h"

namespace signature 
{
    enum SignatureType {
        kEmptySignature = 0,
        kCosmicSignature,
        kNeutrinoSignature,
        kPrimaryMuonSignature,
        kChargedKaonSignature,
        kKaonShortSignature,
        kLambdaSignature,
        kChargedSigmaSignature
    };

    using Signature = std::vector<art::Ptr<simb::MCParticle>>;
    using Pattern = std::vector<std::pair<SignatureType, Signature>>;
    using ParticleMap = std::map<int, art::Ptr<simb::MCParticle>>;

    class SignatureToolBase {
    public:
        virtual ~SignatureToolBase() noexcept = default;

        virtual void configure(fhicl::ParameterSet const& pset) {
            _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
            _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
        }

        bool constructSignature(art::Event const& evt, Signature& signature) {
            signature.clear();
            auto const& truth_handle = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
            if (truth_handle->empty() || truth_handle->size() != 1) 
                return false;
            
            bool signature_found = false;
            this->findSignature(evt, signature, signature_found);
            if (!signature_found) 
                signature.clear();
            
            return signature_found;
        }

        virtual SignatureType getSignatureType() const {
            return kEmptySignature;
        }

        virtual bool isDetectable(art::Event const& evt, Signature const& signature) const {
            return !signature.empty();
        }

    protected:
        art::InputTag _MCPproducer;
        art::InputTag _MCTproducer;

        void fillSignature(const art::Ptr<simb::MCParticle>& mcp, Signature& signature) {
            signature.push_back(mcp);
        }

        virtual void findSignature(art::Event const& evt, Signature& signature, bool& signature_found) = 0;

        ParticleMap buildParticleMap(const art::Event& evt) const {
            auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
            ParticleMap mcp_map;
            
            for (size_t mcp_i = 0; mcp_i < mcp_h->size(); mcp_i++) {
                const art::Ptr<simb::MCParticle> mcp(mcp_h, mcp_i);
                mcp_map[mcp->TrackId()] = mcp;
            }
            
            return mcp_map;
        }

        template <typename Func>
        void traverseDecayChain(const art::Ptr<simb::MCParticle>& particle, 
                               const ParticleMap& mcp_map, 
                               Func callback) const {
            auto daughters = common::GetDaughters(particle, mcp_map);
            for (const auto& daughter : daughters) {
                callback(daughter);
                traverseDecayChain(daughter, mcp_map, callback);
            }
        }

        bool matchesDecayProducts(const art::Ptr<simb::MCParticle>& particle,
                                 const ParticleMap& mcp_map,
                                 const std::vector<int>& expected_pdgs) const {
            auto daughters = common::GetDaughters(particle, mcp_map);
            if (daughters.size() != expected_pdgs.size()) {
                return false;
            }

            std::vector<int> found_pdgs;
            for (const auto& daughter : daughters) {
                found_pdgs.push_back(daughter->PdgCode());
            }

            std::sort(found_pdgs.begin(), found_pdgs.end());
            
            std::vector<int> sorted_expected = expected_pdgs;
            std::sort(sorted_expected.begin(), sorted_expected.end());
            
            return found_pdgs == sorted_expected;
        }

        std::optional<TVector3> findDecayVertex(const art::Event& evt, 
                                              int parent_pdg, 
                                              const std::vector<int>& daughter_pdgs) const {
            auto mcp_map = buildParticleMap(evt);
            auto const& mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
            
            for (const auto& mcp : *mcp_h) {
                if (mcp.PdgCode() == parent_pdg && mcp.EndProcess() == "Decay") {
                    if (matchesDecayProducts(mcp_map.at(mcp.TrackId()), mcp_map, daughter_pdgs)) {
                        const TLorentzVector& end_position = mcp.EndPosition();
                        return TVector3(end_position.X(), end_position.Y(), end_position.Z());
                    }
                }
            }
            
            return std::nullopt;
        }
    };
} 

#endif