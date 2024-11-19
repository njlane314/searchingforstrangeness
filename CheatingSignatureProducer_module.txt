#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/Utilities/AssociationUtil.h"  

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h" 

#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/DecayVertexProvider.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>

class CheatingSignatureProducer : public art::EDProducer {
public:
    explicit CheatingSignatureProducer(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt) override;
    void beginJob() override;
    void endJob() override;

private:
    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag, _PFPproducer;
    art::InputTag _OutputTag; 

    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;
};

CheatingSignatureProducer::CheatingSignatureProducer(fhicl::ParameterSet const& pset)
    : EDProducer{pset},
      _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")},
      _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")},
      _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")},
      _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")},
      _PFPproducer{pset.get<art::InputTag>("PFPproducer", "pandora")}, 
      _OutputTag{pset.get<art::InputTag>("OutputTag", "pandoracheating")}
{
    produces<std::vector<recob::Vertex>>(_OutputTag.instance());
    produces<std::vector<recob::Hit>>(_OutputTag.instance());
    produces<std::vector<recob::Slice>>(_OutputTag.instance());
    produces<art::Assns<recob::Slice, recob::Vertex>>(_OutputTag.instance());
    produces<art::Assns<recob::Slice, recob::Hit>>(_OutputTag.instance());

    const auto& tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (const auto& tool_pset_label : tool_psets.get_pset_names()) {
        const auto tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }
}

void CheatingSignatureProducer::produce(art::Event& evt) 
{
    std::vector<signature::Signature> signature_coll;
    bool found_all_signatures = true;
    for (auto& signatureTool : _signatureToolsVec) {
        bool found_signature = signatureTool->identifySignalParticles(evt, signature_coll);
        if (!found_signature) 
            found_all_signatures = false;
    }

    auto outputVertex = std::make_unique<std::vector<recob::Vertex>>();
    auto outputHitVector = std::make_unique<std::vector<recob::Hit>>();
    auto outputSlice = std::make_unique<std::vector<recob::Slice>>();
    auto outputPFParticleVector = std::make_unique<std::vector<int>>();
    auto outputSliceVertexAssoc = std::make_unique<art::Assns<recob::Slice, recob::Vertex>>();
    auto outputSliceHitAssoc = std::make_unique<art::Assns<recob::Slice, recob::Hit>>();

    if (!found_all_signatures) {
        evt.put(std::move(outputVertex));
        evt.put(std::move(outputHitVector));
        evt.put(std::move(outputSlice));
        evt.put(std::move(outputPFParticleVector));
        evt.put(std::move(outputSliceVertexAssoc));
        evt.put(std::move(outputSliceHitAssoc));
        return;
    }

    std::vector<art::Ptr<recob::Hit>> all_hits;
    art::Handle<std::vector<recob::Hit>> hit_handle;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> mcp_bkth_assoc;
    if (evt.getByLabel(_HitProducer, hit_handle)) {
        art::fill_ptr_vector(all_hits, hit_handle);
        mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, evt, _BacktrackTag);
    }

    std::vector<art::Ptr<recob::Hit>> signature_hits;
    for (const auto& hit : all_hits) {
        bool is_associated_with_signature = false;

        if (mcp_bkth_assoc) {
            const auto& associated_particles = mcp_bkth_assoc->at(hit.key());
            const auto& match_data = mcp_bkth_assoc->data(hit.key());

            for (size_t i = 0; i < associated_particles.size(); ++i) {
                if (match_data[i]->isMaxIDE == 1) {
                    int track_id = associated_particles[i]->TrackId();
                    for (const auto& signature : signature_coll) {
                        if (track_id == signature.trckid) {
                            is_associated_with_signature = true;
                            break;
                        }
                    }
                }
                if (is_associated_with_signature)
                    break;
            }
        }

        if (is_associated_with_signature) 
            signature_hits.push_back(hit);
    }

    for (const auto& hit : signature_hits) 
        outputHitVector->push_back(*hit);

    int sliceID = 0;  
    recob::Slice::Point_t sliceCenter(0.0, 0.0, 0.0);  
    recob::Slice::Vector_t sliceDirection(0.0, 0.0, 0.0);  
    recob::Slice::Point_t end0Pos(0.0, 0.0, 0.0);  
    recob::Slice::Point_t end1Pos(0.0, 0.0, 0.0);  
    float aspectRatio = 0.0;  
    float totalCharge = 0.0; 
    recob::Slice newSlice(sliceID, sliceCenter, sliceDirection, end0Pos, end1Pos, aspectRatio, totalCharge);
    outputSlice->push_back(newSlice);

    for (const auto& signatureTool : _signatureToolsVec) {
        auto* decayTool = dynamic_cast<DecayVertexProvider*>(signatureTool.get());
        if (decayTool) {
            std::optional<TVector3> decay_vertex_opt = decayTool->getDecayVertex(evt);
            if (decay_vertex_opt) {
                TVector3 decay_vertex = *decay_vertex_opt;
                double xyz[3] = {decay_vertex.X(), decay_vertex.Y(), decay_vertex.Z()};
                recob::Vertex decayVertex(xyz, -1);
                outputVertex->push_back(decayVertex);
                util::CreateAssn(evt, *outputSlice, art::Ptr<recob::Vertex>(outputVertex->back()), *outputSliceVertexAssoc);
                break;
            }
        }
    }

    for (const auto& hit : signature_hits) {
        util::CreateAssn(evt, *outputSlice, hit, *outputSliceHitAssoc);
    }

    evt.put(std::move(outputVertex), _OutputTag.instance());
    evt.put(std::move(outputSlice), _OutputTag.instance());
    evt.put(std::move(outputSliceVertexAssoc), _OutputTag.instance());
    evt.put(std::move(outputSliceHitAssoc), _OutputTag.instance());
}

void CheatingSignatureProducer::beginJob() {}

void CheatingSignatureProducer::endJob() {}

DEFINE_ART_MODULE(CheatingSignatureProducer)