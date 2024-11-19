#include "art/Framework/Core/EDFilter.h"
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

#include "CommonFunctions/Region.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>

class PatternRecognitionFilter : public art::EDFilter
{
public:
    explicit PatternRecognitionFilter(fhicl::ParameterSet const &p);

    PatternRecognitionFilter(PatternRecognitionFilter const &) = delete;
    PatternRecognitionFilter(PatternRecognitionFilter &&) = delete;
    PatternRecognitionFilter &operator=(PatternRecognitionFilter const &) = delete;
    PatternRecognitionFilter &operator=(PatternRecognitionFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _MCPproducer, _HitProducer, _BacktrackTag, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer;
    
    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;
};

PatternRecognitionFilter::PatternRecognitionFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _MCPproducer(pset.get<art::InputTag>("MCPproducer", "largeant"))
    , _HitProducer(pset.get<art::InputTag>("Hproducer", "gaushit"))
    , _BacktrackTag(pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch"))
    , _PFPproducer{pset.get<art::InputTag>("PFPproducer", "pandora")}
    , _CLSproducer{pset.get<art::InputTag>("CLSproducer", "pandora")}
    , _SHRproducer{pset.get<art::InputTag>("SHRproducer", "pandora")}
    , _SLCproducer{pset.get<art::InputTag>("SLCproducer", "pandora")}
    , _VTXproducer{pset.get<art::InputTag>("VTXproducer", "pandora")}
    , _PCAproducer{pset.get<art::InputTag>("PCAproducer", "pandora")}
    , _TRKproducer{pset.get<art::InputTag>("TRKproducer", "pandora")}
{
    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }
}

bool PatternRecognitionFilter::filter(art::Event &evt)
{
    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(evt, _PFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(_PFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(_CLSproducer),
                                                        proxy::withAssociated<recob::Slice>(_SLCproducer),
                                                        proxy::withAssociated<recob::Track>(_TRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(_VTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(_PCAproducer),
                                                        proxy::withAssociated<recob::Shower>(_SHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(_PFPproducer));

    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(evt, _CLSproducer,
                                                proxy::withAssociated<recob::Hit>(_CLSproducer));

    std::vector<signature::Signature> sig_coll;
    for (auto& signatureTool : _signatureToolsVec) {
        if (!signatureTool->identifySignalParticles(evt, sig_coll))
            return false;
    }

    auto const &all_hits = evt.getValidHandle<std::vector<recob::Hit>>(_HitProducer);
    auto mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(all_hits, evt, _BacktrackTag);

    std::unordered_map<int, int> sig_mcp_hits;  
    for (unsigned int ih = 0; ih < all_hits->size(); ih++)
    {
        auto assmcp = mcp_bkth_assoc->at(ih);
        auto assmdt = mcp_bkth_assoc->data(ih);
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
        {
            auto mcp = assmcp[ia];
            auto amd = assmdt[ia];
            if (amd->isMaxIDE != 1)
                continue;

            for (const auto &sig : sig_coll) 
            {
                if (mcp->TrackId() == sig.trckid)
                    sig_mcp_hits[sig.trckid]++;
            }
        }
    }

    std::unordered_map<int, int> sig_pfp_hits;
    std::unordered_map<int, int> pfp_total_hits;
    std::unordered_map<int, std::unordered_map<int, int>> pfp_mcp_shared_hits;
    auto [_, nu_slice] = common::getNuSliceHits(pfp_proxy, clus_proxy);
    for (const common::ProxyPfpElem_t &pfp_pxy : nu_slice)
    {
        if (pfp_pxy->IsPrimary())
            continue;

        std::vector<art::Ptr<recob::Hit>> pfp_hits;
        auto clus_pxy_v = pfp_pxy.get<recob::Cluster>();

        for (auto ass_clus : clus_pxy_v)
        {
            const auto &clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();

            for (const auto &hit : clus_hit_v)
                pfp_hits.push_back(hit);
        } 

        pfp_total_hits[pfp_pxy->Self()] = pfp_hits.size();

        for (auto hit : pfp_hits)
        {
            auto assmcp = mcp_bkth_assoc->at(hit.key());       
            auto assmdt = mcp_bkth_assoc->data(hit.key());    
            for (size_t i = 0; i < assmcp.size(); i++)
            {
                if (assmdt[i]->isMaxIDE != 1) 
                    continue;

                for (const auto &sig : sig_coll)
                {
                    if (assmcp[i]->TrackId() == sig.trckid)
                        pfp_mcp_shared_hits[pfp_pxy->Self()][sig.trckid]++;
                }
            }
        }

        std::unordered_map<int, int> sig_match;
        for (const auto &sig : sig_coll) 
        {
            int match_pfp = -1;
            float match_score = -1.0;

            for (const auto &[pfp_id, mcp_hits_map] : pfp_mcp_shared_hits) 
            {
                if (mcp_hits_map.find(sig.trckid) == mcp_hits_map.end())
                    continue;

                int shared_hits = mcp_hits_map.at(sig.trckid);
                float purity = (float)shared_hits / pfp_total_hits[pfp_id];
                float completeness = (float)shared_hits / sig_mcp_hits[sig.trckid];
                if (purity <= 0.5 || completeness <= 0.1)
                    continue;

                float score = std::sqrt(purity * purity + completeness * completeness);

                if (score > match_score) 
                {
                    match_score = score;
                    match_pfp = pfp_id;
                }
            }

            if (match_pfp != -1) 
                sig_match[sig.trckid] = match_pfp;
        }

        std::unordered_set<int> unique_pfps;
        for (const auto &[_, pfp_id] : sig_match) 
        {
            if (unique_pfps.find(pfp_id) != unique_pfps.end())
                return false; 

            unique_pfps.insert(pfp_id);
        }

        if (sig_match.size() != sig_coll.size())
            return false;
    }

    return true;
}

DEFINE_ART_MODULE(PatternRecognitionFilter)