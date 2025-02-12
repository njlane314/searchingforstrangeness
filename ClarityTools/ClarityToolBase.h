#ifndef CLARITY_TOOLBASE_H
#define CLARITY_TOOLBASE_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"
#include "CommonFunctions/BadChannels.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/VertexToolBase.h"

#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "TDatabasePDG.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <chrono>

namespace claritytools {

class ClarityToolBase {

public:

    ClarityToolBase(fhicl::ParameterSet const& pset) :
      _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _DeadChannelTag{pset.get<art::InputTag>("DeadChannelTag")}
    , _verbose{pset.get<bool>("Verbose",false)}
    {
    }   
 
    virtual ~ClarityToolBase() noexcept = default;
    
    virtual void configure(fhicl::ParameterSet const& pset)
    {
    }

    std::map<common::PandoraView,std::vector<bool>> filter3Plane(const art::Event &e, const signature::Pattern& patt);
    std::vector<bool> filter(const art::Event &e, const signature::Pattern& patt, common::PandoraView view);
    virtual bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view) = 0;

private:

    const art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag, _DeadChannelTag;

protected:

    std::vector<bool> _bad_channel_mask;

    const geo::GeometryCore* _geo = art::ServiceHandle<geo::Geometry>()->provider();
     
    bool loadEventHandles(const art::Event &e, common::PandoraView targetDetectorPlane);
    std::vector<art::Ptr<recob::Hit>> _evt_hits;
    std::vector<art::Ptr<recob::Hit>> _mc_hits;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

    const bool _verbose;

};

bool ClarityToolBase::loadEventHandles(const art::Event &e, common::PandoraView targetDetectorPlane){

    common::SetBadChannelMask(e,_DeadChannelTag,_bad_channel_mask);

    _evt_hits.clear();
    _mc_hits.clear();

    art::Handle<std::vector<recob::Hit>> hit_h;
    if (!e.getByLabel(_HitProducer, hit_h)) 
        return false;

    art::fill_ptr_vector(_evt_hits, hit_h);
    _mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_h, e, _BacktrackTag);

    for (const auto& hit : _evt_hits) {
        if (_bad_channel_mask[hit->Channel()]) 
            continue; 

        const geo::WireID& wire_id = hit->WireID(); 
        if (wire_id.Plane != static_cast<unsigned int>(targetDetectorPlane))
            continue;

        auto assmcp = _mcp_bkth_assoc->at(hit.key());
        auto assmdt = _mcp_bkth_assoc->data(hit.key());
        for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
            auto amd = assmdt[ia];
            if (amd->isMaxIDEN != 1)
                continue;
            
            _mc_hits.push_back(hit);
        }
    }

    return true;

}

std::map<common::PandoraView,std::vector<bool>> ClarityToolBase::filter3Plane(const art::Event &e, const signature::Pattern& patt){

  std::map<common::PandoraView,std::vector<bool>> result;

  for(int view = common::TPC_VIEW_U;view != common::N_VIEWS; view++){ 
    result[static_cast<common::PandoraView>(view)] = this->filter(e,patt,static_cast<common::PandoraView>(view));
  }

  return result;

}

std::vector<bool> ClarityToolBase::filter(const art::Event &e, const signature::Pattern& patt, common::PandoraView view){

  std::vector<bool> result;
  for (const auto& sig : patt) {
    result.push_back(this->filter(e,sig,view));
  }

  return result;

}

}

#endif
