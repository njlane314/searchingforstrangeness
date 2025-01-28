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
    {
    }   
 
    virtual ~ClarityToolBase() noexcept = default;
    
    virtual void configure(fhicl::ParameterSet const& pset)
    {
        _bad_channel_file = pset.get<std::string>("BadChannelFile");
        loadBadChannelMap();
    }

    //virtual bool filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc) = 0;
    std::map<common::PandoraView,bool> filter3Plane(const art::Event &e, const signature::Pattern& patt);
    virtual bool filter(const art::Event &e, const signature::Pattern& patt, common::PandoraView view) = 0;

protected:

    std::vector<bool> _bad_channel_mask;
    std::string _bad_channel_file;

    const geo::GeometryCore* _geo = art::ServiceHandle<geo::Geometry>()->provider();
     
    bool loadEventHandles(const art::Event &e, common::PandoraView targetDetectorPlane);
    std::vector<art::Ptr<recob::Hit>> _evt_hits;
    std::vector<art::Ptr<recob::Hit>> _mc_hits;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

private:

    const art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag;
    void loadBadChannelMap();

};

void ClarityToolBase::loadBadChannelMap(){

    std::cout << "Loading Bad Channel Map" << std::endl;

    size_t num_channels = _geo->Nchannels();
    _bad_channel_mask.resize(num_channels, false);

    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) 
            throw cet::exception("PatternClarityFilter") << "-- Bad channel file not found: " << _bad_channel_file;

        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        while (std::getline(inFile, line)) {
            if (line.find("#") != std::string::npos) continue;
            std::istringstream ss(line);
            int ch1, ch2;
            ss >> ch1;
            if (!(ss >> ch2)) ch2 = ch1;
            for (int i = ch1; i <= ch2; ++i) {
                _bad_channel_mask[i] = true;
            }
        }
    }

}

bool ClarityToolBase::loadEventHandles(const art::Event &e, common::PandoraView targetDetectorPlane){

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

std::map<common::PandoraView,bool> ClarityToolBase::filter3Plane(const art::Event &e, const signature::Pattern& patt){


  std::map<common::PandoraView,bool> result;

  for(int view = common::TPC_VIEW_U;view != common::N_VIEWS; view++){ 

    result[static_cast<common::PandoraView>(view)] = this->filter(e,patt,static_cast<common::PandoraView>(view));
  }

  return result;

}

}

#endif
