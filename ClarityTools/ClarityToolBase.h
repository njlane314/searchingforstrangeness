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
    , _DeadChannelTag{pset.get<art::InputTag>("DeadChannelTag")}
    {
    }   
 
    virtual ~ClarityToolBase() noexcept = default;
    
    virtual void configure(fhicl::ParameterSet const& pset)
    {
        //_bad_channel_file = pset.get<std::string>("BadChannelFile");
        //loadBadChannelMap();
    }

    std::map<common::PandoraView,std::vector<bool>> filter3Plane(const art::Event &e, const signature::Pattern& patt);
    std::vector<bool> filter(const art::Event &e, const signature::Pattern& patt, common::PandoraView view);
    virtual bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view) = 0;

protected:

    std::vector<bool> _bad_channel_mask;
    //std::string _bad_channel_file;

    const geo::GeometryCore* _geo = art::ServiceHandle<geo::Geometry>()->provider();
     
    bool loadEventHandles(const art::Event &e, common::PandoraView targetDetectorPlane);
    std::vector<art::Ptr<recob::Hit>> _evt_hits;
    std::vector<art::Ptr<recob::Hit>> _mc_hits;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

private:

    const art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag, _DeadChannelTag;

    //void loadBadChannelMap();

};
/*
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

  size_t num_channels = _geo->Nchannels();
  _bad_channel_mask.resize(num_channels+1,false);

  std::vector<int> bad_ch;
  std::ifstream inFile(_bad_channel_file, std::ios::in);
  std::string line;
  int ch;
  while (inFile >> ch) {
     bad_ch.push_back(ch);
  } 

  for(size_t i=0;i<num_channels+1;i++)
    if(std::find(bad_ch.begin(),bad_ch.end(),i) != bad_ch.end())
      _bad_channel_mask.at(i) = true; 

}
*/
bool ClarityToolBase::loadEventHandles(const art::Event &e, common::PandoraView targetDetectorPlane){

    size_t num_channels = _geo->Nchannels();
    _bad_channel_mask.resize(num_channels+1,false);

    art::Handle<std::vector<int>> bad_ch_h;
    std::vector<art::Ptr<int>> bad_ch_v;

    if(!e.getByLabel(_DeadChannelTag,bad_ch_h))
        return false;

    art::fill_ptr_vector(bad_ch_v,bad_ch_h);
   
    for(auto ch : bad_ch_v){
      //std::cout << *ch << std::endl;
      _bad_channel_mask.at(*ch) = true; 
    }


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
