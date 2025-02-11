#ifndef _BAD_CHANNELS_H_
#define _BAD_CHANNELS_H_

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

namespace common {

void SetBadChannelMask(const art::Event& e, const art::InputTag& bad_channel_tag, std::vector<bool>& bad_channel_mask)
{

  //std::cout << "bad_channel_tag = " << bad_channel_tag << std::endl;

  const geo::GeometryCore* geo = art::ServiceHandle<geo::Geometry>()->provider();

  size_t num_channels = geo->Nchannels();
  //std::cout << "num_channels = " <<  num_channels << std::endl;
  bad_channel_mask.resize(num_channels, false);

  art::Handle<std::vector<int>> bad_ch_h;
  std::vector<art::Ptr<int>> bad_ch_v;

  if(!e.getByLabel(bad_channel_tag,bad_ch_h))
    throw cet::exception("SetBadChannelMask:") << "Bad channel product with tag" << bad_channel_tag << " not found" << std::endl;
    

  art::fill_ptr_vector(bad_ch_v,bad_ch_h);

  for(auto ch : bad_ch_v){
    //std::cout << *ch << std::endl;
    bad_channel_mask.at(*ch) = true; 
  }

}

}
#endif
