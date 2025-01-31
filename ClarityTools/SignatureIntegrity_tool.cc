#ifndef SIGNATURE_INTEGRITY_H
#define SIGNATURE_INTEGRITY_H

#include "ClarityToolBase.h" 

namespace claritytools {

class SignatureIntegrity : ClarityToolBase {

public:
    explicit SignatureIntegrity(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)} 
    //, _chan_act_reg{pset.get<int>("ChannelActiveRegion", 3)}
    , _chan_act_reg{pset.get<int>("ChannelActiveRegion", 3)}
    {
        configure(pset);
    }

    ~SignatureIntegrity() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    //bool filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view);

private:

   const int _chan_act_reg;

};


bool SignatureIntegrity::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
//bool SignatureIntegrity::filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
{

  std::cout << "Checking Signature Integrity" << std::endl;
  this->loadEventHandles(e,view);

  //geo::PlaneID plane = view;

  auto isChannelRegionActive = [&](const TVector3& point) -> bool {
    for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
      //std::cout << "plane = " << plane.Plane << "  " << view << std::endl;
      //std::cout << "Casting: " << static_cast<unsigned int>(plane.Plane) << "  " << static_cast<unsigned int>(view) << std::endl;
      if(static_cast<unsigned int>(plane.Plane) != static_cast<unsigned int>(view)) continue;
      try {
        geo::WireID wire = _geo->NearestWireID(point, plane);
        raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

        for (int offset = -_chan_act_reg; offset <= _chan_act_reg; ++offset) {
          raw::ChannelID_t neighboring_channel = central_channel + offset;
          //std::cout << "Checking channel " << neighboring_channel << std::endl;
          if (neighboring_channel < 0 || static_cast<size_t>(neighboring_channel) >= _geo->Nchannels())
            continue; 

          if (_bad_channel_mask[neighboring_channel]){
            //std::cout << "Bad channel " << neighboring_channel << std::endl;
            return false; 
          }

        }
      } catch (const cet::exception&) {
        return false; 
      }
    }
    return true;
  };

  for (const auto& mcp_s : sig) {
    //std::cout << "Checking particle " << mcp_s->PdgCode() << "  " << mcp_s->TrackId() << std::endl;
   // std::cout << "Start " << mcp_s->Vx() << "  " <<  mcp_s->Vy() << "  " <<  mcp_s->Vz() << std::endl;
   // std::cout << "End " << mcp_s->EndX() << "  " <<  mcp_s->EndY() << "  " <<  mcp_s->EndZ() << std::endl;
    TVector3 start(mcp_s->Vx(), mcp_s->Vy(), mcp_s->Vz());
    TVector3 end(mcp_s->EndX(), mcp_s->EndY(), mcp_s->EndZ());

    if (!isChannelRegionActive(start))
      return false;

    if (std::abs(mcp_s->PdgCode()) != 13) { 
      if (!isChannelRegionActive(end))
        return false;
    }
  }

  return true;
}


DEFINE_ART_CLASS_TOOL(SignatureIntegrity)

}

#endif
