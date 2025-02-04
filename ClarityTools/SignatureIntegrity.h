#ifndef SIGNATURE_INTEGRITY_H
#define SIGNATURE_INTEGRITY_H

#include "ClarityToolBase.h" 

namespace claritytools {

class SignatureIntegrity : public ClarityToolBase {

public:
    explicit SignatureIntegrity(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)} 
    , _chan_act_reg{pset.get<int>("ChannelActiveRegion", 2)}
    {
        configure(pset);
    }

    ~SignatureIntegrity() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    virtual bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view) = 0;

private:

   const int _chan_act_reg;

protected:
   bool isChannelRegionActive(const TVector3& point, const common::PandoraView& view) const; 
   bool checkStart(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const; 
   bool checkEnd(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const; 

};

bool SignatureIntegrity::isChannelRegionActive(const TVector3& point, const common::PandoraView& view) const
{
  for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
    //std::cout << "plane = " << plane.Plane << "  " << view << std::endl;
    //std::cout << "Casting: " << static_cast<unsigned int>(plane.Plane) << "  " << static_cast<unsigned int>(view) << std::endl;
    if(static_cast<unsigned int>(plane.Plane) != static_cast<unsigned int>(view)) continue;
    try {
      geo::WireID wire = _geo->NearestWireID(point, plane);
      raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

      for (int offset = -_chan_act_reg; offset <= _chan_act_reg; ++offset) {
        raw::ChannelID_t neighboring_channel = central_channel + offset;
        std::cout << "Checking channel " << neighboring_channel << std::endl;
        if (neighboring_channel < 0 || static_cast<size_t>(neighboring_channel) >= _geo->Nchannels())
          continue; 

        if (_bad_channel_mask[neighboring_channel]){
          std::cout << "Bad channel " << neighboring_channel << std::endl;
          return false; 
        }

      }
    } catch (const cet::exception&) {
      return false; 
    }
  }
  return true;
}
/*
bool SignatureIntegrity::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{

  std::cout << "Checking General Signature Integrity" << std::endl;
  this->loadEventHandles(e,view);

  for (const auto& mcp_s : sig.second) {
    if(!checkStart(mcp_s,view) || !checkEnd(mcp_s,view)) return false;
  }

  return true;
}
*/
bool SignatureIntegrity::checkStart(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const
{
    TVector3 start(part->Vx(), part->Vy(), part->Vz());
    std::cout << "Start: " << start.X() << "  " << start.Y() << "  " << start.Z() << std::endl;
    return isChannelRegionActive(start,view);
}


bool SignatureIntegrity::checkEnd(const art::Ptr<simb::MCParticle>& part, common::PandoraView view) const
{
    TVector3 end(part->EndX(), part->EndY(), part->EndZ());
    std::cout << "End: " << end.X() << "  " << end.Y() << "  " << end.Z() << std::endl;
    return isChannelRegionActive(end,view);
}

}

#endif
