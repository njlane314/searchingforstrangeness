#ifndef SIGNATURE_INTEGRITY_H
#define SIGNATURE_INTEGRITY_H

#include "ClarityToolBase.h" 
#include "CommonFunctions/Corrections.h"

namespace claritytools {

class SignatureIntegrity : public ClarityToolBase {

public:
    explicit SignatureIntegrity(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)} 
    , _chan_act_reg{pset.get<int>("ChannelActiveRegion", 2)}
    , _max_bad_channel_frac{pset.get<double>("MaxBadChannelFraction", 0.3)}
    , _max_consecutive_bad_channel{pset.get<int>("MaxConsecutiveBadChannels",8)}
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
   const double _max_bad_channel_frac;
   const int _max_consecutive_bad_channel;

protected:

   bool isChannelRegionActive(const TVector3& point, const common::PandoraView& view) const; 
   bool checkStart(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const; 
   bool checkEnd(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const; 
   bool checkDeadChannelFrac(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const;

};

bool SignatureIntegrity::isChannelRegionActive(const TVector3& point, const common::PandoraView& view) const
{
  for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
    std::cout << "plane = " << plane.Plane << "  " << view << std::endl;
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

bool SignatureIntegrity::checkStart(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{
    /*
    TVector3 start(part->Vx(), part->Vy(), part->Vz());
    std::cout << "Start: " << start.X() << "  " << start.Y() << "  " << start.Z() << std::endl;
    common::ApplySCEMappingXYZ(start[0],start[1],start[2]);
    std::cout << "Corrected Start: " << start.X() << "  " << start.Y() << "  " << start.Z() << std::endl;
    return isChannelRegionActive(start,view);
    */
    float x = part->Vx();
    float y = part->Vy();
    float z = part->Vz();
    //std::cout << "Start: " << x << "  " << y << "  " << z << std::endl;
    common::ApplySCEMappingXYZ(x,y,z);
    //std::cout << "Corrected Start: " << x << "  " << y << "  " << z << std::endl;
    return isChannelRegionActive(TVector3(x,y,z),view);
}


bool SignatureIntegrity::checkEnd(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{
    /*
    TVector3 end(part->EndX(), part->EndY(), part->EndZ());
    std::cout << "End: " << end.X() << "  " << end.Y() << "  " << end.Z() << std::endl;
    return isChannelRegionActive(end,view);
    */
    float x = part->EndX();
    float y = part->EndY();
    float z = part->EndZ();
    common::ApplySCEMappingXYZ(x,y,z);    
    return isChannelRegionActive(TVector3(x,y,z),view);
}

// Check how many dead channels lie between the start and end of the track
bool SignatureIntegrity::checkDeadChannelFrac(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{

  std::cout << "checkDeadChannelFrac" << std::endl;

  float startx = part->Vx();
  float starty = part->Vy();
  float startz = part->Vz();
  common::ApplySCEMappingXYZ(startx,starty,startz);    
  float endx = part->EndX();
  float endy = part->EndY();
  float endz = part->EndZ();
  common::ApplySCEMappingXYZ(endx,endy,endz);    

  for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
    if(static_cast<unsigned int>(plane.Plane) != static_cast<unsigned int>(view)) continue;
    try {
      geo::WireID start_wire = _geo->NearestWireID(TVector3(startx,starty,startz), plane);
      raw::ChannelID_t start_channel = _geo->PlaneWireToChannel(start_wire);
      geo::WireID end_wire = _geo->NearestWireID(TVector3(endx,endy,endz), plane);
      raw::ChannelID_t end_channel = _geo->PlaneWireToChannel(end_wire);
      double channels = abs(static_cast<int>(start_channel) - static_cast<int>(end_channel));
      double bad_channels = 0;
      int cons_bad_ch = 0;
      int last_bad_ch = -1000;
      std::cout << "Channel range: " << std::min(start_channel,end_channel) << "  " << std::max(start_channel,end_channel) << std::endl;
      for(int ch=std::min(start_channel,end_channel);ch<=std::max(start_channel,end_channel);ch++){
        if(_bad_channel_mask[ch]){
          bad_channels++;
          std::cout << "Bad channel: " << ch << std::endl;
          if(last_bad_ch == ch - 1) cons_bad_ch++;
          else cons_bad_ch = 1; 
          last_bad_ch = ch;
          std::cout << "cons_bad_ch = " << cons_bad_ch << std::endl;
          if(cons_bad_ch > _max_consecutive_bad_channel) return false; 
        }
      }

      if(bad_channels/channels > _max_bad_channel_frac) return false; 

    }
    catch (const cet::exception&) {
      return false; 
    }
  }

  return true;

}

}

#endif
