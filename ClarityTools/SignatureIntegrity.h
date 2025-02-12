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
    , _dist_to_scan{pset.get<double>("DistToScan",3)}
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
   const double _dist_to_scan;

protected:

   bool isChannelRegionActive(const TVector3& point, const common::PandoraView& view, int act_reg) const; 
   bool checkStart(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const; 
   bool checkStart2(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const; 
   bool checkEnd(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const; 
   bool checkEnd2(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const; 
   bool checkDeadChannelFrac(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const;

};

bool SignatureIntegrity::isChannelRegionActive(const TVector3& point, const common::PandoraView& view, int act_reg) const
{
  for (geo::PlaneID const& plane : _geo->IteratePlaneIDs()) {
    if(static_cast<unsigned int>(plane.Plane) != static_cast<unsigned int>(view)) continue;
    try {
      geo::WireID wire = _geo->NearestWireID(point, plane);
      raw::ChannelID_t central_channel = _geo->PlaneWireToChannel(wire);

      for (int offset = -act_reg; offset <= act_reg; ++offset) {
        raw::ChannelID_t neighboring_channel = central_channel + offset;
        if (neighboring_channel < 0 || static_cast<size_t>(neighboring_channel) >= _geo->Nchannels())
          continue; 

        if (_bad_channel_mask[neighboring_channel]){
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
    float x = part->Vx();
    float y = part->Vy();
    float z = part->Vz();
    common::ApplySCEMappingXYZ(x,y,z);
    bool pass = isChannelRegionActive(TVector3(x,y,z),view,_chan_act_reg);

    if(_verbose){
        if(pass) std::cout << "Track start ok" << std::endl;
        else std::cout << "Track start bad" << std::endl;
    }
    
   return pass;  
}

// Modified CheckStart method - check the first few points of the track, but only the exact channel they land on
// Intended to stop rejecting events with dead channels upstream of displaced vertex (which don't matter)
bool SignatureIntegrity::checkStart2(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{

    TVector3 start(part->Vx(),part->Vy(),part->Vz());

    for(int i_p=0;i_p<=part->NumberTrajectoryPoints();i_p++){

      float x = part->Vx(i_p);
      float y = part->Vy(i_p);
      float z = part->Vz(i_p);

      double dist = (TVector3(x,y,z) - start).Mag(); 
      if(dist > _dist_to_scan) break;

      common::ApplySCEMappingXYZ(x,y,z);
      bool pass = isChannelRegionActive(TVector3(x,y,z),view,0);

      if(_verbose && !pass){
        std::cout << "Track start bad" << std::endl;
        return false;
      }

    }
   
    std::cout << "Track start ok" << std::endl;
 
    return true;  
}

bool SignatureIntegrity::checkEnd(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{
    float x = part->EndX();
    float y = part->EndY();
    float z = part->EndZ();
    common::ApplySCEMappingXYZ(x,y,z);    
    bool pass = isChannelRegionActive(TVector3(x,y,z),view,_chan_act_reg);

    if(_verbose){
        if(pass) std::cout << "Track end ok" << std::endl;
        else std::cout << "Track end bad" << std::endl;
    }

    return pass;
}

bool SignatureIntegrity::checkEnd2(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{

    TVector3 end(part->EndX(),part->EndY(),part->EndZ());

    for(int i_p=part->NumberTrajectoryPoints();i_p>=0;i_p--){

      float x = part->Vx(i_p);
      float y = part->Vy(i_p);
      float z = part->Vz(i_p);

      double dist = (TVector3(x,y,z) - end).Mag(); 
      if(dist > _dist_to_scan) break;

      common::ApplySCEMappingXYZ(x,y,z);
      bool pass = isChannelRegionActive(TVector3(x,y,z),view,0);

      if(_verbose && !pass){
        std::cout << "Track end bad" << std::endl;
        return false;
      }

    }
   
    std::cout << "Track end ok" << std::endl;
 
    return true;  
}


// Check how many dead channels lie between the start and end of the track
bool SignatureIntegrity::checkDeadChannelFrac(const art::Ptr<simb::MCParticle>& part, const common::PandoraView view) const
{

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
      for(int ch=std::min(start_channel,end_channel);ch<=std::max(start_channel,end_channel);ch++){
        if(_bad_channel_mask[ch]){
          bad_channels++;
          if(last_bad_ch == ch - 1) cons_bad_ch++;
          else cons_bad_ch = 1; 
          last_bad_ch = ch;
          if(cons_bad_ch > _max_consecutive_bad_channel){
            if(_verbose)
              std::cout << "Tool many consecutive dead channels" << std::endl;
            return false; 
          }
        }
      }

      if(_verbose){
        std::cout << "Bad channels = " << bad_channels << std::endl;
        std::cout << "Bad channel fraction = " << bad_channels/channels << std::endl;
      }

      if(bad_channels/channels > _max_bad_channel_frac){
        if(_verbose) 
          std::cout << "Failed bad channel fraction" << std::endl;
        return false; 
      }
    }
    catch (const cet::exception&) {
      return false; 
    }
  }

  return true;

}

}

#endif
