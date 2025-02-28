#ifndef PATTERN_COMPLETENESS_H
#define PATTERN_COMPLETENESS_H

#include "ClarityToolBase.h" 
#include "CommonFunctions/Geometry.h"

namespace claritytools {

std::vector<TLorentzVector> PadOutTrajPoints(const art::Ptr<simb::MCParticle>& part, double cutoff=1000){
 
  std::vector<TLorentzVector> points;
  double d=0; 

  for(int i_p=1;i_p<part->NumberTrajectoryPoints();i_p++){
    TLorentzVector last = part->Position(i_p-1);
    TLorentzVector current = part->Position(i_p);
    d += (current.Vect() - last.Vect()).Mag();
    TLorentzVector move = 0.2*(current - last);
    for(int i=0;i<5;i++)
      points.push_back(last+i*move);
    if(d > cutoff) break;
  }

  return points;
 
}

class PatternCompleteness : ClarityToolBase {

public:
    explicit PatternCompleteness(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)} 
      , _dist_to_scan{pset.get<double>("DistToScan",15)}
      , _part_comp_thresh{pset.get<double>("ParticleCompletenessThreshold",0.7)}
      , _max_dist{pset.get<double>("MaxDistance",0.6)}  
    {
        configure(pset);
    }

    ~PatternCompleteness() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    bool filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view);
    

private:

   const double _dist_to_scan;
   const double _part_comp_thresh;
   const double _max_dist;

   const double _at = 1.0/18.2148;
   const double _aw = 1.0/3.33328;

};

bool PatternCompleteness::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{

  if(_verbose)
    std::cout << "Checking PatternCompleteness in view " << view << " for signature " << signature::GetSignatureName(sig) << std::endl;

  this->loadEventHandles(e,view);

  for(const art::Ptr<simb::MCParticle>& part : sig.second){

    if(_verbose)
      std::cout << "Checking particle pdg=" << part->PdgCode() << " trackid=" << part->TrackId() << std::endl;


    std::vector<TLorentzVector> points = PadOutTrajPoints(part,_dist_to_scan);

    TVector3 start(part->Vx(),part->Vy(),part->Vz());
 
    int points_scanned = 0;
    int good_points = 0;

    for(TLorentzVector pos : points){

      if(_verbose)
        std::cout << "Checking point at dist: " << (pos.Vect() - start).Mag() << std::endl;

      double a_pos[3] = {pos.X(),pos.Y(),pos.Z()};
      if(!common::point_inside_fv(a_pos)){
        if(_verbose) 
          std::cout << "Track goes out of FV" << std::endl;
        return false;
      }

      common::True2RecoMapping(pos);

      points_scanned++;

      // Check if there is a hit near this point that truth matches to the particle
      double nearest_hit_dist2 = 10000; 
      art::Ptr<recob::Hit> nearest_hit;
      for (art::Ptr<recob::Hit> hit : _mc_hits) {
        if(hit->View() != view) continue;
        float d2 = common::HitPtDistance(pos.Vect(),hit,_aw,_at);
        if(d2 < nearest_hit_dist2){
          nearest_hit_dist2 = d2;
          nearest_hit = hit; 
        }
      }

      if(nearest_hit_dist2 == 10000) continue;

      if(_verbose)
        std::cout << "Nearest hit separation = " << sqrt(nearest_hit_dist2) << std::endl;

      if(sqrt(nearest_hit_dist2) > _max_dist) continue;
        
      auto assmcp = _mcp_bkth_assoc->at(nearest_hit.key());
      auto assmdt = _mcp_bkth_assoc->data(nearest_hit.key());

      assert(assmcp.size() == assmdt.size());

      bool found_hit = false;
      for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
        auto amd = assmdt[ia];
        if (assmcp[ia]->TrackId() == part->TrackId() && amd->isMaxIDEN == 1) {
          found_hit = true;
          break;
        }
      }
       
      if(!found_hit && _verbose){
        std::cout << "Can't find any nearby hit " << (pos.Vect() - start).Mag() << "cm along track" << std::endl;
      }  

      if(found_hit)
        good_points++;
        
    }
    

    if(points_scanned == 0) return false;

    if(_verbose)
      std::cout << "Good points/points scanned = " << good_points << "/" << points_scanned << "=" << (double)good_points/points_scanned << std::endl;

    if((double)good_points/points_scanned < _part_comp_thresh){
      if(_verbose)
        std::cout << "Track failed" << std::endl; 
      return false;
    }
       
  }

  return true;

}

DEFINE_ART_CLASS_TOOL(PatternCompleteness)

}

#endif
