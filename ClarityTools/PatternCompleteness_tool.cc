#ifndef PATTERN_COMPLETENESS_H
#define PATTERN_COMPLETENESS_H

#include "ClarityToolBase.h" 

namespace claritytools {

class PatternCompleteness : ClarityToolBase {

public:
    explicit PatternCompleteness(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)} 
      , _sig_hit_thresh{pset.get<int>("SignatureHitThreshold", 10)}
      , _sig_hit_comp_thresh{pset.get<double>("SignatureHitCompletenessThreshold", 0.05)}
      , _part_hit_thresh{pset.get<int>("ParticleHitThreshold", 4)}
      , _part_hit_frac_thresh{pset.get<double>("ParticleHitFractionThreshold", 0.05)}
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

   const int _sig_hit_thresh;
   const double _sig_hit_comp_thresh;
   const int _part_hit_thresh;
   const double _part_hit_frac_thresh;

};

bool PatternCompleteness::filter(const art::Event &e, const signature::Signature& sig, common::PandoraView view)
{

  if(_verbose)
    std::cout << "Checking PatternCompleteness in view " << view << " for signature " << signature::GetSignatureName(sig) << std::endl;

    this->loadEventHandles(e,view);

    std::unordered_map<int, int> sig_hit_map;
    double tot_sig_hit = 0; 

    std::vector<art::Ptr<recob::Hit>> sig_hits;
    for (const auto& mcp_s : sig.second) {
      double sig_hit = 0;

      for (const auto& hit : _mc_hits) {
        auto assmcp = _mcp_bkth_assoc->at(hit.key());
        auto assmdt = _mcp_bkth_assoc->data(hit.key());

        for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
          auto amd = assmdt[ia];
          if (assmcp[ia]->TrackId() == mcp_s->TrackId() && amd->isMaxIDEN == 1) {
            sig_hits.push_back(hit);
            sig_hit += 1; 
          }
        }
      }

      sig_hit_map[mcp_s->TrackId()] += sig_hit;
      tot_sig_hit += sig_hit;
      if(_verbose)
        std::cout << "Particle pdg=" << mcp_s->PdgCode() << " trackid=" << mcp_s->TrackId() << " hits = " << sig_hit << std::endl;
    }

    if (_mc_hits.empty() || sig_hits.empty()) 
        return false;

    for (const auto& [trackid, num_hits] : sig_hit_map) 
    {
      if(_verbose)
        std::cout << "Signature trackid=" << trackid << " num_hits = " << num_hits << " num_hits/tot_sig_hit = " << num_hits/tot_sig_hit << std::endl;
      if(num_hits < _part_hit_thresh){
        if(_verbose)
          std::cout << "Signature trackid=" << trackid << " failed PatternCompleteness in plane " << view << " with too few hits" << std::endl;
        return false;
      }
      if (num_hits / tot_sig_hit < _part_hit_frac_thresh){
        if(_verbose)
          std::cout << "Signature trackid=" << trackid << " failed PatternCompleteness in plane " << view << " with too small frac of hits" << std::endl;
        return false;       
      }
    }

    if (tot_sig_hit < _sig_hit_thresh){
      if(_verbose)
        std::cout << "Signature " << signature::GetSignatureName(sig) << " failed PatternCompleteness with too few hits" << std::endl;
      return false;
    }

    if(_verbose)
      std::cout << "Signature " << signature::GetSignatureName(sig) << " passed PatternCompleteness in plane " << view << std::endl;
    
    return true;
}

DEFINE_ART_CLASS_TOOL(PatternCompleteness)

}

#endif
