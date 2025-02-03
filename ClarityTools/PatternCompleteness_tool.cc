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
      , _part_hit_thresh{pset.get<int>("ParticleHitThreshold", 5)}
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
    std::cout << "Checking Pattern Completeness" << std::endl; 
    this->loadEventHandles(e,view);

    std::unordered_map<int, int> sig_hit_map;
    double tot_sig_hit = 0; 

    std::vector<art::Ptr<recob::Hit>> sig_hits;
    for (const auto& mcp_s : sig) {
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
    }

    if (_mc_hits.empty() || sig_hits.empty()) 
        return false;

    //double sig_comp = static_cast<double>(sig_hits.size()) / _mc_hits.size();

    //std::cout << "Signature completeness " << sig_comp << std::endl;
    std::cout << "Total signature hits " << tot_sig_hit << std::endl;

    for (const auto& [trackid, num_hits] : sig_hit_map) 
    {
        std::cout << "Particle " << trackid << " hits " << num_hits << std::endl;
        if(num_hits < _part_hit_thresh) return false;
        if (num_hits / tot_sig_hit < _part_hit_frac_thresh) return false;       
    }

    if (/*sig_comp < _sig_hit_comp_thresh ||*/ tot_sig_hit < _sig_hit_thresh)
        return false;

    return true;
}

DEFINE_ART_CLASS_TOOL(PatternCompleteness)

}

#endif
