#ifndef PATTERN_COMPLETENESS_H
#define PATTERN_COMPLETENESS_H

#include "ClarityToolBase.h" 

namespace claritytools {

class PatternCompleteness : ClarityToolBase {

public:
    explicit PatternCompleteness(const fhicl::ParameterSet& pset) :
      ClarityToolBase{(pset)} 
      , _patt_hit_comp_thresh{pset.get<double>("PatternHitCompletenessThreshold", 0.5)}
      , _patt_hit_thresh{pset.get<int>("PatternHitThreshold", 100)}
      , _sig_hit_comp_thresh{pset.get<double>("SignatureHitCompletenessThreshold", 0.1)}
    {
        configure(pset);
    }

    ~PatternCompleteness() override = default;
    
    void configure(fhicl::ParameterSet const& pset) override
    {
        ClarityToolBase::configure(pset);
    }

    //bool filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filter(const art::Event &e, const signature::Pattern& patt, common::PandoraView view);

private:

   const double _patt_hit_comp_thresh;
   const int _patt_hit_thresh;
   const double _sig_hit_comp_thresh;

};

//bool PatternCompleteness::filter(const art::Event &e, const signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc)
bool PatternCompleteness::filter(const art::Event &e, const signature::Pattern& patt, common::PandoraView view)
{

    //std::cout << "Testing PatternCompleteness" << std::endl;
    this->loadEventHandles(e,view);

    std::unordered_map<int, int> sig_hit_map;
    double tot_patt_hit = 0; 

    std::vector<art::Ptr<recob::Hit>> patt_hits;
    for (const auto& sig : patt) {
        for (const auto& mcp_s : sig) {
            double sig_hit = 0;

            for (const auto& hit : _mc_hits) {
                auto assmcp = _mcp_bkth_assoc->at(hit.key());
                auto assmdt = _mcp_bkth_assoc->data(hit.key());

                for (unsigned int ia = 0; ia < assmcp.size(); ++ia){
                    auto amd = assmdt[ia];
                    if (assmcp[ia]->TrackId() == mcp_s->TrackId() && amd->isMaxIDEN == 1) {
                        patt_hits.push_back(hit);
                        sig_hit += 1; 
                    }
                }
            }

            sig_hit_map[mcp_s->TrackId()] += sig_hit;
            tot_patt_hit += sig_hit;
        }
    }

    if (_mc_hits.empty() || patt_hits.empty()) 
        return false;

    double patt_comp = static_cast<double>(patt_hits.size()) / _mc_hits.size();
    //std::cout << "Pattern completeness " << patt_comp << std::endl;
    //std::cout << "Total pattern hits " << tot_patt_hit << std::endl;
    if (patt_comp < _patt_hit_comp_thresh || tot_patt_hit < _patt_hit_thresh)
        return false;

    for (const auto& [_, num_hits] : sig_hit_map) 
    {
        //std::cout << "Signature hit " << num_hits << std::endl;
        if (num_hits / tot_patt_hit < _sig_hit_comp_thresh) 
            return false;       
    }

    return true;
}

DEFINE_ART_CLASS_TOOL(PatternCompleteness)

}

#endif
