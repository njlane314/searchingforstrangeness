#ifndef SELECTION_EMPTYSELECTION_CXX
#define SELECTION_EMPTYSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{
class StrangenessEmptySelection : public SelectionToolBase {

public:

    StrangenessEmptySelection(const fhicl::ParameterSet& pset);
    ~StrangenessEmptySelection(){};
    
    void configure(fhicl::ParameterSet const & pset);

    bool selectEvent(art::Event const& e,
            const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);


    void setBranches(TTree* _tree){};

    void resetTTree(TTree* _tree){};
    
private:
    
};

StrangenessEmptySelection::StrangenessEmptySelection(const fhicl::ParameterSet& pset)
{
}

void StrangenessEmptySelection::configure(fhicl::ParameterSet const & pset)
{
}

bool StrangenessEmptySelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v)
{
    
    return true;
}


DEFINE_ART_CLASS_TOOL(StrangenessEmptySelection)
} 

#endif