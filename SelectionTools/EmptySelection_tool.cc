#ifndef SELECTION_EMPTYSELECTION_CXX
#define SELECTION_EMPTYSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{    
    class EmptySelection : public SelectionToolBase {
    public:
        EmptySelection(const fhicl::ParameterSet& pset);
        ~EmptySelection(){};
        
        void configure(fhicl::ParameterSet const & pset);

        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);

        void setBranches(TTree* _tree){};
        
        void resetTTree(TTree* _tree){};   
    };
   
    EmptySelection::EmptySelection(const fhicl::ParameterSet& pset){}
    
    void EmptySelection::configure(fhicl::ParameterSet const & pset) {}
    
    bool EmptySelection::selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) { 
        return true;
    }
    
    DEFINE_ART_CLASS_TOOL(EmptySelection)
} 

#endif