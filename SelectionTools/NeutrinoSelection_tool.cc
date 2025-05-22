#ifndef SELECTION_EMPTYSELECTION_CXX
#define SELECTION_EMPTYSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"

namespace selection
{    
    class NeutrinoSelection : public SelectionToolBase {
    public:
        NeutrinoSelection(const fhicl::ParameterSet& pset);
        ~NeutrinoSelection(){};
        
        void configure(fhicl::ParameterSet const & pset);

        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);

        void setBranches(TTree* _tree){};
        
        void resetTTree(TTree* _tree){};   
    };
   
    NeutrinoSelection::NeutrinoSelection(const fhicl::ParameterSet& pset){}
    
    void NeutrinoSelection::configure(fhicl::ParameterSet const & pset) {}
    
    bool NeutrinoSelection::selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) { 
        return true;
    }
    
    DEFINE_ART_CLASS_TOOL(NeutrinoSelection)
} 

#endif