#ifndef SELECTION_EMPTYSELECTION_CXX
#define SELECTION_EMPTYSELECTION_CXX

#include "SelectionToolBase.h"

#include <iostream>

namespace selection
{    
    class NeutrinoSelection : public SelectionToolBase {
    public:
        NeutrinoSelection(const fhicl::ParameterSet& pset);
        ~NeutrinoSelection(){};
        
        void configure(fhicl::ParameterSet const & pset);

        bool selectEvent(const art::Event& event, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_vec);

        void setBranches(TTree* _tree){};
        
        void resetTTree(TTree* _tree){};   
    };
   
    NeutrinoSelection::NeutrinoSelection(const fhicl::ParameterSet& pset){}
    
    void NeutrinoSelection::configure(fhicl::ParameterSet const & pset) {}
    
    bool NeutrinoSelection::selectEvent(const art::Event& event,
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_vec) {
        return true;
    }
    
    DEFINE_ART_CLASS_TOOL(NeutrinoSelection)
} 

#endif