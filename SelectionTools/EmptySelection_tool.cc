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
        bool selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image>& calo_images, 
                                const std::vector<image::Image>& reco_images, 
                                const std::vector<image::Image>& label_images);
        void setBranches(TTree* _tree){};
        void resetTTree(TTree* _tree){};   
    };
   
    EmptySelection::EmptySelection(const fhicl::ParameterSet& pset){}
    
    void EmptySelection::configure(fhicl::ParameterSet const & pset) {}
    
    bool EmptySelection::selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image>& calo_images, 
                                const std::vector<image::Image>& reco_images, 
                                const std::vector<image::Image>& label_images) { 
        return true;
    }
    
    DEFINE_ART_CLASS_TOOL(EmptySelection)
} 

#endif