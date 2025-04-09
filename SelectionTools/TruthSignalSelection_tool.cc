#ifndef SELECTION_TRUTHSIGNALSELECTION_CXX
#define SELECTION_TRUTHSIGNALSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "EventClassifier.h"

namespace selection
{
    class TruthSignalSelection : public SelectionToolBase {
    public:
        TruthSignalSelection(const fhicl::ParameterSet& pset);
        ~TruthSignalSelection() = default;
        
        void configure(fhicl::ParameterSet const & pset);
        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;
        
    private:
        std::unique_ptr<signature::EventClassifier> _eventClassifier;
    };
    
    TruthSignalSelection::TruthSignalSelection(const fhicl::ParameterSet& pset) {
        configure(pset);
    }
    
    void TruthSignalSelection::configure(fhicl::ParameterSet const& pset) {
        const auto& classifierPset = pset.get<fhicl::ParameterSet>("EventClassifier");
        _eventClassifier = std::make_unique<signature::EventClassifier>(classifierPset);
    }
    
    bool TruthSignalSelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        if (fData) return false;
        return _eventClassifier->isSignal(e);
    }
    
    void TruthSignalSelection::setBranches(TTree* _tree) {
    }
    
    void TruthSignalSelection::resetTTree(TTree* _tree) {
    }
    
    DEFINE_ART_CLASS_TOOL(TruthSignalSelection)
}

#endif