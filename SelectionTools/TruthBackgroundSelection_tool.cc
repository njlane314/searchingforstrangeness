#ifndef SELECTION_TRUTHBACKGROUNDSELECTION_CXX
#define SELECTION_TRUTHBACKGROUNDSELECTION_CXX

#include <iostream>
#include "SelectionToolBase.h"
#include "../EventClassifier.h"
#include "lardataobj/MCBase/MCShower.h"

namespace selection
{
    class TruthBackgroundSelection : public SelectionToolBase {
    public:
        TruthBackgroundSelection(const fhicl::ParameterSet& pset);
        ~TruthBackgroundSelection() = default;
        
        void configure(fhicl::ParameterSet const & pset);
        bool selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;
        
    private:
        std::unique_ptr<signature::EventClassifier> _eventClassifier;
        bool _selectBeamNeutrino;
        bool _selectCosmicRay;
        bool _selectOther;
    };
    
    TruthBackgroundSelection::TruthBackgroundSelection(const fhicl::ParameterSet& pset) {
        configure(pset);
    }
    
    void TruthBackgroundSelection::configure(fhicl::ParameterSet const& pset) {
        const auto& classifierPset = pset.get<fhicl::ParameterSet>("EventClassifier");
        _eventClassifier = std::make_unique<signature::EventClassifier>(classifierPset);
        
        _selectBeamNeutrino = pset.get<bool>("SelectBeamNeutrino", true);
        _selectCosmicRay = pset.get<bool>("SelectCosmicRay", true);
        _selectOther = pset.get<bool>("SelectOther", true);
    }
    
    bool TruthBackgroundSelection::selectEvent(art::Event const& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        if (fData) return false;
        
        if (_eventClassifier->isSignal(e)) return false;
        
        signature::EventType eventType = _eventClassifier->classifyEvent(e);
        
        switch (eventType) {
            case signature::kBeamNeutrino: return _selectBeamNeutrino;
            case signature::kCosmicRay:    return _selectCosmicRay;
            case signature::kOther:        return _selectOther;
            case signature::kSignal:       return false;
            default:                       return false;
        }
    }
    
    void TruthBackgroundSelection::setBranches(TTree* _tree) {
    }
    
    void TruthBackgroundSelection::resetTTree(TTree* _tree) {
    }
    
    DEFINE_ART_CLASS_TOOL(TruthBackgroundSelection)
}

#endif