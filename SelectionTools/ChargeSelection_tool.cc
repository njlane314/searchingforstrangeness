#include <iostream>
#include <vector>
#include "SelectionToolBase.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"

namespace selection
{
    class ChargeSelection : public SelectionToolBase {
    public:
        ChargeSelection(const fhicl::ParameterSet& pset);
        ~ChargeSelection() {}

        void configure(fhicl::ParameterSet const & pset);
        bool selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& pfp_pxy_v);
        void setBranches(TTree* _tree);
        void resetTTree(TTree* _tree);

    private:
        double _charge_threshold;
        double m_neutrino_slice_charge; 
    };

    ChargeSelection::ChargeSelection(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ChargeSelection::configure(fhicl::ParameterSet const & pset) {
        _charge_threshold = pset.get<double>("charge_threshold");
    }

    void ChargeSelection::setBranches(TTree* _tree) {
        _tree->Branch("neutrino_slice_charge", &m_neutrino_slice_charge, "neutrino_slice_charge/D");
    }

    void ChargeSelection::resetTTree(TTree* _tree) {
        m_neutrino_slice_charge = std::numeric_limits<double>::min();
    }

    bool ChargeSelection::selectEvent(art::Event const& e, const std::vector<ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> slice_hits;
        for (const auto& pfp_pxy : pfp_pxy_v) {
            auto hits = pfp_pxy.get<recob::Hit>();
            slice_hits.insert(slice_hits.end(), hits.begin(), hits.end());
        }

        if (slice_hits.empty()) {
            m_neutrino_slice_charge = 0.0; 
            return false;
        }

        double total_charge = 0.0;
        for (const auto& hit : slice_hits) {
            double charge = hit->Integral();
            total_charge += charge;
        }

        m_neutrino_slice_charge = total_charge;  
        return total_charge > _charge_threshold;
    }
}

DEFINE_ART_CLASS_TOOL(selection::ChargeSelection)