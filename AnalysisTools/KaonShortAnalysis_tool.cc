#ifndef ANALYSIS_KAONSHORTANALYSIS_TOOL_CXX
#define ANALYSIS_KAONSHORTANALYSIS_TOOL_CXX

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"

#include "TTree.h"

#include <map>
#include <vector>

namespace analysis {

class KaonShortAnalysis_tool : public AnalysisToolBase {
    public:
    explicit KaonShortAnalysis_tool(fhicl::ParameterSet const &p) {
        this->configure(p);
    }
    ~KaonShortAnalysis_tool() override = default;

    void configure(const fhicl::ParameterSet &pset) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &, std::vector<common::ProxyPfpElem_t> &,
                      bool, bool) override {}

    private:
    art::InputTag fMCPproducer;

    int _n_kshort_truth;
    int _n_kshort_to_pipi;
    bool _has_kshort_truth;
    bool _has_kshort_to_pipi;
};

void KaonShortAnalysis_tool::configure(const fhicl::ParameterSet &p) {
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
}

void KaonShortAnalysis_tool::setBranches(TTree *t) {
    t->Branch("n_kshort_truth", &_n_kshort_truth, "n_kshort_truth/I");
    t->Branch("has_kshort_truth", &_has_kshort_truth, "has_kshort_truth/O");
    t->Branch("n_kshort_to_pipi", &_n_kshort_to_pipi, "n_kshort_to_pipi/I");
    t->Branch("has_kshort_to_pipi", &_has_kshort_to_pipi,
              "has_kshort_to_pipi/O");
}

void KaonShortAnalysis_tool::resetTTree(TTree *) {
    _n_kshort_truth = 0;
    _n_kshort_to_pipi = 0;
    _has_kshort_truth = false;
    _has_kshort_to_pipi = false;
}

void KaonShortAnalysis_tool::analyseEvent(const art::Event &event,
                                          bool is_data) {
    this->resetTTree(nullptr);
    if (is_data)
        return;

    auto const &mcp_h =
        event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mp;
    for (size_t i = 0; i < mcp_h->size(); ++i)
        mp[mcp_h->at(i).TrackId()] = art::Ptr<simb::MCParticle>(mcp_h, i);

    for (size_t i = 0; i < mcp_h->size(); ++i) {
        const art::Ptr<simb::MCParticle> ks_ptr(mcp_h, i);
        const auto &ks = *ks_ptr;
        if (ks.PdgCode() != 310)
            continue;

        ++_n_kshort_truth;

        bool has_pip = false;
        bool has_pim = false;
        for (int d = 0; d < ks.NumberDaughters(); ++d) {
            const int d_tid = ks.Daughter(d);
            auto it = mp.find(d_tid);
            if (it == mp.end())
                continue;
            const int pdg = it->second->PdgCode();
            if (pdg == 211)
                has_pip = true;
            if (pdg == -211)
                has_pim = true;
        }
        if (has_pip && has_pim)
            ++_n_kshort_to_pipi;
    }

    _has_kshort_truth = (_n_kshort_truth > 0);
    _has_kshort_to_pipi = (_n_kshort_to_pipi > 0);
}

DEFINE_ART_CLASS_TOOL(KaonShortAnalysis_tool)

} // namespace analysis

#endif
