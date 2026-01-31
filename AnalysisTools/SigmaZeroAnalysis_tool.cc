#ifndef ANALYSIS_SIGMAZEROANALYSIS_TOOL_CXX
#define ANALYSIS_SIGMAZEROANALYSIS_TOOL_CXX

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"

#include "TTree.h"

#include <cmath>
#include <map>
#include <vector>

namespace analysis {

class SigmaZeroAnalysis_tool : public AnalysisToolBase {
    public:
    explicit SigmaZeroAnalysis_tool(fhicl::ParameterSet const &p) {
        this->configure(p);
    }
    ~SigmaZeroAnalysis_tool() override = default;

    void configure(const fhicl::ParameterSet &pset) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &, std::vector<common::ProxyPfpElem_t> &,
                      bool, bool) override {}

    private:
    art::InputTag fMCPproducer;

    int _n_sigma0_truth;
    int _n_sigma0_to_lambda;
    int _n_sigma0_to_lambda_ppi;
    bool _has_sigma0_truth;
    bool _has_sigma0_to_lambda;
    bool _has_sigma0_to_lambda_ppi;

    struct DecayMatch {
        bool ok = false;
        int p_trkid = -1;
        int pi_trkid = -1;
    };

    static bool IsSigma0(int pdg) { return std::abs(pdg) == 3212; }
    static bool IsLambda(int pdg) { return std::abs(pdg) == 3122; }

    DecayMatch
    MatchLambdaToPPi(const art::Ptr<simb::MCParticle> &lam,
                     const std::map<int, art::Ptr<simb::MCParticle>> &mp) const;
};

void SigmaZeroAnalysis_tool::configure(const fhicl::ParameterSet &p) {
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
}

void SigmaZeroAnalysis_tool::setBranches(TTree *t) {
    t->Branch("n_sigma0_truth", &_n_sigma0_truth, "n_sigma0_truth/I");
    t->Branch("has_sigma0_truth", &_has_sigma0_truth, "has_sigma0_truth/O");
    t->Branch("n_sigma0_to_lambda", &_n_sigma0_to_lambda,
              "n_sigma0_to_lambda/I");
    t->Branch("has_sigma0_to_lambda", &_has_sigma0_to_lambda,
              "has_sigma0_to_lambda/O");
    t->Branch("n_sigma0_to_lambda_ppi", &_n_sigma0_to_lambda_ppi,
              "n_sigma0_to_lambda_ppi/I");
    t->Branch("has_sigma0_to_lambda_ppi", &_has_sigma0_to_lambda_ppi,
              "has_sigma0_to_lambda_ppi/O");
}

void SigmaZeroAnalysis_tool::resetTTree(TTree *) {
    _n_sigma0_truth = 0;
    _n_sigma0_to_lambda = 0;
    _n_sigma0_to_lambda_ppi = 0;
    _has_sigma0_truth = false;
    _has_sigma0_to_lambda = false;
    _has_sigma0_to_lambda_ppi = false;
}

SigmaZeroAnalysis_tool::DecayMatch SigmaZeroAnalysis_tool::MatchLambdaToPPi(
    const art::Ptr<simb::MCParticle> &lam,
    const std::map<int, art::Ptr<simb::MCParticle>> &mp) const {
    DecayMatch ret;
    const int lam_pdg = lam->PdgCode();
    if (!IsLambda(lam_pdg))
        return ret;

    const int want_p = (lam_pdg > 0) ? 2212 : -2212;
    const int want_pi = (lam_pdg > 0) ? -211 : 211;

    int p_tid = -1;
    int pi_tid = -1;
    for (int i = 0; i < lam->NumberDaughters(); ++i) {
        const int d_tid = lam->Daughter(i);
        auto it = mp.find(d_tid);
        if (it == mp.end())
            continue;
        const auto &d = it->second;
        if (d->PdgCode() == want_p)
            p_tid = d->TrackId();
        if (d->PdgCode() == want_pi)
            pi_tid = d->TrackId();
    }

    if (p_tid >= 0 && pi_tid >= 0) {
        ret.ok = true;
        ret.p_trkid = p_tid;
        ret.pi_trkid = pi_tid;
    }
    return ret;
}

void SigmaZeroAnalysis_tool::analyseEvent(const art::Event &event,
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
        const art::Ptr<simb::MCParticle> sz_ptr(mcp_h, i);
        const auto &sz = *sz_ptr;
        if (!IsSigma0(sz.PdgCode()))
            continue;

        ++_n_sigma0_truth;

        art::Ptr<simb::MCParticle> lam_ptr;
        for (int d = 0; d < sz.NumberDaughters(); ++d) {
            const int d_tid = sz.Daughter(d);
            auto it = mp.find(d_tid);
            if (it == mp.end())
                continue;
            if (IsLambda(it->second->PdgCode())) {
                lam_ptr = it->second;
                break;
            }
        }

        if (lam_ptr.isNull())
            continue;
        ++_n_sigma0_to_lambda;
        if (MatchLambdaToPPi(lam_ptr, mp).ok)
            ++_n_sigma0_to_lambda_ppi;
    }

    _has_sigma0_truth = (_n_sigma0_truth > 0);
    _has_sigma0_to_lambda = (_n_sigma0_to_lambda > 0);
    _has_sigma0_to_lambda_ppi = (_n_sigma0_to_lambda_ppi > 0);
}

DEFINE_ART_CLASS_TOOL(SigmaZeroAnalysis_tool)

} // namespace analysis

#endif
