#ifndef ANALYSIS_NEUTRINOANALYSIS_TOOL_CXX
#define ANALYSIS_NEUTRINOANALYSIS_TOOL_CXX

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"

#include "TTree.h"
#include "TVector3.h"

#include <cmath>
#include <limits>

namespace analysis {

class NeutrinoAnalysis_tool : public AnalysisToolBase {
    public:
    explicit NeutrinoAnalysis_tool(fhicl::ParameterSet const &p) {
        this->configure(p);
    }
    ~NeutrinoAnalysis_tool() override = default;

    void configure(const fhicl::ParameterSet &pset) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &, std::vector<common::ProxyPfpElem_t> &,
                      bool, bool) override {}

    private:
    art::InputTag fMCTproducer;

    int _nu_pdg;
    int _ccnc;
    int _interaction_mode;
    int _interaction_type;
    float _nu_E;
    float _nu_v[3];
    float _nu_dir[3];
    bool _is_nu_mu_cc;

    template <class T> static T nan() {
        return std::numeric_limits<T>::quiet_NaN();
    }
};

void NeutrinoAnalysis_tool::configure(const fhicl::ParameterSet &p) {
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
}

void NeutrinoAnalysis_tool::setBranches(TTree *t) {
    t->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
    t->Branch("ccnc", &_ccnc, "ccnc/I");
    t->Branch("interaction_mode", &_interaction_mode, "interaction_mode/I");
    t->Branch("interaction_type", &_interaction_type, "interaction_type/I");
    t->Branch("nu_E", &_nu_E, "nu_E/F");
    t->Branch("nu_v", _nu_v, "nu_v[3]/F");
    t->Branch("nu_dir", _nu_dir, "nu_dir[3]/F");
    t->Branch("is_nu_mu_cc", &_is_nu_mu_cc, "is_nu_mu_cc/O");
}

void NeutrinoAnalysis_tool::resetTTree(TTree *) {
    _nu_pdg = 0;
    _ccnc = -1;
    _interaction_mode = -1;
    _interaction_type = -1;
    _nu_E = nan<float>();
    _nu_v[0] = _nu_v[1] = _nu_v[2] = nan<float>();
    _nu_dir[0] = _nu_dir[1] = _nu_dir[2] = nan<float>();
    _is_nu_mu_cc = false;
}

void NeutrinoAnalysis_tool::analyseEvent(const art::Event &event, bool is_data) {
    this->resetTTree(nullptr);
    if (is_data)
        return;

    TVector3 nu_dir(0.f, 0.f, 0.f);
    auto const &mct_h =
        event.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    if (!mct_h.isValid() || mct_h->empty())
        return;
    const auto &mct = mct_h->front();
    if (mct.NeutrinoSet()) {
        const auto nuinfo = mct.GetNeutrino();
        const auto &nu = nuinfo.Nu();
        _nu_pdg = nu.PdgCode();
        _ccnc = nuinfo.CCNC();
        _interaction_mode = nuinfo.Mode();
        _interaction_type = nuinfo.InteractionType();
        _nu_E = nu.Trajectory().E(0);
        _nu_v[0] = nu.Vx();
        _nu_v[1] = nu.Vy();
        _nu_v[2] = nu.Vz();

        const double pnu = std::sqrt(nu.Px() * nu.Px() + nu.Py() * nu.Py() +
                                     nu.Pz() * nu.Pz());
        if (pnu > 0.0) {
            nu_dir.SetXYZ(nu.Px() / pnu, nu.Py() / pnu, nu.Pz() / pnu);
            _nu_dir[0] = nu_dir.X();
            _nu_dir[1] = nu_dir.Y();
            _nu_dir[2] = nu_dir.Z();
        }
    }
    _is_nu_mu_cc = (std::abs(_nu_pdg) == 14) && (_ccnc == 0);
}

DEFINE_ART_CLASS_TOOL(NeutrinoAnalysis_tool)

} // namespace analysis

#endif
