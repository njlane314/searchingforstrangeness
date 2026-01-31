#ifndef ANALYSIS_LAMBDAANALYSIS_TOOL_CXX
#define ANALYSIS_LAMBDAANALYSIS_TOOL_CXX

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"
#include "Common/PatternRecognitionUtils.h"

#include "TTree.h"
#include "TVector3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {
} // namespace

namespace analysis {

class LambdaAnalysis_tool : public AnalysisToolBase {
    public:
    explicit LambdaAnalysis_tool(fhicl::ParameterSet const &p) {
        this->configure(p);
    }
    ~LambdaAnalysis_tool() override = default;

    void configure(const fhicl::ParameterSet &pset) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event,
                      std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
                      bool is_data, bool is_selected) override;

    private:
    art::InputTag fMCTproducer;
    art::InputTag fMCPproducer;
    art::InputTag fCLSproducer;
    art::InputTag fHITproducer;
    art::InputTag fBKTproducer;

    int _mu_truth_trackid;
    int _mu_truth_pdg;
    float _mu_p;
    float _mu_theta;
    int _lam_trackid;
    int _lam_pdg;
    float _lam_E;
    float _lam_mass;
    float _lam_p[3];
    float _lam_v[3];
    float _lam_end[3];
    float _lam_ct;
    float _lam_decay_sep;

    int _p_trackid;
    int _pi_trackid;
    float _p_p;
    float _pi_p;
    float _ppi_opening_angle;

    bool _pr_valid_assignment;
    float _pr_mu_purity;
    float _pr_mu_completeness;
    float _pr_p_purity;
    float _pr_p_completeness;
    float _pr_pi_purity;
    float _pr_pi_completeness;
    int _mu_true_hits;
    int _p_true_hits;
    int _pi_true_hits;

    template <class T> static T nan() {
        return std::numeric_limits<T>::quiet_NaN();
    }
    static float Mag3(const float x, const float y, const float z) {
        return std::sqrt(x * x + y * y + z * z);
    }
    static float ThreeDistance(float x1, float y1, float z1, float x2, float y2,
                               float z2) {
        return Mag3(x1 - x2, y1 - y2, z1 - z2);
    }

    struct DecayMatch {
        bool ok = false;
        int p_trkid = -1;
        int pi_trkid = -1;
    };

    DecayMatch
    MatchLambdaToPPi(const art::Ptr<simb::MCParticle> &lam,
                     const std::map<int, art::Ptr<simb::MCParticle>> &mp) const;

    void
    FindTruthMuon(const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h,
                  const TVector3 &nu_dir);
};

void LambdaAnalysis_tool::configure(const fhicl::ParameterSet &p) {
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fBKTproducer = p.get<art::InputTag>("BKTproducer");
}

void LambdaAnalysis_tool::setBranches(TTree *t) {
    t->Branch("mu_truth_trackid", &_mu_truth_trackid, "mu_truth_trackid/I");
    t->Branch("mu_truth_pdg", &_mu_truth_pdg, "mu_truth_pdg/I");
    t->Branch("mu_p", &_mu_p, "mu_p/F");
    t->Branch("mu_theta", &_mu_theta, "mu_theta/F");
    t->Branch("lam_trackid", &_lam_trackid, "lam_trackid/I");
    t->Branch("lam_pdg", &_lam_pdg, "lam_pdg/I");
    t->Branch("lam_E", &_lam_E, "lam_E/F");
    t->Branch("lam_mass", &_lam_mass, "lam_mass/F");
    t->Branch("lam_p", _lam_p, "lam_p[3]/F");
    t->Branch("lam_v", _lam_v, "lam_v[3]/F");
    t->Branch("lam_end", _lam_end, "lam_end[3]/F");
    t->Branch("lam_ct", &_lam_ct, "lam_ct/F");
    t->Branch("lam_decay_sep", &_lam_decay_sep, "lam_decay_sep/F");

    t->Branch("p_trackid", &_p_trackid, "p_trackid/I");
    t->Branch("pi_trackid", &_pi_trackid, "pi_trackid/I");
    t->Branch("p_p", &_p_p, "p_p/F");
    t->Branch("pi_p", &_pi_p, "pi_p/F");
    t->Branch("ppi_opening_angle", &_ppi_opening_angle, "ppi_opening_angle/F");

    t->Branch("pr_valid_assignment", &_pr_valid_assignment,
              "pr_valid_assignment/O");
    t->Branch("pr_mu_purity", &_pr_mu_purity, "pr_mu_purity/F");
    t->Branch("pr_mu_completeness", &_pr_mu_completeness,
              "pr_mu_completeness/F");
    t->Branch("pr_p_purity", &_pr_p_purity, "pr_p_purity/F");
    t->Branch("pr_p_completeness", &_pr_p_completeness,
              "pr_p_completeness/F");
    t->Branch("pr_pi_purity", &_pr_pi_purity, "pr_pi_purity/F");
    t->Branch("pr_pi_completeness", &_pr_pi_completeness,
              "pr_pi_completeness/F");
    t->Branch("mu_true_hits", &_mu_true_hits, "mu_true_hits/I");
    t->Branch("p_true_hits", &_p_true_hits, "p_true_hits/I");
    t->Branch("pi_true_hits", &_pi_true_hits, "pi_true_hits/I");
}

void LambdaAnalysis_tool::resetTTree(TTree *) {
    _mu_truth_trackid = -1;
    _mu_truth_pdg = 0;
    _mu_p = nan<float>();
    _mu_theta = nan<float>();
    _lam_trackid = -1;
    _lam_pdg = 0;
    _lam_E = nan<float>();
    _lam_mass = nan<float>();
    _lam_p[0] = _lam_p[1] = _lam_p[2] = nan<float>();
    _lam_v[0] = _lam_v[1] = _lam_v[2] = nan<float>();
    _lam_end[0] = _lam_end[1] = _lam_end[2] = nan<float>();
    _lam_ct = nan<float>();
    _lam_decay_sep = nan<float>();


    _p_trackid = -1;
    _pi_trackid = -1;
    _p_p = nan<float>();
    _pi_p = nan<float>();
    _ppi_opening_angle = nan<float>();

    _pr_valid_assignment = false;
    _pr_mu_purity = nan<float>();
    _pr_mu_completeness = nan<float>();
    _pr_p_purity = nan<float>();
    _pr_p_completeness = nan<float>();
    _pr_pi_purity = nan<float>();
    _pr_pi_completeness = nan<float>();
    _mu_true_hits = -1;
    _p_true_hits = -1;
    _pi_true_hits = -1;
}

LambdaAnalysis_tool::DecayMatch LambdaAnalysis_tool::MatchLambdaToPPi(
    const art::Ptr<simb::MCParticle> &lam,
    const std::map<int, art::Ptr<simb::MCParticle>> &mp) const {
    DecayMatch ret;
    const int lam_pdg = lam->PdgCode();
    if (std::abs(lam_pdg) != 3122)
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

void LambdaAnalysis_tool::FindTruthMuon(
    const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h,
    const TVector3 &nu_dir) {
    _mu_truth_trackid = -1;
    _mu_truth_pdg = 0;
    _mu_p = nan<float>();
    _mu_theta = nan<float>();
    double bestE = -1.0;
    float best_px = 0.f, best_py = 0.f, best_pz = 0.f;
    for (size_t i = 0; i < mcp_h->size(); ++i) {
        const auto &p = mcp_h->at(i);
        if (std::abs(p.PdgCode()) != 13)
            continue;
        if (!(p.Process() == "primary" && p.StatusCode() == 1))
            continue;
        if (p.E() <= bestE)
            continue;
        bestE = p.E();
        _mu_truth_trackid = p.TrackId();
        _mu_truth_pdg = p.PdgCode();
        best_px = p.Px();
        best_py = p.Py();
        best_pz = p.Pz();
    }

    if (_mu_truth_trackid < 0)
        return;

    _mu_p = Mag3(best_px, best_py, best_pz);
    if (nu_dir.Mag() > 0.f && _mu_p > 0.f) {
        TVector3 mu_p(best_px, best_py, best_pz);
        _mu_theta = nu_dir.Angle(mu_p);
    }
}

void LambdaAnalysis_tool::analyseEvent(const art::Event &event, bool is_data) {
    this->resetTTree(nullptr);
    if (is_data)
        return;

    std::array<float, 3> nu_v = {nan<float>(), nan<float>(), nan<float>()};
    TVector3 nu_dir(0.f, 0.f, 0.f);
    auto const &mct_h =
        event.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
    if (!mct_h.isValid() || mct_h->empty())
        return;
    const auto &mct = mct_h->front();
    if (mct.NeutrinoSet()) {
        const auto nuinfo = mct.GetNeutrino();
        const auto &nu = nuinfo.Nu();
        nu_v[0] = nu.Vx();
        nu_v[1] = nu.Vy();
        nu_v[2] = nu.Vz();

        const double pnu = std::sqrt(nu.Px() * nu.Px() + nu.Py() * nu.Py() +
                                     nu.Pz() * nu.Pz());
        if (pnu > 0.0) {
            nu_dir.SetXYZ(nu.Px() / pnu, nu.Py() / pnu, nu.Pz() / pnu);
        }
    }

    auto const &mcp_h =
        event.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);
    std::map<int, art::Ptr<simb::MCParticle>> mp;
    for (size_t i = 0; i < mcp_h->size(); ++i)
        mp[mcp_h->at(i).TrackId()] = art::Ptr<simb::MCParticle>(mcp_h, i);

    FindTruthMuon(mcp_h, nu_dir);

    bool have_sel = false;
    float best_sep = -1.f;
    float best_E = -std::numeric_limits<float>::max();

    for (size_t i = 0; i < mcp_h->size(); ++i) {
        const art::Ptr<simb::MCParticle> lam_ptr(mcp_h, i);
        const auto &lam = *lam_ptr;
        if (std::abs(lam.PdgCode()) != 3122)
            continue;
        const DecayMatch dm = MatchLambdaToPPi(lam_ptr, mp);
        if (!dm.ok)
            continue;

        const simb::MCParticle *p =
            (mp.count(dm.p_trkid) ? mp.at(dm.p_trkid).get() : nullptr);
        const simb::MCParticle *pi =
            (mp.count(dm.pi_trkid) ? mp.at(dm.pi_trkid).get() : nullptr);
        const float pp = p ? Mag3(p->Px(), p->Py(), p->Pz()) : nan<float>();
        const float pip =
            pi ? Mag3(pi->Px(), pi->Py(), pi->Pz()) : nan<float>();

        float decay_sep = nan<float>();
        if (std::isfinite(nu_v[0]) && std::isfinite(nu_v[1]) &&
            std::isfinite(nu_v[2])) {
            decay_sep = ThreeDistance(lam.EndX(), lam.EndY(), lam.EndZ(),
                                      nu_v[0], nu_v[1], nu_v[2]);
        }
        const float sep_rank = std::isfinite(decay_sep) ? decay_sep : -1.f;
        const float E = lam.E();

        const bool better =
            (!have_sel) || (sep_rank > best_sep) ||
            (sep_rank == best_sep && E > best_E);

        if (!better)
            continue;

        have_sel = true;
        best_sep = sep_rank;
        best_E = E;

        _lam_trackid = lam.TrackId();
        _lam_pdg = lam.PdgCode();
        _lam_E = lam.E();
        _lam_mass = lam.Mass();
        _lam_p[0] = lam.Px();
        _lam_p[1] = lam.Py();
        _lam_p[2] = lam.Pz();
        _lam_v[0] = lam.Vx();
        _lam_v[1] = lam.Vy();
        _lam_v[2] = lam.Vz();
        _lam_end[0] = lam.EndX();
        _lam_end[1] = lam.EndY();
        _lam_end[2] = lam.EndZ();
        _lam_decay_sep = decay_sep;

        const float L = ThreeDistance(lam.EndX(), lam.EndY(), lam.EndZ(),
                                      lam.Vx(), lam.Vy(), lam.Vz());
        const float p_mag = Mag3(lam.Px(), lam.Py(), lam.Pz());
        _lam_ct = (p_mag > 0.f) ? (L * lam.Mass() / p_mag) : nan<float>();

        _p_trackid = dm.p_trkid;
        _pi_trackid = dm.pi_trkid;
        _p_p = pp;
        _pi_p = pip;

        _ppi_opening_angle = nan<float>();
        if (p && pi) {
            TVector3 vp(p->Px(), p->Py(), p->Pz());
            TVector3 vpi(pi->Px(), pi->Py(), pi->Pz());
            if (vp.Mag() > 0.f && vpi.Mag() > 0.f)
                _ppi_opening_angle = vp.Angle(vpi);
        }
    }
}

void LambdaAnalysis_tool::analyseSlice(
    const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
    bool is_data, bool) {
    if (is_data)
        return;

    _pr_valid_assignment = false;

    _pr_mu_purity = nan<float>();
    _pr_mu_completeness = nan<float>();
    _pr_p_purity = nan<float>();
    _pr_p_completeness = nan<float>();
    _pr_pi_purity = nan<float>();
    _pr_pi_completeness = nan<float>();

    if (_mu_truth_trackid < 0 || _p_trackid < 0 || _pi_trackid < 0)
        return;

    auto const &cluster_h =
        event.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
    art::FindManyP<recob::Hit> assocHits(cluster_h, event, fCLSproducer);

    std::vector<art::Ptr<recob::Hit>> inputHits;
    inputHits.reserve(1024);

    for (size_t ip = 0; ip < slice_pfp_vec.size(); ++ip) {
        auto clusters = slice_pfp_vec[ip].get<recob::Cluster>();
        for (auto const &c : clusters) {
            auto hits = assocHits.at(c.key());
            inputHits.insert(inputHits.end(), hits.begin(), hits.end());
        }
    }

    auto const &hit_h =
        event.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    auto assocMCPart = std::make_unique<
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
        hit_h, event, fBKTproducer);

    _mu_true_hits = static_cast<int>(common::CountTruthHitsInSlice(
        _mu_truth_trackid, inputHits, assocMCPart));
    _p_true_hits = static_cast<int>(
        common::CountTruthHitsInSlice(_p_trackid, inputHits, assocMCPart));
    _pi_true_hits = static_cast<int>(
        common::CountTruthHitsInSlice(_pi_trackid, inputHits, assocMCPart));

    const std::vector<int> tids = {_mu_truth_trackid, _p_trackid, _pi_trackid};
    const auto metrics = common::ComputePRMetrics(
        event, fCLSproducer, fHITproducer, fBKTproducer, slice_pfp_vec, tids);
    if (!metrics.valid)
        return;

    _pr_valid_assignment = metrics.valid;
    _pr_mu_purity = metrics.purity[0];
    _pr_p_purity = metrics.purity[1];
    _pr_pi_purity = metrics.purity[2];
    _pr_mu_completeness = metrics.completeness[0];
    _pr_p_completeness = metrics.completeness[1];
    _pr_pi_completeness = metrics.completeness[2];

}

DEFINE_ART_CLASS_TOOL(LambdaAnalysis_tool)

} // namespace analysis

#endif
