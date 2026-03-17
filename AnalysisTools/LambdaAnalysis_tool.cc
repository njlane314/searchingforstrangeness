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

    struct TruthSummary {
        bool is_nc = false;
        bool has_strange_fs = false;
        bool has_nc_strangeness = false;

        bool has_fs_lambda0 = false;
        bool has_fs_sigma0 = false;
        bool has_g4_lambda0_to_ppi = false;

        bool has_ks0 = false;
        bool has_pi0 = false;

        int n_fs_lambda0 = 0;
        int n_fs_sigma0 = 0;
        int n_g4_lambda0_to_ppi = 0;
        int n_ks0 = 0;
        int n_pi0 = 0;

        std::vector<float> fs_lambda0_p;
        std::vector<float> fs_sigma0_p;

        void clear() { *this = TruthSummary{}; }
    };

    struct G4LambdaDecaySummary {
        std::vector<int> trackid;
        std::vector<int> pdg;
        std::vector<float> E;
        std::vector<float> mass;

        std::vector<float> px;
        std::vector<float> py;
        std::vector<float> pz;
        std::vector<float> p_mag;

        std::vector<float> vx;
        std::vector<float> vy;
        std::vector<float> vz;

        std::vector<float> endx;
        std::vector<float> endy;
        std::vector<float> endz;

        std::vector<float> ct;
        std::vector<float> decay_sep;

        std::vector<int> p_trackid;
        std::vector<int> pi_trackid;
        std::vector<float> p_p;
        std::vector<float> pi_p;
        std::vector<float> ppi_opening_angle;

        std::vector<int> pr_valid_assignment;
        std::vector<float> pr_mu_purity;
        std::vector<float> pr_mu_completeness;
        std::vector<float> pr_p_purity;
        std::vector<float> pr_p_completeness;
        std::vector<float> pr_pi_purity;
        std::vector<float> pr_pi_completeness;
        std::vector<int> mu_true_hits;
        std::vector<int> p_true_hits;
        std::vector<int> pi_true_hits;

        void clear() { *this = G4LambdaDecaySummary{}; }
    };

    TruthSummary _truth;
    G4LambdaDecaySummary _g4_lambda;

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
    static bool IsGeneratorFinalState(const simb::MCParticle &p) {
        return p.StatusCode() == 1;
    }
    static bool IsStrangeHadron(const int apdg) {
        switch (apdg) {
        case 130:  // K_L0
        case 310:  // K_S0
        case 311:  // K0 / K0bar
        case 321:  // K+/-
        case 3122: // Lambda0 / anti-Lambda0
        case 3112: // Sigma-
        case 3212: // Sigma0 / anti-Sigma0
        case 3222: // Sigma+
        case 3312: // Xi-
        case 3322: // Xi0
        case 3334: // Omega-
            return true;
        default:
            return false;
        }
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

    void FillGeneratorTruthSummary(const simb::MCTruth &mct);
    void FillG4LambdaDecaySummary(
        const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h,
        const std::map<int, art::Ptr<simb::MCParticle>> &mp,
        const std::array<float, 3> &nu_v);
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

    t->Branch("truth_is_nc", &_truth.is_nc, "truth_is_nc/O");
    t->Branch("truth_has_strange_fs", &_truth.has_strange_fs,
              "truth_has_strange_fs/O");
    t->Branch("truth_has_nc_strangeness", &_truth.has_nc_strangeness,
              "truth_has_nc_strangeness/O");
    t->Branch("truth_has_fs_lambda0", &_truth.has_fs_lambda0,
              "truth_has_fs_lambda0/O");
    t->Branch("truth_has_fs_sigma0", &_truth.has_fs_sigma0,
              "truth_has_fs_sigma0/O");
    t->Branch("truth_has_g4_lambda0_to_ppi", &_truth.has_g4_lambda0_to_ppi,
              "truth_has_g4_lambda0_to_ppi/O");
    t->Branch("truth_has_ks0", &_truth.has_ks0, "truth_has_ks0/O");
    t->Branch("truth_has_pi0", &_truth.has_pi0, "truth_has_pi0/O");

    t->Branch("truth_n_fs_lambda0", &_truth.n_fs_lambda0,
              "truth_n_fs_lambda0/I");
    t->Branch("truth_n_fs_sigma0", &_truth.n_fs_sigma0,
              "truth_n_fs_sigma0/I");
    t->Branch("truth_n_g4_lambda0_to_ppi", &_truth.n_g4_lambda0_to_ppi,
              "truth_n_g4_lambda0_to_ppi/I");
    t->Branch("truth_n_ks0", &_truth.n_ks0, "truth_n_ks0/I");
    t->Branch("truth_n_pi0", &_truth.n_pi0, "truth_n_pi0/I");

    t->Branch("truth_fs_lambda0_p", &_truth.fs_lambda0_p);
    t->Branch("truth_fs_sigma0_p", &_truth.fs_sigma0_p);

    t->Branch("g4_lambda_trackid", &_g4_lambda.trackid);
    t->Branch("g4_lambda_pdg", &_g4_lambda.pdg);
    t->Branch("g4_lambda_E", &_g4_lambda.E);
    t->Branch("g4_lambda_mass", &_g4_lambda.mass);

    t->Branch("g4_lambda_px", &_g4_lambda.px);
    t->Branch("g4_lambda_py", &_g4_lambda.py);
    t->Branch("g4_lambda_pz", &_g4_lambda.pz);
    t->Branch("g4_lambda_p_mag", &_g4_lambda.p_mag);

    t->Branch("g4_lambda_vx", &_g4_lambda.vx);
    t->Branch("g4_lambda_vy", &_g4_lambda.vy);
    t->Branch("g4_lambda_vz", &_g4_lambda.vz);

    t->Branch("g4_lambda_endx", &_g4_lambda.endx);
    t->Branch("g4_lambda_endy", &_g4_lambda.endy);
    t->Branch("g4_lambda_endz", &_g4_lambda.endz);

    t->Branch("g4_lambda_ct", &_g4_lambda.ct);
    t->Branch("g4_lambda_decay_sep", &_g4_lambda.decay_sep);

    t->Branch("g4_lambda_p_trackid", &_g4_lambda.p_trackid);
    t->Branch("g4_lambda_pi_trackid", &_g4_lambda.pi_trackid);
    t->Branch("g4_lambda_p_p", &_g4_lambda.p_p);
    t->Branch("g4_lambda_pi_p", &_g4_lambda.pi_p);
    t->Branch("g4_lambda_ppi_opening_angle", &_g4_lambda.ppi_opening_angle);

    t->Branch("g4_lambda_pr_valid_assignment", &_g4_lambda.pr_valid_assignment);
    t->Branch("g4_lambda_pr_mu_purity", &_g4_lambda.pr_mu_purity);
    t->Branch("g4_lambda_pr_mu_completeness", &_g4_lambda.pr_mu_completeness);
    t->Branch("g4_lambda_pr_p_purity", &_g4_lambda.pr_p_purity);
    t->Branch("g4_lambda_pr_p_completeness", &_g4_lambda.pr_p_completeness);
    t->Branch("g4_lambda_pr_pi_purity", &_g4_lambda.pr_pi_purity);
    t->Branch("g4_lambda_pr_pi_completeness", &_g4_lambda.pr_pi_completeness);
    t->Branch("g4_lambda_mu_true_hits", &_g4_lambda.mu_true_hits);
    t->Branch("g4_lambda_p_true_hits", &_g4_lambda.p_true_hits);
    t->Branch("g4_lambda_pi_true_hits", &_g4_lambda.pi_true_hits);
}

void LambdaAnalysis_tool::resetTTree(TTree *) {
    _mu_truth_trackid = -1;
    _mu_truth_pdg = 0;
    _mu_p = nan<float>();
    _mu_theta = nan<float>();

    _truth.clear();
    _g4_lambda.clear();
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

void LambdaAnalysis_tool::FillGeneratorTruthSummary(const simb::MCTruth &mct) {
    _truth.is_nc =
        mct.NeutrinoSet() && (mct.GetNeutrino().CCNC() == simb::kNC);

    for (int i = 0; i < mct.NParticles(); ++i) {
        const auto &p = mct.GetParticle(i);
        if (!IsGeneratorFinalState(p))
            continue;

        const int apdg = std::abs(p.PdgCode());
        const float pmag = Mag3(p.Px(), p.Py(), p.Pz());

        if (IsStrangeHadron(apdg))
            _truth.has_strange_fs = true;

        if (apdg == 3122) {
            _truth.has_fs_lambda0 = true;
            ++_truth.n_fs_lambda0;
            _truth.fs_lambda0_p.push_back(pmag);
        } else if (apdg == 3212) {
            _truth.has_fs_sigma0 = true;
            ++_truth.n_fs_sigma0;
            _truth.fs_sigma0_p.push_back(pmag);
        } else if (apdg == 310) {
            _truth.has_ks0 = true;
            ++_truth.n_ks0;
        } else if (apdg == 111) {
            _truth.has_pi0 = true;
            ++_truth.n_pi0;
        }
    }

    _truth.has_nc_strangeness = _truth.is_nc && _truth.has_strange_fs;
}

void LambdaAnalysis_tool::FillG4LambdaDecaySummary(
    const art::ValidHandle<std::vector<simb::MCParticle>> &mcp_h,
    const std::map<int, art::Ptr<simb::MCParticle>> &mp,
    const std::array<float, 3> &nu_v) {

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

        float ct = nan<float>();
        const float L = ThreeDistance(lam.EndX(), lam.EndY(), lam.EndZ(),
                                      lam.Vx(), lam.Vy(), lam.Vz());
        const float p_mag = Mag3(lam.Px(), lam.Py(), lam.Pz());
        if (p_mag > 0.f)
            ct = L * lam.Mass() / p_mag;

        float opening_angle = nan<float>();
        if (p && pi) {
            TVector3 vp(p->Px(), p->Py(), p->Pz());
            TVector3 vpi(pi->Px(), pi->Py(), pi->Pz());
            if (vp.Mag() > 0.f && vpi.Mag() > 0.f)
                opening_angle = vp.Angle(vpi);
        }

        _truth.has_g4_lambda0_to_ppi = true;
        ++_truth.n_g4_lambda0_to_ppi;

        _g4_lambda.trackid.push_back(lam.TrackId());
        _g4_lambda.pdg.push_back(lam.PdgCode());
        _g4_lambda.E.push_back(lam.E());
        _g4_lambda.mass.push_back(lam.Mass());

        _g4_lambda.px.push_back(lam.Px());
        _g4_lambda.py.push_back(lam.Py());
        _g4_lambda.pz.push_back(lam.Pz());
        _g4_lambda.p_mag.push_back(p_mag);

        _g4_lambda.vx.push_back(lam.Vx());
        _g4_lambda.vy.push_back(lam.Vy());
        _g4_lambda.vz.push_back(lam.Vz());

        _g4_lambda.endx.push_back(lam.EndX());
        _g4_lambda.endy.push_back(lam.EndY());
        _g4_lambda.endz.push_back(lam.EndZ());

        _g4_lambda.ct.push_back(ct);
        _g4_lambda.decay_sep.push_back(decay_sep);

        _g4_lambda.p_trackid.push_back(dm.p_trkid);
        _g4_lambda.pi_trackid.push_back(dm.pi_trkid);
        _g4_lambda.p_p.push_back(pp);
        _g4_lambda.pi_p.push_back(pip);
        _g4_lambda.ppi_opening_angle.push_back(opening_angle);

        _g4_lambda.pr_valid_assignment.push_back(0);
        _g4_lambda.pr_mu_purity.push_back(nan<float>());
        _g4_lambda.pr_mu_completeness.push_back(nan<float>());
        _g4_lambda.pr_p_purity.push_back(nan<float>());
        _g4_lambda.pr_p_completeness.push_back(nan<float>());
        _g4_lambda.pr_pi_purity.push_back(nan<float>());
        _g4_lambda.pr_pi_completeness.push_back(nan<float>());
        _g4_lambda.mu_true_hits.push_back(-1);
        _g4_lambda.p_true_hits.push_back(-1);
        _g4_lambda.pi_true_hits.push_back(-1);
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
    FillGeneratorTruthSummary(mct);
    FillG4LambdaDecaySummary(mcp_h, mp, nu_v);
}

void LambdaAnalysis_tool::analyseSlice(
    const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
    bool is_data, bool) {
    if (is_data)
        return;
    if (_g4_lambda.trackid.empty())
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

    int mu_true_hits = -1;
    if (_mu_truth_trackid >= 0) {
        mu_true_hits = static_cast<int>(common::CountTruthHitsInSlice(
            _mu_truth_trackid, inputHits, assocMCPart));
    }

    for (size_t i = 0; i < _g4_lambda.trackid.size(); ++i) {
        const int p_tid = _g4_lambda.p_trackid[i];
        const int pi_tid = _g4_lambda.pi_trackid[i];

        _g4_lambda.mu_true_hits[i] = mu_true_hits;
        if (p_tid >= 0) {
            _g4_lambda.p_true_hits[i] = static_cast<int>(
                common::CountTruthHitsInSlice(p_tid, inputHits, assocMCPart));
        }
        if (pi_tid >= 0) {
            _g4_lambda.pi_true_hits[i] = static_cast<int>(
                common::CountTruthHitsInSlice(pi_tid, inputHits, assocMCPart));
        }

        if (_mu_truth_trackid < 0 || p_tid < 0 || pi_tid < 0)
            continue;

        const std::vector<int> tids = {_mu_truth_trackid, p_tid, pi_tid};
        const auto metrics = common::ComputePRMetrics(
            event, fCLSproducer, fHITproducer, fBKTproducer, slice_pfp_vec,
            tids);
        if (!metrics.valid)
            continue;
        if (metrics.purity.size() < 3 || metrics.completeness.size() < 3)
            continue;

        _g4_lambda.pr_valid_assignment[i] = 1;
        _g4_lambda.pr_mu_purity[i] = metrics.purity[0];
        _g4_lambda.pr_mu_completeness[i] = metrics.completeness[0];
        _g4_lambda.pr_p_purity[i] = metrics.purity[1];
        _g4_lambda.pr_p_completeness[i] = metrics.completeness[1];
        _g4_lambda.pr_pi_purity[i] = metrics.purity[2];
        _g4_lambda.pr_pi_completeness[i] = metrics.completeness[2];
    }
}

DEFINE_ART_CLASS_TOOL(LambdaAnalysis_tool)

} // namespace analysis

#endif
