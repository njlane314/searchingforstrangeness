#ifndef ANALYSIS_VERTEXTOPOLOGY_CXX
#define ANALYSIS_VERTEXTOPOLOGY_CXX

#include "AnalysisToolBase.h"
#include "Common/ProxyTypes.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "TVector3.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>

namespace analysis {

class VertexTopology : public AnalysisToolBase {
public:
    explicit VertexTopology(const fhicl::ParameterSet &pset);
    ~VertexTopology() { };

    void configure(fhicl::ParameterSet const &pset) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event,
                      std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
                      bool is_data, bool is_selected) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;

private:
    // ---------- Original "VertexCleanness" metrics ----------
    void compute_vtx_scores(const std::vector<TVector3>& dirs,
                            const std::vector<float>& weights,
                            const TVector3& beam_dir,
                            float& back_frac,
                            float& off_frac,
                            float& score);

    TVector3 fBNBdir;
    TVector3 fNuMIdir;

    float fVtxRadius;
    float fFwdCos;
    float fBackRMin;
    float fBackRMax;
    float fBackMargin;
    float fBackQmax;
    float fVtxQmax;
    float fAlpha;
    float fBeta;

    float _vtx_score_bnb;
    float _vtx_score_numi;
    float _vtx_backfrac_bnb;
    float _vtx_backfrac_numi;
    float _vtx_offfrac_bnb;
    float _vtx_offfrac_numi;

    // ---------- HadFlow+ event-shape around the vertex ----------
    // Config
    float fKernelR;       // cm, soft radial kernel scale
    float fRhoRef;        // reference density for normalization
    float fMuParTarget;   // target <|cosθ|> before penalizing (≈0.4)
    float fFwdLambda;     // penalty strength (≈0.2)
    float fWThrust;       // weight for thrust deficit (≈0.35)
    float fWSpher;        // weight for sphericity (≈0.35)
    float fWRho;          // weight for density term (≈0.30)
    int   fThrustIters;   // iterations for thrust maximization

    // Outputs (beam-independent primitives)
    float _had_thrust_def;    // 1 - thrust
    float _had_sphericity;    // sphericity
    float _had_rho_term;      // normalized local activity

    // Beam-dependent pieces and final scores
    float _had_mu_parallel_bnb;
    float _had_mu_parallel_numi;
    float _had_fwd_penalty_bnb;
    float _had_fwd_penalty_numi;
    float _hadflow_bnb;
    float _hadflow_numi;

private:
    // HadFlow core
    void compute_hadflow(const std::vector<TVector3>& u,
                         const std::vector<float>& wtilde,
                         const TVector3& beam_dir,
                         float& mu_parallel,
                         float& fwd_penalty,
                         float& score) const;

    // Helpers
    inline float kernel(float r) const {
        const float x = r / std::max(1e-6f, fKernelR);
        return std::exp(-x*x);
    }
    static inline float clamp01(float x) {
        return (x < 0.f) ? 0.f : (x > 1.f ? 1.f : x);
    }

    float thrust_deficit(const std::vector<TVector3>& u,
                         const std::vector<float>& w,
                         int iters) const;

    float sphericity(const std::vector<TVector3>& u,
                     const std::vector<float>& w) const;

    float rho_term(const std::vector<float>& w) const;
};

// ===== Implementation =====

VertexTopology::VertexTopology(const fhicl::ParameterSet &pset) { configure(pset); }

void VertexTopology::configure(fhicl::ParameterSet const &pset)
{
    // Beam directions
    auto bnb = pset.get<std::vector<float>>("BNBBeamDir", {1,0,0});
    if (bnb.size() == 3) fBNBdir = TVector3(bnb[0], bnb[1], bnb[2]).Unit();
    auto numi = pset.get<std::vector<float>>("NuMIBeamDir", {0,0,1});
    if (numi.size() == 3) fNuMIdir = TVector3(numi[0], numi[1], numi[2]).Unit();

    // --- Original cleanness defaults ---
    fVtxRadius = pset.get<float>("VertexRadius", 5.f);
    float theta = pset.get<float>("ForwardAngleDeg", 25.f);
    fFwdCos = std::cos(theta * TMath::DegToRad());
    fBackRMin = pset.get<float>("BackwardRadiusMin", 3.f);
    fBackRMax = pset.get<float>("BackwardRadiusMax", 7.f);
    fBackMargin = pset.get<float>("BackwardCosMargin", 0.05f);
    fBackQmax = pset.get<float>("BackwardMax", std::numeric_limits<float>::max());
    fVtxQmax = pset.get<float>("VertexResidualMax", std::numeric_limits<float>::max());
    fAlpha = pset.get<float>("Alpha", 1.f);
    fBeta  = pset.get<float>("Beta",  1.f);

    // --- HadFlow+ defaults ---
    fKernelR     = pset.get<float>("KernelR",        25.f);   // cm
    fRhoRef      = pset.get<float>("RhoRef",         50.f);   // tune on sidebands
    fMuParTarget = pset.get<float>("MuParTarget",     0.4f);  // ⟨|cosθ|⟩ target
    fFwdLambda   = pset.get<float>("FwdLambda",       0.2f);  // penalty strength
    fWThrust     = pset.get<float>("WThrust",         0.35f);
    fWSpher      = pset.get<float>("WSpher",          0.35f);
    fWRho        = pset.get<float>("WRho",            0.30f);
    fThrustIters = pset.get<int>("ThrustIters",       6);
}

void VertexTopology::analyseEvent(const art::Event &, bool) { }

void VertexTopology::analyseSlice(const art::Event &,
                                  std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
                                  bool, bool)
{
    // --- Locate neutrino vertex (same as your original) ---
    TVector3 vtx;
    bool has_vtx = false;
    for (auto const &pfp : slice_pfp_vec) {
        int pdg = std::abs(pfp->PdgCode());
        if (pdg == 12 || pdg == 14) {
            auto vtx_v = pfp.get<recob::Vertex>();
            if (vtx_v.size() == 1) {
                double xyz[3]; vtx_v[0]->XYZ(xyz);
                vtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
                has_vtx = true;
            }
            break;
        }
    }
    if (!has_vtx) return;

    // --- Gather displacement vectors and weights (1.0 each; replace with charge if desired) ---
    std::vector<TVector3> dirs;
    std::vector<float>    weights;
    dirs.reserve(1024); weights.reserve(1024);

    for (auto const &pfp : slice_pfp_vec) {
        int pdg = std::abs(pfp->PdgCode());
        if (pdg == 12 || pdg == 14) continue; // skip neutrino PFP
        auto sp_v = pfp.get<recob::SpacePoint>();
        for (auto const &sp : sp_v) {
            double xyz[3]; sp->XYZ(xyz);
            TVector3 r(xyz[0]-vtx.X(), xyz[1]-vtx.Y(), xyz[2]-vtx.Z());
            if (r.Mag() < 1e-6) continue;
            dirs.push_back(r);
            weights.push_back(1.f); // TODO: sum associated hit charge instead, if available
        }
    }

    // --- Original vertex-cleanness summaries (BNB / NuMI) ---
    compute_vtx_scores(dirs, weights, fBNBdir,
                       _vtx_backfrac_bnb, _vtx_offfrac_bnb, _vtx_score_bnb);
    compute_vtx_scores(dirs, weights, fNuMIdir,
                       _vtx_backfrac_numi, _vtx_offfrac_numi, _vtx_score_numi);

    // --- HadFlow+ primitives (make unit directions + kernel-weighted weights) ---
    std::vector<TVector3> u;        u.reserve(dirs.size());
    std::vector<float>    wtilde;   wtilde.reserve(dirs.size());

    for (size_t i = 0; i < dirs.size(); ++i) {
        const TVector3& r = dirs[i];
        const double R = r.Mag();
        const float  k = kernel(R);
        u.emplace_back(r.Unit());
        // multiply by any upstream weight (currently 1.f)
        wtilde.emplace_back(k * weights[i]);
    }

    if (wtilde.empty()) {
        _had_thrust_def = 0.f;
        _had_sphericity = 0.f;
        _had_rho_term   = 0.f;

        _had_mu_parallel_bnb  = 0.f;
        _had_mu_parallel_numi = 0.f;
        _had_fwd_penalty_bnb  = 0.f;
        _had_fwd_penalty_numi = 0.f;
        _hadflow_bnb = 0.f;
        _hadflow_numi = 0.f;
        return;
    }

    _had_thrust_def = thrust_deficit(u, wtilde, fThrustIters);
    _had_sphericity = sphericity(u, wtilde);
    _had_rho_term   = rho_term(wtilde);

    // Beam-aware HadFlow scores
    compute_hadflow(u, wtilde, fBNBdir,
                    _had_mu_parallel_bnb, _had_fwd_penalty_bnb, _hadflow_bnb);
    compute_hadflow(u, wtilde, fNuMIdir,
                    _had_mu_parallel_numi, _had_fwd_penalty_numi, _hadflow_numi);
}

// ---------- Original vertex-cleanness computation ----------
void VertexTopology::compute_vtx_scores(const std::vector<TVector3>& dirs,
                                        const std::vector<float>& weights,
                                        const TVector3& beam_dir,
                                        float& back_frac,
                                        float& off_frac,
                                        float& score)
{
    float sum_ann = 0.f, sum_back = 0.f;
    float sum_vtx = 0.f, sum_vtx_off = 0.f;

    for (size_t i = 0; i < dirs.size(); ++i) {
        const TVector3& r = dirs[i];
        const float w = weights[i];
        const TVector3 u = r.Unit();
        const float proj = u.Dot(beam_dir);
        const float mag = r.Mag();

        if (mag >= fBackRMin && mag <= fBackRMax) {
            sum_ann += w;
            if (proj < -fBackMargin) sum_back += w;
        }

        if (mag <= fVtxRadius) {
            sum_vtx += w;
            if (proj < fFwdCos) sum_vtx_off += w;
        }
    }

    back_frac = (sum_ann > 0) ? sum_back / sum_ann : 0.f;
    if (sum_back > fBackQmax) back_frac = 1.f;

    off_frac = (sum_vtx > 0) ? sum_vtx_off / sum_vtx : 0.f;
    if (sum_vtx_off > fVtxQmax) off_frac = 1.f;

    score = fAlpha * back_frac + fBeta * off_frac;
}

// ---------- HadFlow helpers ----------
float VertexTopology::thrust_deficit(const std::vector<TVector3>& u,
                                     const std::vector<float>& w,
                                     int iters) const
{
    double sumW = 0.0; TVector3 n0(0,0,0);
    for (size_t i = 0; i < u.size(); ++i) { sumW += w[i]; n0 += w[i]*u[i]; }
    if (sumW <= 0.0) return 0.f;
    TVector3 n = (n0.Mag2() > 0) ? n0.Unit() : TVector3(1,0,0);

    for (int k = 0; k < iters; ++k) {
        TVector3 s(0,0,0);
        for (size_t i = 0; i < u.size(); ++i) {
            const double sign = (u[i].Dot(n) >= 0.0) ? 1.0 : -1.0;
            s += (w[i]*sign) * u[i];
        }
        if (s.Mag2() < 1e-12) break;
        n = s.Unit();
    }

    double numer = 0.0;
    for (size_t i = 0; i < u.size(); ++i)
        numer += w[i] * std::abs(u[i].Dot(n));
    const double T = (sumW > 0.0) ? (numer / sumW) : 1.0;
    return static_cast<float>(clamp01(1.0 - T));
}

float VertexTopology::sphericity(const std::vector<TVector3>& u,
                                 const std::vector<float>& w) const
{
    double sumW = 0.0;
    TMatrixDSym M(3); M.Zero();
    for (size_t i = 0; i < u.size(); ++i) {
        const double wi = w[i];
        sumW += wi;
        const double x = u[i].X(), y = u[i].Y(), z = u[i].Z();
        M(0,0) += wi*x*x; M(0,1) += wi*x*y; M(0,2) += wi*x*z;
        M(1,1) += wi*y*y; M(1,2) += wi*y*z;
        M(2,2) += wi*z*z;
    }
    if (sumW <= 0.0) return 0.f;
    M(1,0) = M(0,1); M(2,0) = M(0,2); M(2,1) = M(1,2);
    M *= (1.0 / sumW);

    TMatrixDSymEigen eig(M);
    TVectorD vals = eig.GetEigenValues(); // ascending
    std::array<double,3> lam{vals[0], vals[1], vals[2]};
    std::sort(lam.begin(), lam.end(), std::greater<double>());

    const double S = 1.5 * (lam[1] + lam[2]); // [0,1]
    return static_cast<float>(clamp01(S));
}

float VertexTopology::rho_term(const std::vector<float>& w) const
{
    double sumW = 0.0;
    for (float wi : w) sumW += wi;
    const double vol = (4.0/3.0) * TMath::Pi() * std::pow(fKernelR, 3); // effective ball volume at scale R
    const double rho = (vol > 0.0) ? (sumW / vol) : 0.0;
    const double term = (fRhoRef > 0.0) ? std::min(rho / fRhoRef, 1.0) : 0.0;
    return static_cast<float>(clamp01(term));
}

void VertexTopology::compute_hadflow(const std::vector<TVector3>& u,
                                     const std::vector<float>& wtilde,
                                     const TVector3& beam_dir,
                                     float& mu_parallel,
                                     float& fwd_penalty,
                                     float& score) const
{
    if (wtilde.empty()) {
        mu_parallel = fwd_penalty = score = 0.f;
        return;
    }
    double sumW = 0.0, sumWabs = 0.0;
    for (size_t i = 0; i < u.size(); ++i) {
        const double wi = wtilde[i];
        sumW += wi;
        sumWabs += wi * std::abs(u[i].Dot(beam_dir));
    }
    mu_parallel = (sumW > 0.0) ? static_cast<float>(sumWabs / sumW) : 0.f;

    const float P = std::max(0.f, mu_parallel - fMuParTarget);
    fwd_penalty = P;

    const float raw =
        fWThrust * _had_thrust_def +
        fWSpher  * _had_sphericity +
        fWRho    * _had_rho_term -
        fFwdLambda * fwd_penalty;

    score = clamp01(raw);
}

// ---------- I/O ----------
void VertexTopology::setBranches(TTree *t)
{
    // Original cleanness (kept for compatibility)
    t->Branch("vtx_score_bnb",     &_vtx_score_bnb,     "vtx_score_bnb/F");
    t->Branch("vtx_score_numi",    &_vtx_score_numi,    "vtx_score_numi/F");
    t->Branch("vtx_backfrac_bnb",  &_vtx_backfrac_bnb,  "vtx_backfrac_bnb/F");
    t->Branch("vtx_backfrac_numi", &_vtx_backfrac_numi, "vtx_backfrac_numi/F");
    t->Branch("vtx_offfrac_bnb",   &_vtx_offfrac_bnb,   "vtx_offfrac_bnb/F");
    t->Branch("vtx_offfrac_numi",  &_vtx_offfrac_numi,  "vtx_offfrac_numi/F");

    // HadFlow+ primitives and scores
    t->Branch("had_thrust_def",      &_had_thrust_def,      "had_thrust_def/F");
    t->Branch("had_sphericity",      &_had_sphericity,      "had_sphericity/F");
    t->Branch("had_rho_term",        &_had_rho_term,        "had_rho_term/F");

    t->Branch("had_mu_parallel_bnb", &_had_mu_parallel_bnb, "had_mu_parallel_bnb/F");
    t->Branch("had_mu_parallel_numi",&_had_mu_parallel_numi,"had_mu_parallel_numi/F");
    t->Branch("had_fwd_penalty_bnb", &_had_fwd_penalty_bnb, "had_fwd_penalty_bnb/F");
    t->Branch("had_fwd_penalty_numi",&_had_fwd_penalty_numi,"had_fwd_penalty_numi/F");

    t->Branch("hadflow_bnb",         &_hadflow_bnb,         "hadflow_bnb/F");
    t->Branch("hadflow_numi",        &_hadflow_numi,        "hadflow_numi/F");
}

void VertexTopology::resetTTree(TTree *)
{
    const float NaN = std::numeric_limits<float>::quiet_NaN();

    _vtx_score_bnb = NaN;   _vtx_score_numi = NaN;
    _vtx_backfrac_bnb = NaN; _vtx_backfrac_numi = NaN;
    _vtx_offfrac_bnb = NaN;  _vtx_offfrac_numi = NaN;

    _had_thrust_def = NaN; _had_sphericity = NaN; _had_rho_term = NaN;
    _had_mu_parallel_bnb = NaN; _had_mu_parallel_numi = NaN;
    _had_fwd_penalty_bnb = NaN; _had_fwd_penalty_numi = NaN;
    _hadflow_bnb = NaN; _hadflow_numi = NaN;
}

DEFINE_ART_CLASS_TOOL(VertexTopology)

} // namespace analysis

#endif

