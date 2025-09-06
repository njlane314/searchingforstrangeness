#ifndef ANALYSIS_VERTEXTOPOLOGY_CXX
#define ANALYSIS_VERTEXTOPOLOGY_CXX

#include "AnalysisToolBase.h"
#include "Common/ProxyTypes.h"

#include "canvas/Utilities/InputTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <iostream>

namespace analysis {

class VertexTopology : public AnalysisToolBase {
public:
    explicit VertexTopology(const fhicl::ParameterSet &pset);
    ~VertexTopology() {};
    void configure(fhicl::ParameterSet const &pset) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event,
                      std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
                      bool is_data, bool is_selected) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;
    void print() const;

private:
    art::InputTag fPFPproducer;
    art::InputTag fSPproducer;
    TVector3 fBNBdir;
    TVector3 fNuMIdir;
    float fVtxRadius;
    float fFwdCos;
    float fBackRMin;
    float fBackRMax;
    float fBackMargin;
    float fKernelR;
    int   fThrustIters;
    float fContrastRc;
    float fContrastRa;
    float fContrastRb;

    float _vtx_backfrac_bnb;
    float _vtx_backfrac_numi;
    float _vtx_offfrac_bnb;
    float _vtx_offfrac_numi;

    float _had_thrust_def;
    float _had_sphericity;
    float _had_rho_term;
    float _had_mu_parallel_bnb;
    float _had_mu_parallel_numi;

    float _c_fb_bnb;
    float _c_fb_numi;
    float _c_fwd_bnb;
    float _c_fwd_numi;
    float _c_mu_bnb;
    float _c_mu_numi;
    float _c_thrust;
    float _c_sphericity;
    float _c_rho;
    float _c_vtxdens;

    inline float kernel(float r) const {
        const float x = r / std::max(1e-6f, fKernelR);
        return std::exp(-x * x);
    }

    static inline float clamp01(float x) {
        return (x < 0.f) ? 0.f : (x > 1.f ? 1.f : x);
    }

    void backward_off_axis_fractions(const std::vector<TVector3> &dirs,
                                     const std::vector<float> &weights,
                                     const TVector3 &beam_dir,
                                     float &backward_frac,
                                     float &off_axis_frac);

    float mu_parallel(const std::vector<TVector3> &u,
                      const std::vector<float> &wtilde,
                      const TVector3 &beam_dir) const;

    float thrust_deficit(const std::vector<TVector3> &u,
                         const std::vector<float> &w, int iters) const;

    float sphericity(const std::vector<TVector3> &u,
                     const std::vector<float> &w) const;

    float density_contrast(const std::vector<TVector3> &dirs,
                           const std::vector<float> &wtilde,
                           float Rc, float Ra, float Rb) const;

    float fb_contrast(const std::vector<TVector3> &dirs,
                      const std::vector<float> &weights,
                      const TVector3 &bhat,
                      float rmin, float rmax, float delta) const;

    float rho_contrast(const std::vector<float> &wtilde) const;

    float fwd_contrast(const std::vector<TVector3> &dirs,
                       const std::vector<float> &weights,
                       const TVector3 &bhat,
                       float rvtx, float cf) const;

    float mu_contrast(float mu) const;
};

VertexTopology::VertexTopology(const fhicl::ParameterSet &pset) { configure(pset); }

void VertexTopology::configure(fhicl::ParameterSet const &pset) {
    fPFPproducer = pset.get<art::InputTag>("PFPproducer");
    fSPproducer = pset.get<art::InputTag>("SPproducer", fPFPproducer);
    
    auto bnb = pset.get<std::vector<float>>("BNBBeamDir", {0.f, 0.f, 47000.f});
    if (bnb.size() == 3) fBNBdir = TVector3(bnb[0], bnb[1], bnb[2]).Unit();
    auto numi = pset.get<std::vector<float>>("NuMIBeamDir", {5502.f, 7259.f, 67270.f});
    if (numi.size() == 3) fNuMIdir = TVector3(numi[0], numi[1], numi[2]).Unit();
    
    fVtxRadius = pset.get<float>("VertexRadius", 5.f);
    float theta = pset.get<float>("ForwardAngleDeg", 25.f);
    fFwdCos = std::cos(theta * TMath::DegToRad());
    fBackRMin = pset.get<float>("BackwardRadiusMin", 3.f);
    fBackRMax = pset.get<float>("BackwardRadiusMax", 7.f);
    fBackMargin = pset.get<float>("BackwardCosMargin", 0.05f);
    fKernelR = pset.get<float>("KernelR", 25.f);
    fThrustIters = pset.get<int>("ThrustIters", 6);
    fContrastRc = pset.get<float>("ContrastRc", 5.f);
    fContrastRa = pset.get<float>("ContrastRa", 10.f);
    fContrastRb = pset.get<float>("ContrastRb", 25.f);
}

void VertexTopology::analyseEvent(const art::Event &, bool) { print(); }

void VertexTopology::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool, bool) {
    TVector3 vtx;
    bool has_vtx = false;
    for (auto const &pfp : slice_pfp_vec) {
        int pdg = std::abs(pfp->PdgCode());
        if (pdg == 12 || pdg == 14) {
            auto vtx_v = pfp.get<recob::Vertex>();
            if (vtx_v.size() == 1) {
                double xyz[3];
                vtx_v[0]->XYZ(xyz);
                vtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
                has_vtx = true;
            }
            break;
        }
    }
    if (!has_vtx) {
        print();
        return;
    }

    auto const &pfp_h =
        event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn(pfp_h, event, fSPproducer);

    auto const &sp_h =
        event.getValidHandle<std::vector<recob::SpacePoint>>(fSPproducer);
    art::FindManyP<recob::Hit> sp_hit_assn(sp_h, event, fSPproducer);

    std::vector<TVector3> dirs;
    std::vector<float> weights;
    const size_t nspacepoints = sp_h->size();
    dirs.reserve(nspacepoints);
    weights.reserve(nspacepoints);

    for (auto const &pfp : slice_pfp_vec) {
        int pdg = std::abs(pfp->PdgCode());
        if (pdg == 12 || pdg == 14) continue;
        const std::vector<art::Ptr<recob::SpacePoint>> &sp_v = pfp_spacepoint_assn.at(pfp.index());
        for (auto const &sp : sp_v) {
            double xyz[3];
            sp->XYZ(xyz);
            TVector3 r(xyz[0] - vtx.X(), xyz[1] - vtx.Y(), xyz[2] - vtx.Z());
            if (r.Mag() < 1e-6) continue;
            float weight = 0.f;
            const std::vector<art::Ptr<recob::Hit>> &hits = sp_hit_assn.at(sp.key());
            for (auto const &h : hits) weight += h->Integral();
            dirs.push_back(r);
            weights.push_back(weight);
        }
    }

    backward_off_axis_fractions(dirs, weights, fBNBdir, _vtx_backfrac_bnb,
                                _vtx_offfrac_bnb);
    backward_off_axis_fractions(dirs, weights, fNuMIdir, _vtx_backfrac_numi,
                                _vtx_offfrac_numi);

    std::vector<TVector3> u;
    u.reserve(dirs.size());
    std::vector<float> wtilde;
    wtilde.reserve(dirs.size());
    for (size_t i = 0; i < dirs.size(); ++i) {
        const TVector3 &r = dirs[i];
        const double R = r.Mag();
        const float k = kernel(R);
        u.emplace_back(r.Unit());
        wtilde.emplace_back(k * weights[i]);
    }
    if (wtilde.empty()) {
        print();
        return;
    }

    _had_thrust_def = thrust_deficit(u, wtilde, fThrustIters);
    _had_sphericity = sphericity(u, wtilde);
    _had_rho_term   = rho_contrast(wtilde);
    _had_mu_parallel_bnb  = mu_parallel(u, wtilde, fBNBdir);
    _had_mu_parallel_numi = mu_parallel(u, wtilde, fNuMIdir);

    _c_fb_bnb   = fb_contrast(dirs, weights, fBNBdir, fBackRMin, fBackRMax, fBackMargin);
    _c_fb_numi  = fb_contrast(dirs, weights, fNuMIdir, fBackRMin, fBackRMax, fBackMargin);
    _c_fwd_bnb  = fwd_contrast(dirs, weights, fBNBdir, fVtxRadius, fFwdCos);
    _c_fwd_numi = fwd_contrast(dirs, weights, fNuMIdir, fVtxRadius, fFwdCos);
    _c_mu_bnb   = mu_contrast(_had_mu_parallel_bnb);
    _c_mu_numi  = mu_contrast(_had_mu_parallel_numi);
    _c_thrust   = clamp01(_had_thrust_def);
    _c_sphericity = clamp01(_had_sphericity);
    _c_rho = clamp01(_had_rho_term);
    _c_vtxdens = density_contrast(dirs, wtilde, fContrastRc, fContrastRa, fContrastRb);
    print();
}

void VertexTopology::backward_off_axis_fractions(
    const std::vector<TVector3> &dirs, const std::vector<float> &weights,
    const TVector3 &beam_dir, float &backward_frac,
    float &off_axis_frac) {
    float sum_ann = 0.f, sum_back = 0.f;
    float sum_vtx = 0.f, sum_vtx_off = 0.f;
    for (size_t i = 0; i < dirs.size(); ++i) {
        const TVector3 &r = dirs[i];
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
    backward_frac = (sum_ann > 0) ? sum_back / sum_ann : 0.f;
    off_axis_frac = (sum_vtx > 0) ? sum_vtx_off / sum_vtx : 0.f;
}

float VertexTopology::thrust_deficit(const std::vector<TVector3> &u,
                                     const std::vector<float> &w,
                                     int iters) const {
    double sumW = 0.0;
    TVector3 n0(0, 0, 0);
    for (size_t i = 0; i < u.size(); ++i) {
        sumW += w[i];
        n0 += w[i] * u[i];
    }
    if (sumW <= 0.0) return 0.f;
    TVector3 n = (n0.Mag2() > 0) ? n0.Unit() : TVector3(1, 0, 0);
    for (int k = 0; k < iters; ++k) {
        TVector3 s(0, 0, 0);
        for (size_t i = 0; i < u.size(); ++i) {
            const double sign = (u[i].Dot(n) >= 0.0) ? 1.0 : -1.0;
            s += (w[i] * sign) * u[i];
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

float VertexTopology::sphericity(const std::vector<TVector3> &u,
                                 const std::vector<float> &w) const {
    double sumW = 0.0;
    TMatrixDSym M(3); M.Zero();
    for (size_t i = 0; i < u.size(); ++i) {
        const double wi = w[i];
        sumW += wi;
        const double x = u[i].X(), y = u[i].Y(), z = u[i].Z();
        M(0,0) += wi * x * x;
        M(0,1) += wi * x * y;
        M(0,2) += wi * x * z;
        M(1,1) += wi * y * y;
        M(1,2) += wi * y * z;
        M(2,2) += wi * z * z;
    }
    if (sumW <= 0.0) return 0.f;
    M(1,0) = M(0,1);
    M(2,0) = M(0,2);
    M(2,1) = M(1,2);
    M *= (1.0 / sumW);
    TMatrixDSymEigen eig(M);
    TVectorD vals = eig.GetEigenValues();
    std::array<double,3> lam{vals[0], vals[1], vals[2]};
    std::sort(lam.begin(), lam.end(), std::greater<double>());
    const double S = 1.5 * (lam[1] + lam[2]);
    return static_cast<float>(clamp01(S));
}

float VertexTopology::density_contrast(const std::vector<TVector3> &dirs,
                                       const std::vector<float> &wtilde,
                                       float Rc, float Ra, float Rb) const {
    double sumC=0.0, sumA=0.0;
    const double volC = (4.0/3.0) * TMath::Pi() * Rc*Rc*Rc;
    const double volA = (4.0/3.0) * TMath::Pi() * (Rb*Rb*Rb - Ra*Ra*Ra);
    for (size_t i = 0; i < dirs.size(); ++i) {
        const double R = dirs[i].Mag();
        if (R <= Rc) sumC += wtilde[i];
        else if (R > Ra && R <= Rb) sumA += wtilde[i];
    }
    const double rhoC = (volC>0.0) ? (sumC/volC) : 0.0;
    const double rhoA = (volA>0.0) ? (sumA/volA) : 0.0;
    const double denom = rhoC + rhoA;
    return (denom>0.0) ? static_cast<float>(rhoC / denom) : 0.f;
}

float VertexTopology::fb_contrast(const std::vector<TVector3> &dirs,
                                  const std::vector<float> &weights,
                                  const TVector3 &bhat,
                                  float rmin, float rmax, float delta) const {
    double WF=0.0, WB=0.0;
    for (size_t i=0;i<dirs.size();++i){
        const double R = dirs[i].Mag();
        if (R < rmin || R > rmax) continue;
        const double c = dirs[i].Unit().Dot(bhat);
        if (c >  delta) WF += weights[i];
        if (c < -delta) WB += weights[i];
    }
    const double denom = WF + WB;
    const double afb = (denom>0.0) ? (WF - WB) / denom : 0.0;
    return static_cast<float>(0.5 * (afb + 1.0));
}

float VertexTopology::rho_contrast(const std::vector<float> &wtilde) const {
    double sumW = 0.0;
    for (float wi : wtilde) sumW += wi;
    const double vol = std::pow(TMath::Pi(), 1.5) * std::pow(fKernelR, 3);
    const double rho = (vol > 0.0) ? (sumW / vol) : 0.0;
    const double z   = std::log1p(rho);
    const double y   = (2.0 / TMath::Pi()) * std::atan(z);
    return static_cast<float>(clamp01(y));
}

float VertexTopology::fwd_contrast(const std::vector<TVector3> &dirs,
                                   const std::vector<float> &weights,
                                   const TVector3 &bhat,
                                   float rvtx, float cf) const {
    double Win=0.0, Wout=0.0;
    for (size_t i=0;i<dirs.size();++i){
        const double R = dirs[i].Mag();
        if (R > rvtx) continue;
        const double c = dirs[i].Unit().Dot(bhat);
        if (c >= cf) Win += weights[i];
        else Wout += weights[i];
    }
    const double denom = Win + Wout;
    const double base = (cf<=1.0 && cf>=-1.0) ? (1.0 - cf)/2.0 : 0.0;
    double cval = (denom>0.0) ? (Win/denom) : 0.0;
    if (base < 1.0) cval = (cval - base) / (1.0 - base);
    if (cval < 0.0) cval = 0.0; if (cval > 1.0) cval = 1.0;
    return static_cast<float>(cval);
}

float VertexTopology::mu_contrast(float mu) const {
    float v = (mu - 0.5f) / 0.5f;
    return clamp01(v);
}

float VertexTopology::mu_parallel(const std::vector<TVector3> &u,
                                  const std::vector<float> &wtilde,
                                  const TVector3 &beam_dir) const {
    if (wtilde.empty()) return 0.f;
    double sumW = 0.0, sumWabs = 0.0;
    for (size_t i = 0; i < u.size(); ++i) {
        const double wi = wtilde[i];
        sumW += wi;
        sumWabs += wi * std::abs(u[i].Dot(beam_dir));
    }
    return (sumW > 0.0) ? static_cast<float>(sumWabs / sumW) : 0.f;
}

void VertexTopology::setBranches(TTree *t) {
    t->Branch("vtx_backfrac_bnb",  &_vtx_backfrac_bnb,  "vtx_backfrac_bnb/F");
    t->Branch("vtx_backfrac_numi", &_vtx_backfrac_numi, "vtx_backfrac_numi/F");
    t->Branch("vtx_offfrac_bnb",   &_vtx_offfrac_bnb,   "vtx_offfrac_bnb/F");
    t->Branch("vtx_offfrac_numi",  &_vtx_offfrac_numi,  "vtx_offfrac_numi/F");
    t->Branch("had_thrust_def",    &_had_thrust_def,    "had_thrust_def/F");
    t->Branch("had_sphericity",    &_had_sphericity,    "had_sphericity/F");
    t->Branch("had_rho_term",      &_had_rho_term,      "had_rho_term/F");
    t->Branch("had_mu_parallel_bnb",  &_had_mu_parallel_bnb,  "had_mu_parallel_bnb/F");
    t->Branch("had_mu_parallel_numi", &_had_mu_parallel_numi, "had_mu_parallel_numi/F");

    t->Branch("c_fb_bnb",   &_c_fb_bnb,   "c_fb_bnb/F");
    t->Branch("c_fb_numi",  &_c_fb_numi,  "c_fb_numi/F");
    t->Branch("c_fwd_bnb",  &_c_fwd_bnb,  "c_fwd_bnb/F");
    t->Branch("c_fwd_numi", &_c_fwd_numi, "c_fwd_numi/F");
    t->Branch("c_mu_bnb",   &_c_mu_bnb,   "c_mu_bnb/F");
    t->Branch("c_mu_numi",  &_c_mu_numi,  "c_mu_numi/F");
    t->Branch("c_thrust",   &_c_thrust,   "c_thrust/F");
    t->Branch("c_sphericity", &_c_sphericity, "c_sphericity/F");
    t->Branch("c_rho", &_c_rho, "c_rho/F");
    t->Branch("c_vtxdens", &_c_vtxdens, "c_vtxdens/F");
}

void VertexTopology::resetTTree(TTree *) {
    _vtx_backfrac_bnb = std::numeric_limits<float>::quiet_NaN();
    _vtx_backfrac_numi = std::numeric_limits<float>::quiet_NaN();
    _vtx_offfrac_bnb = std::numeric_limits<float>::quiet_NaN();
    _vtx_offfrac_numi = std::numeric_limits<float>::quiet_NaN();
    _had_thrust_def = std::numeric_limits<float>::quiet_NaN();
    _had_sphericity = std::numeric_limits<float>::quiet_NaN();
    _had_rho_term = std::numeric_limits<float>::quiet_NaN();
    _had_mu_parallel_bnb = std::numeric_limits<float>::quiet_NaN();
    _had_mu_parallel_numi = std::numeric_limits<float>::quiet_NaN();
    _c_fb_bnb = std::numeric_limits<float>::quiet_NaN();
    _c_fb_numi = std::numeric_limits<float>::quiet_NaN();
    _c_fwd_bnb = std::numeric_limits<float>::quiet_NaN();
    _c_fwd_numi = std::numeric_limits<float>::quiet_NaN();
    _c_mu_bnb = std::numeric_limits<float>::quiet_NaN();
    _c_mu_numi = std::numeric_limits<float>::quiet_NaN();
    _c_thrust = std::numeric_limits<float>::quiet_NaN();
    _c_sphericity = std::numeric_limits<float>::quiet_NaN();
    _c_rho = std::numeric_limits<float>::quiet_NaN();
    _c_vtxdens = std::numeric_limits<float>::quiet_NaN();
}

void VertexTopology::print() const {
    std::cout << "Branch values:\n"
              << "  vtx_backfrac_bnb: " << _vtx_backfrac_bnb << '\n'
              << "  vtx_backfrac_numi: " << _vtx_backfrac_numi << '\n'
              << "  vtx_offfrac_bnb: " << _vtx_offfrac_bnb << '\n'
              << "  vtx_offfrac_numi: " << _vtx_offfrac_numi << '\n'
              << "  had_thrust_def: " << _had_thrust_def << '\n'
              << "  had_sphericity: " << _had_sphericity << '\n'
              << "  had_rho_term: " << _had_rho_term << '\n'
              << "  had_mu_parallel_bnb: " << _had_mu_parallel_bnb << '\n'
              << "  had_mu_parallel_numi: " << _had_mu_parallel_numi << '\n'
              << "  c_fb_bnb: " << _c_fb_bnb << '\n'
              << "  c_fb_numi: " << _c_fb_numi << '\n'
              << "  c_fwd_bnb: " << _c_fwd_bnb << '\n'
              << "  c_fwd_numi: " << _c_fwd_numi << '\n'
              << "  c_mu_bnb: " << _c_mu_bnb << '\n'
              << "  c_mu_numi: " << _c_mu_numi << '\n'
              << "  c_thrust: " << _c_thrust << '\n'
              << "  c_sphericity: " << _c_sphericity << '\n'
              << "  c_rho: " << _c_rho << '\n'
              << "  c_vtxdens: " << _c_vtxdens << std::endl;
}

DEFINE_ART_CLASS_TOOL(VertexTopology)

} // namespace analysis

#endif
