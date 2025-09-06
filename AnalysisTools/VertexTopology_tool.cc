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
#include <cmath>
#include <limits>
#include <vector>

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

private:
  art::InputTag fPFPproducer;
  art::InputTag fSpacePointproducer;
  // ---------- Original "VertexCleanness" metrics ----------
  void compute_back_off_fractions(const std::vector<TVector3> &dirs,
                                  const std::vector<float> &weights,
                                  const TVector3 &beam_dir, float &back_frac,
                                  float &off_frac);

  TVector3 fBNBdir;
  TVector3 fNuMIdir;

  float fVtxRadius;
  float fFwdCos;
  float fBackRMin;
  float fBackRMax;
  float fBackMargin;
  float fBackQmax;
  float fVtxQmax;

  // Composite vertex score parameters removed; only primitives retained
  float _vtx_backfrac_bnb;
  float _vtx_backfrac_numi;
  float _vtx_offfrac_bnb;
  float _vtx_offfrac_numi;

  // ---------- HadFlow+ event-shape around the vertex ----------
  // Config
  float fKernelR; // cm, soft radial kernel scale
  float fRhoRef;  // reference density for normalization
  // Composite score parameters removed; downstream tuning expected
  int fThrustIters; // iterations for thrust maximization

  // Outputs (beam-independent primitives)
  float _had_thrust_def; // 1 - thrust
  float _had_sphericity; // sphericity
  float _had_rho_term;   // normalized local activity

  // Beam-dependent primitives
  float _had_mu_parallel_bnb;
  float _had_mu_parallel_numi;

private:
  // HadFlow core (only primitive, no composite score)
  float compute_mu_parallel(const std::vector<TVector3> &u,
                            const std::vector<float> &wtilde,
                            const TVector3 &beam_dir) const;

  // Helpers
  inline float kernel(float r) const {
    const float x = r / std::max(1e-6f, fKernelR);
    return std::exp(-x * x);
  }
  static inline float clamp01(float x) {
    return (x < 0.f) ? 0.f : (x > 1.f ? 1.f : x);
  }

  float thrust_deficit(const std::vector<TVector3> &u,
                       const std::vector<float> &w, int iters) const;

  float sphericity(const std::vector<TVector3> &u,
                   const std::vector<float> &w) const;

  float rho_term(const std::vector<float> &w) const;
};

// ===== Implementation =====

VertexTopology::VertexTopology(const fhicl::ParameterSet &pset) {
  configure(pset);
}

void VertexTopology::configure(fhicl::ParameterSet const &pset) {
  fPFPproducer = pset.get<art::InputTag>("PFPproducer");
  fSpacePointproducer =
      pset.get<art::InputTag>("SpacePointproducer", fPFPproducer);

  // Beam directions can be provided via FHiCL; otherwise use fixed defaults
  // Nominal BNB beam direction from target to detector
  auto bnb = pset.get<std::vector<float>>("BNBBeamDir", {0.f, 0.f, 47000.f});
  if (bnb.size() == 3)
    fBNBdir = TVector3(bnb[0], bnb[1], bnb[2]).Unit();

  // Nominal NuMI beam direction from target to detector
  auto numi =
      pset.get<std::vector<float>>("NuMIBeamDir", {5502.f, 7259.f, 67270.f});
  if (numi.size() == 3)
    fNuMIdir = TVector3(numi[0], numi[1], numi[2]).Unit();

  // --- Original cleanness defaults ---
  fVtxRadius = pset.get<float>("VertexRadius", 5.f);
  float theta = pset.get<float>("ForwardAngleDeg", 25.f);
  fFwdCos = std::cos(theta * TMath::DegToRad());
  fBackRMin = pset.get<float>("BackwardRadiusMin", 3.f);
  fBackRMax = pset.get<float>("BackwardRadiusMax", 7.f);
  fBackMargin = pset.get<float>("BackwardCosMargin", 0.05f);
  fBackQmax = pset.get<float>("BackwardMax", std::numeric_limits<float>::max());
  fVtxQmax =
      pset.get<float>("VertexResidualMax", std::numeric_limits<float>::max());

  // --- HadFlow+ defaults (primitives only) ---
  fKernelR = pset.get<float>("KernelR", 25.f); // cm
  fRhoRef = pset.get<float>("RhoRef", 50.f);   // tune on sidebands
  fThrustIters = pset.get<int>("ThrustIters", 6);
}

void VertexTopology::analyseEvent(const art::Event &, bool) {}

void VertexTopology::analyseSlice(
    const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
    bool /*is_data*/, bool) {
  // Beam directions are supplied by configuration; no per-event recalculation

  // --- Locate neutrino vertex (same as your original) ---
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
  if (!has_vtx)
    return;

  // --- Get PFParticle -> SpacePoint associations and SpacePoint -> Hit ---
  auto const &pfp_h =
      event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
  art::FindManyP<recob::SpacePoint> pfp_spacepoint_assn(pfp_h, event,
                                                        fSpacePointproducer);

  auto const &sp_h =
      event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointproducer);
  art::FindManyP<recob::Hit> sp_hit_assn(sp_h, event, fSpacePointproducer);

  // --- Gather displacement vectors and weights (sum charge of associated hits)
  // ---
  std::vector<TVector3> dirs;
  std::vector<float> weights;
  // Reserve enough space for all space points in the event
  size_t nspacepoints = sp_h->size();
  dirs.reserve(nspacepoints);
  weights.reserve(nspacepoints);

  for (auto const &pfp : slice_pfp_vec) {
    int pdg = std::abs(pfp->PdgCode());
    if (pdg == 12 || pdg == 14)
      continue; // skip neutrino PFP
    const std::vector<art::Ptr<recob::SpacePoint>> &sp_v =
        pfp_spacepoint_assn.at(pfp.index());
    for (auto const &sp : sp_v) {
      double xyz[3];
      sp->XYZ(xyz);
      TVector3 r(xyz[0] - vtx.X(), xyz[1] - vtx.Y(), xyz[2] - vtx.Z());
      if (r.Mag() < 1e-6)
        continue;

      float weight = 1.f;
      const std::vector<art::Ptr<recob::Hit>> &hits = sp_hit_assn.at(sp.key());
      if (!hits.empty()) {
        weight = 0.f;
        for (auto const &h : hits)
          weight += h->Integral();
      }

      dirs.push_back(r);
      weights.push_back(weight);
    }
  }

  // --- Original vertex-cleanness summaries (BNB / NuMI) ---
  compute_back_off_fractions(dirs, weights, fBNBdir, _vtx_backfrac_bnb,
                             _vtx_offfrac_bnb);
  compute_back_off_fractions(dirs, weights, fNuMIdir, _vtx_backfrac_numi,
                             _vtx_offfrac_numi);

  // --- HadFlow+ primitives (make unit directions + kernel-weighted weights)
  // ---
  std::vector<TVector3> u;
  u.reserve(dirs.size());
  std::vector<float> wtilde;
  wtilde.reserve(dirs.size());

  for (size_t i = 0; i < dirs.size(); ++i) {
    const TVector3 &r = dirs[i];
    const double R = r.Mag();
    const float k = kernel(R);
    u.emplace_back(r.Unit());
    // multiply by any upstream weight (currently 1.f)
    wtilde.emplace_back(k * weights[i]);
  }

  if (wtilde.empty()) {
    _had_thrust_def = 0.f;
    _had_sphericity = 0.f;
    _had_rho_term = 0.f;

    _had_mu_parallel_bnb = 0.f;
    _had_mu_parallel_numi = 0.f;
    return;
  }

  _had_thrust_def = thrust_deficit(u, wtilde, fThrustIters);
  _had_sphericity = sphericity(u, wtilde);
  _had_rho_term = rho_term(wtilde);

  // Beam-aware HadFlow primitive
  _had_mu_parallel_bnb = compute_mu_parallel(u, wtilde, fBNBdir);
  _had_mu_parallel_numi = compute_mu_parallel(u, wtilde, fNuMIdir);
}

// ---------- Original vertex-cleanness computation ----------
void VertexTopology::compute_back_off_fractions(
    const std::vector<TVector3> &dirs, const std::vector<float> &weights,
    const TVector3 &beam_dir, float &back_frac, float &off_frac) {
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
      if (proj < -fBackMargin)
        sum_back += w;
    }

    if (mag <= fVtxRadius) {
      sum_vtx += w;
      if (proj < fFwdCos)
        sum_vtx_off += w;
    }
  }

  back_frac = (sum_ann > 0) ? sum_back / sum_ann : 0.f;
  if (sum_back > fBackQmax)
    back_frac = 1.f;

  off_frac = (sum_vtx > 0) ? sum_vtx_off / sum_vtx : 0.f;
  if (sum_vtx_off > fVtxQmax)
    off_frac = 1.f;
}

// ---------- HadFlow helpers ----------
float VertexTopology::thrust_deficit(const std::vector<TVector3> &u,
                                     const std::vector<float> &w,
                                     int iters) const {
  double sumW = 0.0;
  TVector3 n0(0, 0, 0);
  for (size_t i = 0; i < u.size(); ++i) {
    sumW += w[i];
    n0 += w[i] * u[i];
  }
  if (sumW <= 0.0)
    return 0.f;
  TVector3 n = (n0.Mag2() > 0) ? n0.Unit() : TVector3(1, 0, 0);

  for (int k = 0; k < iters; ++k) {
    TVector3 s(0, 0, 0);
    for (size_t i = 0; i < u.size(); ++i) {
      const double sign = (u[i].Dot(n) >= 0.0) ? 1.0 : -1.0;
      s += (w[i] * sign) * u[i];
    }
    if (s.Mag2() < 1e-12)
      break;
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
  TMatrixDSym M(3);
  M.Zero();
  for (size_t i = 0; i < u.size(); ++i) {
    const double wi = w[i];
    sumW += wi;
    const double x = u[i].X(), y = u[i].Y(), z = u[i].Z();
    M(0, 0) += wi * x * x;
    M(0, 1) += wi * x * y;
    M(0, 2) += wi * x * z;
    M(1, 1) += wi * y * y;
    M(1, 2) += wi * y * z;
    M(2, 2) += wi * z * z;
  }
  if (sumW <= 0.0)
    return 0.f;
  M(1, 0) = M(0, 1);
  M(2, 0) = M(0, 2);
  M(2, 1) = M(1, 2);
  M *= (1.0 / sumW);

  TMatrixDSymEigen eig(M);
  TVectorD vals = eig.GetEigenValues(); // ascending
  std::array<double, 3> lam{vals[0], vals[1], vals[2]};
  std::sort(lam.begin(), lam.end(), std::greater<double>());

  const double S = 1.5 * (lam[1] + lam[2]); // [0,1]
  return static_cast<float>(clamp01(S));
}

float VertexTopology::rho_term(const std::vector<float> &w) const {
  double sumW = 0.0;
  for (float wi : w)
    sumW += wi;
  const double vol = (4.0 / 3.0) * TMath::Pi() *
                     std::pow(fKernelR, 3); // effective ball volume at scale R
  const double rho = (vol > 0.0) ? (sumW / vol) : 0.0;
  const double term = (fRhoRef > 0.0) ? std::min(rho / fRhoRef, 1.0) : 0.0;
  return static_cast<float>(clamp01(term));
}

float VertexTopology::compute_mu_parallel(const std::vector<TVector3> &u,
                                          const std::vector<float> &wtilde,
                                          const TVector3 &beam_dir) const {
  if (wtilde.empty())
    return 0.f;
  double sumW = 0.0, sumWabs = 0.0;
  for (size_t i = 0; i < u.size(); ++i) {
    const double wi = wtilde[i];
    sumW += wi;
    sumWabs += wi * std::abs(u[i].Dot(beam_dir));
  }
  return (sumW > 0.0) ? static_cast<float>(sumWabs / sumW) : 0.f;
}

// ---------- I/O ----------
void VertexTopology::setBranches(TTree *t) {
  // Original cleanness (kept for compatibility)
  t->Branch("vtx_backfrac_bnb", &_vtx_backfrac_bnb, "vtx_backfrac_bnb/F");
  t->Branch("vtx_backfrac_numi", &_vtx_backfrac_numi, "vtx_backfrac_numi/F");
  t->Branch("vtx_offfrac_bnb", &_vtx_offfrac_bnb, "vtx_offfrac_bnb/F");
  t->Branch("vtx_offfrac_numi", &_vtx_offfrac_numi, "vtx_offfrac_numi/F");

  // HadFlow+ primitives
  t->Branch("had_thrust_def", &_had_thrust_def, "had_thrust_def/F");
  t->Branch("had_sphericity", &_had_sphericity, "had_sphericity/F");
  t->Branch("had_rho_term", &_had_rho_term, "had_rho_term/F");

  t->Branch("had_mu_parallel_bnb", &_had_mu_parallel_bnb,
            "had_mu_parallel_bnb/F");
  t->Branch("had_mu_parallel_numi", &_had_mu_parallel_numi,
            "had_mu_parallel_numi/F");
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
}

DEFINE_ART_CLASS_TOOL(VertexTopology)

} // namespace analysis

#endif
