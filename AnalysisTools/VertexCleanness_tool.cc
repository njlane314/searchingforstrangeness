#ifndef ANALYSIS_VERTEXCLEANNESS_CXX
#define ANALYSIS_VERTEXCLEANNESS_CXX

#include "AnalysisToolBase.h"
#include "Common/ProxyTypes.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "TVector3.h"
#include "TMath.h"

#include <limits>
#include <cmath>
#include <vector>

namespace analysis {

class VertexCleanness : public AnalysisToolBase {
public:
    explicit VertexCleanness(const fhicl::ParameterSet &pset);
    ~VertexCleanness() { };

    void configure(fhicl::ParameterSet const &pset) override;
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event,
                      std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
                      bool is_data, bool is_selected) override;
    void setBranches(TTree *tree) override;
    void resetTTree(TTree *tree) override;

private:
    void compute_scores(const std::vector<TVector3>& dirs,
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
};

VertexCleanness::VertexCleanness(const fhicl::ParameterSet &pset)
{
    configure(pset);
}

void VertexCleanness::configure(fhicl::ParameterSet const &pset)
{
    auto bnb = pset.get<std::vector<float>>("BNBBeamDir", {1,0,0});
    if (bnb.size() == 3)
        fBNBdir = TVector3(bnb[0], bnb[1], bnb[2]).Unit();
    auto numi = pset.get<std::vector<float>>("NuMIBeamDir", {0,0,1});
    if (numi.size() == 3)
        fNuMIdir = TVector3(numi[0], numi[1], numi[2]).Unit();

    fVtxRadius = pset.get<float>("VertexRadius", 5.f);
    float theta = pset.get<float>("ForwardAngleDeg", 25.f);
    fFwdCos = std::cos(theta * TMath::DegToRad());
    fBackRMin = pset.get<float>("BackwardRadiusMin", 3.f);
    fBackRMax = pset.get<float>("BackwardRadiusMax", 7.f);
    fBackMargin = pset.get<float>("BackwardCosMargin", 0.05f);
    fBackQmax = pset.get<float>("BackwardMax", std::numeric_limits<float>::max());
    fVtxQmax = pset.get<float>("VertexResidualMax", std::numeric_limits<float>::max());
    fAlpha = pset.get<float>("Alpha", 1.f);
    fBeta = pset.get<float>("Beta", 1.f);
}

void VertexCleanness::analyseEvent(const art::Event &event, bool is_data) { }

void VertexCleanness::analyseSlice(const art::Event &event,
                                   std::vector<common::ProxyPfpElem_t> &slice_pfp_vec,
                                   bool is_data, bool is_selected)
{
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

    std::vector<TVector3> dirs;
    std::vector<float> weights;
    for (auto const &pfp : slice_pfp_vec) {
        int pdg = std::abs(pfp->PdgCode());
        if (pdg == 12 || pdg == 14)
            continue;
        auto sp_v = pfp.get<recob::SpacePoint>();
        for (auto const &sp : sp_v) {
            double xyz[3];
            sp->XYZ(xyz);
            TVector3 r(xyz[0] - vtx.X(), xyz[1] - vtx.Y(), xyz[2] - vtx.Z());
            if (r.Mag() < 1e-6)
                continue;
            dirs.push_back(r);
            weights.push_back(1.f);
        }
    }

    compute_scores(dirs, weights, fBNBdir,
                   _vtx_backfrac_bnb, _vtx_offfrac_bnb, _vtx_score_bnb);
    compute_scores(dirs, weights, fNuMIdir,
                   _vtx_backfrac_numi, _vtx_offfrac_numi, _vtx_score_numi);
}

void VertexCleanness::compute_scores(const std::vector<TVector3>& dirs,
                                     const std::vector<float>& weights,
                                     const TVector3& beam_dir,
                                     float& back_frac,
                                     float& off_frac,
                                     float& score)
{
    float sum_ann = 0.f;
    float sum_back = 0.f;
    float sum_vtx = 0.f;
    float sum_vtx_off = 0.f;
    for (size_t i = 0; i < dirs.size(); ++i) {
        const auto& r = dirs[i];
        float w = weights[i];
        auto u = r.Unit();
        float proj = u.Dot(beam_dir);
        float mag = r.Mag();

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

    score = fAlpha * back_frac + fBeta * off_frac;
}

void VertexCleanness::setBranches(TTree *tree)
{
    tree->Branch("vtx_score_bnb", &_vtx_score_bnb, "vtx_score_bnb/F");
    tree->Branch("vtx_score_numi", &_vtx_score_numi, "vtx_score_numi/F");
    tree->Branch("vtx_backfrac_bnb", &_vtx_backfrac_bnb, "vtx_backfrac_bnb/F");
    tree->Branch("vtx_backfrac_numi", &_vtx_backfrac_numi, "vtx_backfrac_numi/F");
    tree->Branch("vtx_offfrac_bnb", &_vtx_offfrac_bnb, "vtx_offfrac_bnb/F");
    tree->Branch("vtx_offfrac_numi", &_vtx_offfrac_numi, "vtx_offfrac_numi/F");
}

void VertexCleanness::resetTTree(TTree *tree)
{
    _vtx_score_bnb = std::numeric_limits<float>::quiet_NaN();
    _vtx_score_numi = std::numeric_limits<float>::quiet_NaN();
    _vtx_backfrac_bnb = std::numeric_limits<float>::quiet_NaN();
    _vtx_backfrac_numi = std::numeric_limits<float>::quiet_NaN();
    _vtx_offfrac_bnb = std::numeric_limits<float>::quiet_NaN();
    _vtx_offfrac_numi = std::numeric_limits<float>::quiet_NaN();
}

DEFINE_ART_CLASS_TOOL(VertexCleanness)

} // namespace analysis

#endif
