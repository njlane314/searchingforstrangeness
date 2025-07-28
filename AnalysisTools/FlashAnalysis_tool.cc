#ifndef ANALYSIS_FLASHANALYSIS_CXX
#define ANALYSIS_FLASHANALYSIS_CXX

#include "AnalysisToolBase.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"

namespace analysis {

class FlashAnalysis : public AnalysisToolBase {
public:
    explicit FlashAnalysis(const fhicl::ParameterSet& pset);
    void configure(const fhicl::ParameterSet& pset) override;
    void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;
    void analyseEvent(const art::Event& event, bool is_data) override {}
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    art::InputTag fSliceProducer;
    art::InputTag fFlashMatchProducer;
    art::InputTag fPFPProducer;

    float _t0;
    float _flash_match_score;
    float _flash_total_pe;
    float _flash_time;
    float _flash_z_center;
    float _flash_z_width;
    float _slice_charge;
    float _slice_z_center;
    float _charge_light_ratio;
    float _flash_slice_z_dist;
    float _flash_pe_per_charge;
};

FlashAnalysis::FlashAnalysis(const fhicl::ParameterSet& pset) {
    this->configure(pset);
}

void FlashAnalysis::configure(const fhicl::ParameterSet& pset) {
    fSliceProducer = pset.get<art::InputTag>("SliceProducer", "pandora");
    fFlashMatchProducer = pset.get<art::InputTag>("FlashMatchProducer", "flashmatch");
    fPFPProducer = pset.get<art::InputTag>("PFPProducer", "pandora");
}

void FlashAnalysis::setBranches(TTree* tree) {
    tree->Branch("t0", &_t0, "t0/F");
    tree->Branch("flash_match_score", &_flash_match_score, "flash_match_score/F");
    tree->Branch("flash_total_pe", &_flash_total_pe, "flash_total_pe/F");
    tree->Branch("flash_time", &_flash_time, "flash_time/F");
    tree->Branch("flash_z_center", &_flash_z_center, "flash_z_center/F");
    tree->Branch("flash_z_width", &_flash_z_width, "flash_z_width/F");
    tree->Branch("slice_charge", &_slice_charge, "slice_charge/F");
    tree->Branch("slice_z_center", &_slice_z_center, "slice_z_center/F");
    tree->Branch("charge_light_ratio", &_charge_light_ratio, "charge_light_ratio/F");
    tree->Branch("flash_slice_z_dist", &_flash_slice_z_dist, "flash_slice_z_dist/F");
    tree->Branch("flash_pe_per_charge", &_flash_pe_per_charge, "flash_pe_per_charge/F");
}

void FlashAnalysis::resetTTree(TTree* tree) {
    _t0 = std::numeric_limits<float>::lowest();
    _flash_match_score = std::numeric_limits<float>::lowest();
    _flash_total_pe = std::numeric_limits<float>::lowest();
    _flash_time = std::numeric_limits<float>::lowest();
    _flash_z_center = std::numeric_limits<float>::lowest();
    _flash_z_width = std::numeric_limits<float>::lowest();
    _slice_charge = std::numeric_limits<float>::lowest();
    _slice_z_center = std::numeric_limits<float>::lowest();
    _charge_light_ratio = std::numeric_limits<float>::lowest();
    _flash_slice_z_dist = std::numeric_limits<float>::lowest();
    _flash_pe_per_charge = std::numeric_limits<float>::lowest();
}

void FlashAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {
    art::Ptr<recob::Slice> current_slice;
    for (const auto& pfp : slice_pfp_v) {
        if (pfp->IsPrimary()) {
            auto const& slice_v = pfp.get<recob::Slice>();
            if (!slice_v.empty()) {
                current_slice = slice_v.front();
                break;
            }
        }
    }
    if (!current_slice) return;

    auto const& slice_h = event.getValidHandle<std::vector<recob::Slice>>(fSliceProducer);
    art::FindManyP<recob::OpFlash> fmFlash(slice_h, event, fFlashMatchProducer);
    art::FindManyP<anab::T0> fmT0(slice_h, event, fFlashMatchProducer);

    const std::vector<art::Ptr<anab::T0>>& t0_v = fmT0.at(current_slice.key());
    const std::vector<art::Ptr<recob::OpFlash>>& flash_v = fmFlash.at(current_slice.key());

    if (!t0_v.empty()) {
        const art::Ptr<anab::T0>& t0 = t0_v.front();
        _t0 = t0->Time();
        _flash_match_score = t0->TriggerConfidence();
    }

    if (!flash_v.empty()) {
        const art::Ptr<recob::OpFlash>& matched_flash = flash_v.front();
        _flash_time = matched_flash->Time();
        _flash_total_pe = matched_flash->TotalPE();
        _flash_z_center = matched_flash->ZCenter();
        _flash_z_width = matched_flash->ZWidth();
    }

    _slice_z_center = current_slice->Center().Z();
    _flash_slice_z_dist = std::abs(_slice_z_center - _flash_z_center);

    art::FindManyP<recob::Hit> fmHits(slice_h, event, fSliceProducer);
    const std::vector<art::Ptr<recob::Hit>>& slice_hits = fmHits.at(current_slice.key());

    _slice_charge = 0.0;
    art::ServiceHandle<geo::Geometry> geom;
    for (const auto& hit : slice_hits) {
        if (geom->View(hit->Channel()) == geo::kZ) {
            _slice_charge += hit->Integral();
        }
    }

    if (_flash_total_pe > 0) {
        _charge_light_ratio = _slice_charge / _flash_total_pe;
    }
    
    if (_slice_charge > 0) {
        _flash_pe_per_charge = _flash_total_pe / _slice_charge;
    }
}

DEFINE_ART_CLASS_TOOL(FlashAnalysis)

}
#endif