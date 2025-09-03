#ifndef ANALYSIS_FLASHANALYSIS_CXX
#define ANALYSIS_FLASHANALYSIS_CXX

#include "canvas/Persistency/Common/Ptr.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"

#include "AnalysisToolBase.h"

namespace analysis {

class FlashAnalysis : public AnalysisToolBase {
public:
    explicit FlashAnalysis(const fhicl::ParameterSet& pset);
    void configure(const fhicl::ParameterSet& pset) override;
    void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec, bool is_data, bool is_selected) override;
    void analyseEvent(const art::Event& event, bool is_data) override {}
    void setBranches(TTree* tree) override;
    void resetTTree(TTree* tree) override;

private:
    art::InputTag fPFPproducer;
    art::InputTag fT0producer;
    art::InputTag fOpFlashProducer;

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
    fPFPproducer     = pset.get<art::InputTag>("PFPproducer");
    fT0producer      = pset.get<art::InputTag>("T0producer");
    fOpFlashProducer = pset.get<art::InputTag>("OpFlashProducer");
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

void FlashAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec, bool is_data, bool is_selected) {
    auto const& pfp_h = event.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
    auto const& opflash_h = event.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashProducer);
    auto const& cluster_h = event.getValidHandle<std::vector<recob::Cluster>>(fPFPproducer);

    art::FindManyP<anab::T0> fmT0(pfp_h, event, fT0producer);
    art::FindManyP<recob::Slice> fmSlice(pfp_h, event, fPFPproducer);
    art::FindManyP<recob::Cluster> fmCluster(pfp_h, event, fPFPproducer);
    art::FindManyP<recob::Hit> fmHitFromCluster(cluster_h, event, fPFPproducer);

    art::Ptr<recob::PFParticle> primary_pfp_ptr;
    for (const auto& pfp_pxy : slice_pfp_vec) {
        if (pfp_pxy->IsPrimary()) {
            primary_pfp_ptr = art::Ptr<recob::PFParticle>(pfp_h, pfp_pxy.index());
            break;
        }
    }

    if (!primary_pfp_ptr) return;

    const auto& t0_v = fmT0.at(primary_pfp_ptr.key());
    if (!t0_v.empty()) {
        const art::Ptr<anab::T0>& t0 = t0_v.front();
        _t0 = t0->Time();
        _flash_match_score = t0->TriggerConfidence();

        for (const recob::OpFlash& flash : *opflash_h) {
            if (std::abs(flash.Time() - _t0) < 0.1) {
                _flash_time = flash.Time();
                _flash_total_pe = flash.TotalPE();
                _flash_z_center = flash.ZCenter();
                _flash_z_width = flash.ZWidth();
                break;
            }
        }
    }

    const auto& slice_v = fmSlice.at(primary_pfp_ptr.key());
    if (slice_v.empty()) return;
    const art::Ptr<recob::Slice>& current_slice = slice_v.front();
    
    _slice_z_center = current_slice->Center().Z();
    _flash_slice_z_dist = std::abs(_slice_z_center - _flash_z_center);

    const auto& pfp_clusters = fmCluster.at(primary_pfp_ptr.key());
    _slice_charge = 0.0;
    art::ServiceHandle<geo::Geometry> geom;

    for (const auto& cluster_ptr : pfp_clusters) {
        const auto& cluster_hits = fmHitFromCluster.at(cluster_ptr.key());
        for (const auto& hit_ptr : cluster_hits) {
            if (geom->View(hit_ptr->Channel()) == geo::kZ) {
                _slice_charge += hit_ptr->Integral();
            }
        }
    }
    
    if (_flash_total_pe > 0 && _slice_charge != 0) {
        _charge_light_ratio = _slice_charge / _flash_total_pe;
        _flash_pe_per_charge = _flash_total_pe / _slice_charge;
    }
}

DEFINE_ART_CLASS_TOOL(FlashAnalysis)

}
#endif