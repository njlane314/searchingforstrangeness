#ifndef ANALYSIS_ENERGYANALYSIS_CXX
#define ANALYSIS_ENERGYANALYSIS_CXX

#include "AnalysisToolBase.h" 
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/Types.h"

namespace analysis {

class EnergyAnalysis : public AnalysisToolBase {
public:
    explicit EnergyAnalysis(const fhicl::ParameterSet& pset);
    void configure(fhicl::ParameterSet const & pset) override;
    void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;
    void analyseEvent(art::Event const& e, bool fData) override {}
    void setBranches(TTree* _tree) override;
    void resetTTree(TTree* _tree) override;

private:
    art::InputTag fTrackProducer;
    art::InputTag fCalorimetryProducer;
    art::InputTag fClusterProducer;
    float fTrkShrScore;
    float fShrEnergyBias;
    float fADCtoMeVMIP;

    float _neutrino_energy_plane0;
    float _neutrino_energy_plane1;
    float _neutrino_energy_plane2;
    float _slice_calo_energy_plane0;
    float _slice_calo_energy_plane1;
    float _slice_calo_energy_plane2;
};

EnergyAnalysis::EnergyAnalysis(const fhicl::ParameterSet& pset) {
    this->configure(pset);
}

void EnergyAnalysis::configure(fhicl::ParameterSet const & pset) {
    fTrackProducer       = pset.get<art::InputTag>("TrackProducer");
    fCalorimetryProducer = pset.get<art::InputTag>("CalorimetryProducer");
    fClusterProducer     = pset.get<art::InputTag>("ClusterProducer");
    fTrkShrScore         = pset.get<float>("TrkShrScore");
    fShrEnergyBias       = pset.get<float>("ShrEnergyBias");
    fADCtoMeVMIP         = pset.get<float>("ADCtoMeVMIP");
}

void EnergyAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) {
    auto const& calo_proxy = proxy::getCollection<std::vector<recob::Track>>(e, fTrackProducer,
        proxy::withAssociated<anab::Calorimetry>(fCalorimetryProducer));
    auto const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fClusterProducer,
        proxy::withAssociated<recob::Hit>(fClusterProducer));

    for (const auto& pfp_pxy : slice_pfp_v) {
        auto pdg = std::abs(pfp_pxy->PdgCode());
        if (pdg == 12 || pdg == 14) continue;

        auto trkshrscore = common::GetTrackShowerScore(pfp_pxy);
        auto pxy_cls_v = pfp_pxy.get<recob::Cluster>();

        for (const auto& pxy_cls : pxy_cls_v) {
            const auto& clus = clus_proxy[pxy_cls.key()];
            auto clus_hit_v = clus.get<recob::Hit>();
            float& slice_calo_energy = (clus->Plane().Plane == 0) ? _slice_calo_energy_plane0 :
                                       (clus->Plane().Plane == 1) ? _slice_calo_energy_plane1 : _slice_calo_energy_plane2;
            for (const auto& hit : clus_hit_v) {
                slice_calo_energy += hit->Integral() * fADCtoMeVMIP;
            }
        }

        if (trkshrscore > fTrkShrScore) {
            auto pxy_trk_v = pfp_pxy.get<recob::Track>();
            if (pxy_trk_v.size() != 1) continue;

            auto calo_v = calo_proxy[pxy_trk_v[0].key()].get<anab::Calorimetry>();
            for (const auto& calo : calo_v) {
                auto rr_v = calo->ResidualRange();
                auto dedx = calo->dEdx();
                if (rr_v.empty() || dedx.empty()) continue;

                size_t nmax = std::min(rr_v.size(), dedx.size()) - 1;
                float dE = 0.0f;
                for (size_t n = 0; n < nmax; ++n) {
                    float energy_step = std::abs(rr_v[n] - rr_v[n + 1]) * dedx[n];
                    if (energy_step < 100.0) dE += energy_step;
                }

                float& neutrino_energy = (calo->PlaneID().Plane == 0) ? _neutrino_energy_plane0 :
                                         (calo->PlaneID().Plane == 1) ? _neutrino_energy_plane1 : _neutrino_energy_plane2;
                neutrino_energy += dE;
            }
        } else {
            auto pxy_shr_v = pfp_pxy.get<recob::Shower>();
            if (pxy_shr_v.size() > 0) { 
                _neutrino_energy_plane0 += pxy_shr_v[0]->Energy()[0] / fShrEnergyBias;
                _neutrino_energy_plane1 += pxy_shr_v[0]->Energy()[1] / fShrEnergyBias;
                _neutrino_energy_plane2 += pxy_shr_v[0]->Energy()[2] / fShrEnergyBias;
            }
        }
    }
}

void EnergyAnalysis::setBranches(TTree* _tree) {
    _tree->Branch("neutrino_energy_plane0", &_neutrino_energy_plane0, "neutrino_energy_plane0/F");
    _tree->Branch("neutrino_energy_plane1", &_neutrino_energy_plane1, "neutrino_energy_plane1/F");
    _tree->Branch("neutrino_energy_plane2", &_neutrino_energy_plane2, "neutrino_energy_plane2/F");
    _tree->Branch("slice_calo_energy_plane0", &_slice_calo_energy_plane0, "slice_calo_energy_plane0/F");
    _tree->Branch("slice_calo_energy_plane1", &_slice_calo_energy_plane1, "slice_calo_energy_plane1/F");
    _tree->Branch("slice_calo_energy_plane2", &_slice_calo_energy_plane2, "slice_calo_energy_plane2/F");
}

void EnergyAnalysis::resetTTree(TTree* _tree) {
    _neutrino_energy_plane0 = 0.0;
    _neutrino_energy_plane1 = 0.0;
    _neutrino_energy_plane2 = 0.0;
    _slice_calo_energy_plane0 = 0.0;
    _slice_calo_energy_plane1 = 0.0;
    _slice_calo_energy_plane2 = 0.0;
}

DEFINE_ART_CLASS_TOOL(EnergyAnalysis)
}
#endif