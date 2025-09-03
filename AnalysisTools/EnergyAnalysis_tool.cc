#ifndef ANALYSIS_ENERGY_CXX
#define ANALYSIS_ENERGY_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "Common/BacktrackingFuncs.h"
#include "Common/TrackShowerScoreFuncs.h"

namespace analysis {

class EnergyAnalysis : public AnalysisToolBase {
public:
    EnergyAnalysis(const fhicl::ParameterSet& pset);
    ~EnergyAnalysis(){ };
    void configure(fhicl::ParameterSet const & pset);
    void analyseEvent(art::Event const& e, bool fData) override;
    void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;
    void setBranches(TTree* _tree) override;
    void resetTTree(TTree* _tree) override;
    
private:
    art::InputTag fTRKproducer;
    art::InputTag fCALproducer;
    art::InputTag fCLSproducer;
    float fTrkShrScore;
    float fShrEnergyBias;
    float fADCtoMeVMIP;

    float _neutrino_energy_0;
    float _neutrino_energy_1; 
    float _neutrino_energy_2;
    float _slice_calo_energy_0;
    float _slice_calo_energy_1;
    float _slice_calo_energy_2;
};

EnergyAnalysis::EnergyAnalysis(const fhicl::ParameterSet& p) {
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fTRKproducer = p.get<art::InputTag>("TRKproducer");
    fCALproducer = p.get<art::InputTag>("CALproducer");
    fTrkShrScore = p.get<float>("TrkShrScore");
    fShrEnergyBias = p.get<float>("ShrEnergyBias");
    fADCtoMeVMIP = p.get<float>("ADCtoMeVMIP");
}

void EnergyAnalysis::configure(fhicl::ParameterSet const & p) {}

void EnergyAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) {
    common::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                                                        proxy::withAssociated<anab::Calorimetry>(fCALproducer));
    common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e,fCLSproducer,proxy::withAssociated<recob::Hit>(fCLSproducer));
    
    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++) {
        auto const &pfp_pxy = slice_pfp_v.at(i_pfp);
        auto PDG = fabs(pfp_pxy->PdgCode());
        if ((PDG == 12) || (PDG == 14))
            continue;

        auto trkshrscore = common::GetTrackShowerScore(pfp_pxy);
        auto pxy_cls_v = pfp_pxy.get<recob::Cluster>();
        for (size_t c=0; c < pxy_cls_v.size(); c++) {
            const auto &clus = clus_proxy[(pxy_cls_v[c]).key()];
            auto clus_hit_v = clus.get<recob::Hit>();
            if (clus->Plane().Plane == 0){
                for (size_t h=0; h < clus_hit_v.size(); h++)
                    _slice_calo_energy_0 += clus_hit_v[h]->Integral() * fADCtoMeVMIP;
            }
            if (clus->Plane().Plane == 1){
                for (size_t h=0; h < clus_hit_v.size(); h++)
                    _slice_calo_energy_1 += clus_hit_v[h]->Integral() * fADCtoMeVMIP;
            }
            if (clus->Plane().Plane == 2){
                for (size_t h=0; h < clus_hit_v.size(); h++)
                    _slice_calo_energy_2 += clus_hit_v[h]->Integral() * fADCtoMeVMIP;
            }
        }
        if (trkshrscore > fTrkShrScore) {
            auto pxy_trk_v = pfp_pxy.get<recob::Track>();
            if (pxy_trk_v.size() != 1) continue;
            auto calo_v = calo_proxy[(pxy_trk_v[0]).key()].get<anab::Calorimetry>();
            if (calo_v.size() != 1) continue;
            for (auto const& calo : calo_v) {
                auto rr_v = calo->ResidualRange();
                auto dedx = calo->dEdx();
                if ((rr_v.size() == 0) || (dedx.size() == 0)) continue;
                size_t nmax = rr_v.size() - 1;
                if ((dedx.size() - 1) < nmax){ nmax = dedx.size() - 1; }
                float dE = 0;
                for (size_t n=0; n < nmax; n++) {
                    float energy = fabs(rr_v[n] - rr_v[n+1]) * dedx[n];
                    if (energy < 100)
                        dE += energy;
                }
                auto const& plane = calo->PlaneID().Plane;
                if (plane == 0) 
                    _neutrino_energy_0 += dE;
                if (plane == 1) 
                    _neutrino_energy_1 += dE;
                if (plane == 2) 
                    _neutrino_energy_2 += dE;
            }
        }
        else {
            auto pxy_shr_v = pfp_pxy.get<recob::Shower>();
            if (pxy_shr_v.size() == 1) {
                auto const& shr = pxy_shr_v[0];
                auto const& energy_v = shr->Energy();
                if (energy_v.size() >= 3) {
                    _neutrino_energy_0 += energy_v[0] / fShrEnergyBias;
                    _neutrino_energy_1 += energy_v[1] / fShrEnergyBias;
                    _neutrino_energy_2 += energy_v[2] / fShrEnergyBias;
                }
            }
        }
    }
    return;
}

void EnergyAnalysis::analyseEvent(art::Event const &e, bool fData) {}

void EnergyAnalysis::setBranches(TTree* _tree) {
    _tree->Branch("neutrino_energy_0",&_neutrino_energy_0,"neutrino_energy_0/F");
    _tree->Branch("neutrino_energy_1",&_neutrino_energy_1,"neutrino_energy_1/F");
    _tree->Branch("neutrino_energy_2",&_neutrino_energy_2,"neutrino_energy_2/F");
    _tree->Branch("slice_calo_energy_0",&_slice_calo_energy_0,"slice_calo_energy_0/F");
    _tree->Branch("slice_calo_energy_1",&_slice_calo_energy_1,"slice_calo_energy_1/F");
    _tree->Branch("slice_calo_energy_2",&_slice_calo_energy_2,"slice_calo_energy_2/F");
}

void EnergyAnalysis::resetTTree(TTree* _tree) {
    _neutrino_energy_0 = std::numeric_limits<float>::lowest();
    _neutrino_energy_1 = std::numeric_limits<float>::lowest();
    _neutrino_energy_2 = std::numeric_limits<float>::lowest();
    _slice_calo_energy_0 = std::numeric_limits<float>::lowest();
    _slice_calo_energy_1 = std::numeric_limits<float>::lowest();
    _slice_calo_energy_2 = std::numeric_limits<float>::lowest();
}

DEFINE_ART_CLASS_TOOL(EnergyAnalysis)

}

#endif
