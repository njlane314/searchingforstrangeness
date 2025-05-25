#ifndef ANALYSIS_COSMICIP_CXX
#define ANALYSIS_COSMICIP_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "ubobj/CRT/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/FindManyInChainP.h"

#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"

namespace analysis
{
    class NeutrinoEnergy : public AnalysisToolBase {

    public:
        NeutrinoEnergy(const fhicl::ParameterSet& pset);
        ~NeutrinoEnergy(){ };
        void configure(fhicl::ParameterSet const & pset);
        void analyseEvent(art::Event const& e, bool fData) override;
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool fData, bool selected) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;
        
    private:
        art::InputTag fTRKproducer, fCALOproducer, fPDRproducer;
        float fTrkShrScore;
        float fShrEnergyBias;
        float fADCtoMeVMIP;
        float _NeutrinoEnergy0, _NeutrinoEnergy1, _NeutrinoEnergy2;
        float _SliceCaloEnergy0, _SliceCaloEnergy1, _SliceCaloEnergy2;
    };
    
    NeutrinoEnergy::NeutrinoEnergy(const fhicl::ParameterSet& p)
    {
        fPDRproducer = p.get< art::InputTag >("PDRproducer");
        fTRKproducer = p.get< art::InputTag >("TRKproducer");
        fCALOproducer = p.get< art::InputTag >("CALOproducer");
        fTrkShrScore = p.get< float >("TrkShrScore");
        fShrEnergyBias = p.get< float >("ShrEnergyBias");
        fADCtoMeVMIP = p.get< float >("ADCtoMeVMIP");
    }
    
    void NeutrinoEnergy::configure(fhicl::ParameterSet const & p)
    {
    }
    
    void NeutrinoEnergy::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
    {
        common::ProxyCaloColl_t const& calo_proxy = proxy::getCollection<std::vector<recob::Track> >(e, fTRKproducer,
                                                           proxy::withAssociated<anab::Calorimetry>(fCALOproducer));
        common::ProxyClusColl_t const &clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e,fPDRproducer,proxy::withAssociated<recob::Hit>(fPDRproducer));
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
                        _SliceCaloEnergy0 += clus_hit_v[h]->Integral() * fADCtoMeVMIP;
                }
                if (clus->Plane().Plane == 1){
                    for (size_t h=0; h < clus_hit_v.size(); h++)
                        _SliceCaloEnergy1 += clus_hit_v[h]->Integral() * fADCtoMeVMIP;
                }
                if (clus->Plane().Plane == 2){
                    for (size_t h=0; h < clus_hit_v.size(); h++)
                        _SliceCaloEnergy2 += clus_hit_v[h]->Integral() * fADCtoMeVMIP;
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
                    if ( (rr_v.size() == 0) || (dedx.size() == 0) ) continue;
                    size_t nmax = rr_v.size() - 1;
                    if ( (dedx.size() - 1) < nmax) { nmax = dedx.size() - 1; }
                    float dE = 0;
                    for (size_t n=0; n < nmax; n++) {
                        float energy = fabs(rr_v[n]-rr_v[n+1]) * dedx[n];
                        if (energy < 100)
                            dE += energy;
                    }
                    auto const& plane = calo->PlaneID().Plane;
                    if (plane == 0) 
                        _NeutrinoEnergy0 += dE;
                    if (plane == 1) 
                        _NeutrinoEnergy1 += dE;
                    if (plane == 2) 
                        _NeutrinoEnergy2 += dE;
                }
            }
            else {
                auto pxy_shr_v = pfp_pxy.get<recob::Shower>();
                if (pxy_shr_v.size() == 1) {
                    auto const& shr = pxy_shr_v[0];
                    auto const& energy_v = shr->Energy();
                    if (energy_v.size() >= 3) {
                        _NeutrinoEnergy0 += energy_v[0] / fShrEnergyBias;
                        _NeutrinoEnergy1 += energy_v[1] / fShrEnergyBias;
                        _NeutrinoEnergy2 += energy_v[2] / fShrEnergyBias;
                    }
                }
            }
        }
        return;
    }

    void NeutrinoEnergy::analyseEvent(art::Event const &e, bool fData)
    {
    }

    void NeutrinoEnergy::setBranches(TTree* _tree)
    {
        _tree->Branch("NeutrinoEnergy0",&_NeutrinoEnergy0,"NeutrinoEnergy0/F");
        _tree->Branch("NeutrinoEnergy1",&_NeutrinoEnergy1,"NeutrinoEnergy1/F");
        _tree->Branch("NeutrinoEnergy2",&_NeutrinoEnergy2,"NeutrinoEnergy2/F");
        _tree->Branch("SliceCaloEnergy0",&_SliceCaloEnergy0,"SliceCaloEnergy0/F");
        _tree->Branch("SliceCaloEnergy1",&_SliceCaloEnergy1,"SliceCaloEnergy1/F");
        _tree->Branch("SliceCaloEnergy2",&_SliceCaloEnergy2,"SliceCaloEnergy2/F");
    }
    
    void NeutrinoEnergy::resetTTree(TTree* _tree)
    {
        _NeutrinoEnergy0 = 0;
        _NeutrinoEnergy1 = 0;
        _NeutrinoEnergy2 = 0;
        _SliceCaloEnergy0 = 0;
        _SliceCaloEnergy1 = 0;
        _SliceCaloEnergy2 = 0;
    }
    
    DEFINE_ART_CLASS_TOOL(NeutrinoEnergy)
} // namespace analysis

#endif