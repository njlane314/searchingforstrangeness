#ifndef ANALYSIS_TRACKANALYSIS_CXX
#define ANALYSIS_TRACKANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonFuncs/TypeDefs.h"

#include "../CommonFuncs/BackTrackingFuncs.h"
#include "../CommonFuncs/TrackShowerScoreFuncs.h"
#include "../CommonFuncs/PIDFuncs.h"
#include "../CommonFuncs/SpaceChargeCorrections.h"
#include "../CommonFuncs/Geometry.h"
#include "../CommonFuncs/TypeDefs.h"

namespace analysis
{

class SecondaryVerticesAnalysis : public AnalysisToolBase
{

public:

    SecondaryVerticesAnalysis(const fhicl::ParameterSet &pset);
    ~SecondaryVerticesAnalysis(){};

    void configure(fhicl::ParameterSet const &pset);

    void analyzeEvent(art::Event const &e, bool fData) override;
    void analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

    void SaveTruth(art::Event const &e);
    void fillDefault();

    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:

    art::InputTag fCRTVetoproducer; // producer for CRT veto ass tag [anab::T0 <-> recob::OpFlash]
    art::InputTag fCLSproducer;     // cluster associated to PFP
    art::InputTag fMCTproducer;     // MCTruth from neutrino generator
    art::InputTag fMCPproducer;     // MCParticle from Geant4 stage
    art::InputTag fMCFluxproducer;  // MCFlux producer
    art::InputTag fMCRproducer;
    art::InputTag fSLCproducer; // slice associated to PFP

    art::InputTag fCALOproducer;
    art::InputTag fPIDproducer;
    art::InputTag fTRKproducer;
    art::InputTag fBacktrackTag;
    art::InputTag fHproducer;
    art::InputTag fPFPproducer;
    art::InputTag fVTXproducer;
    art::InputTag fSHRproducer;
    art::InputTag fPCAproducer;

    std::vector<float> _pfp_self_v;
    std::vector<float> _trkshrscore_v;

    std::vector<float> _trk_sep_v;
    std::vector<float> _trk_phi_v;
    std::vector<float> _trk_d_v;
    std::vector<float> _shr_sep_v; 
    std::vector<float> _shr_phi_v;
    std::vector<float> _shr_d_v;
};

SecondaryVerticesAnalysis::SecondaryVerticesAnalysis(const fhicl::ParameterSet &p)
{
    fCALOproducer = p.get<art::InputTag>("CALOproducer");
    fPIDproducer = p.get<art::InputTag>("PIDproducer");
    fTRKproducer = p.get<art::InputTag>("TRKproducer");
    fPFPproducer = p.get<art::InputTag>("PFPproducer");
    fVTXproducer = p.get<art::InputTag>("VTXproducer");
    fBacktrackTag = p.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch");
    fHproducer = p.get<art::InputTag>("Hproducer", "gaushit");
    fSHRproducer = p.get<art::InputTag>("SHRproducer");
    fPCAproducer = p.get<art::InputTag>("PCAproducer");
}

void SecondaryVerticesAnalysis::configure(fhicl::ParameterSet const &p)
{
}

void SecondaryVerticesAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
}

void SecondaryVerticesAnalysis::analyzeSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected)
{
    TVector3 nuvtx;
    for (auto pfp : slice_pfp_v)
    {
        if (pfp->IsPrimary())
        {
            double xyz[3] = {};
            auto vtx = pfp.get<recob::Vertex>();
            if (vtx.size() == 1)
            {
                vtx.at(0)->XYZ(xyz);
                nuvtx.SetXYZ(xyz[0], xyz[1], xyz[2]);
            }

            break;
        }
    }

    /*common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(fCLSproducer),
                                                        proxy::withAssociated<recob::Slice>(fSLCproducer),
                                                        proxy::withAssociated<recob::Track>(fTRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(fVTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(fPCAproducer),
                                                        proxy::withAssociated<recob::Shower>(fSHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(fPFPproducer));*/

    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> assocMCPart;

    if (!fData)
    {
        art::ValidHandle<std::vector<recob::Hit>> inputHits = e.getValidHandle<std::vector<recob::Hit>>(fHproducer);
        assocMCPart = std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(new art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>(inputHits, e, fBacktrackTag));
    }

    float phi_h = std::atan2(nuvtx.Y(), nuvtx.X());

    for (size_t i_pfp = 0; i_pfp < slice_pfp_v.size(); i_pfp++) 
    {
        auto const &pfp_pxy = slice_pfp_v.at(i_pfp);

        auto trk_v = pfp_pxy.get<recob::Track>(); 
        auto shr_v = pfp_pxy.get<recob::Shower>();

        _pfp_self_v.push_back(pfp_pxy->Self());
        _trkshrscore_v.push_back(common::GetTrackShowerScore(pfp_pxy));

        if (trk_v.size() == 1 && shr_v.size() == 1) {
            auto trk = trk_v.at(0);
            auto trk_strt = trk->Start();

            float trk_sep = common::distance3d(nuvtx.X(), nuvtx.Y(), nuvtx.Z(), trk_strt.X(), trk_strt.Y(), trk_strt.Z());
            float trk_phi = std::atan2(trk_strt.X(), trk_strt.Y());
            float trk_d = trk_sep * sin(trk_phi - phi_h);

            _trk_sep_v.push_back(trk_sep); 
            _trk_phi_v.push_back(trk_phi);
            _trk_d_v.push_back(trk_d);

            auto shr = shr_v[0];
            auto shr_strt = shr->ShowerStart();

            float shr_sep = common::distance3d(nuvtx.X(), nuvtx.Y(), nuvtx.Z(), shr_strt.X(), shr_strt.Y(), shr_strt.Z());
            float shr_phi = std::atan2(shr_strt.X(), shr_strt.Y());
            float shr_d = shr_sep * sin(shr_phi - phi_h);

            _shr_sep_v.push_back(shr_sep); 
            _shr_phi_v.push_back(shr_phi);
            _shr_d_v.push_back(shr_d);
        } 
        else {
            std::cout << "Error: track size is " << trk_v.size() << " and shower size is " << shr_v.size() << std::endl;
        }
    }
}

void SecondaryVerticesAnalysis::fillDefault()
{
}

void SecondaryVerticesAnalysis::setBranches(TTree *_tree)
{
    _tree->Branch("secvrtx_pfp_self", "std::vector< float >", &_pfp_self_v);
    _tree->Branch("secvrtx_trkshrscore", "std::vector< float >", &_trkshrscore_v);
    _tree->Branch("secvrtx_trk_sep", "std::vector< float >", &_trk_sep_v);
    _tree->Branch("secvrtx_trk_phi", "std::vector< float >", &_trk_phi_v);
    _tree->Branch("secvrtx_trk_d", "std::vector< float >", &_trk_d_v);
    _tree->Branch("secvrtx_shr_sep", "std::vector< float >", &_shr_sep_v);
    _tree->Branch("secvrtx_shr_phi", "std::vector< float >", &_shr_phi_v);
    _tree->Branch("secvrtx_shr_d", "std::vector< float >", &_shr_d_v);
}

void SecondaryVerticesAnalysis::resetTTree(TTree *_tree)
{
    _pfp_self_v.clear();
    _trkshrscore_v.clear();

    _trk_sep_v.clear();
    _trk_phi_v.clear();
    _trk_d_v.clear();
    _shr_sep_v.clear();
    _shr_phi_v.clear();
    _shr_d_v.clear();
}

DEFINE_ART_CLASS_TOOL(SecondaryVerticesAnalysis)
}

#endif