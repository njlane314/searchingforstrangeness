#ifndef ANALYSIS_BLIPANALYSIS_CXX
#define ANALYSIS_BLIPANALYSIS_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "../CommonDefs/Types.h"

#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/TrackShowerScoreFuncs.h"
#include "../CommonDefs/PIDFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Geometry.h"

#include "../CommonDefs/LLR_PID.h"
#include "../CommonDefs/LLRPID_proton_muon_lookup.h"
#include "../CommonDefs/LLRPID_correction_lookup.h"
#include "../CommonDefs/CalibrationFuncs.h"

#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"

namespace analysis {

class BlipAnalysis : public AnalysisToolBase {
public:
    BlipAnalysis(const fhicl::ParameterSet &pset);
    ~BlipAnalysis(){};
    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(art::Event const &e, bool fData) override;
    void analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;
    void SaveTruth(art::Event const &e);
    void fillDefault();
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    int _run, _sub, _evt;
    blip::BlipRecoAlg* fBlipAlg;
    std::vector<int>   _blip_ID;
    std::vector<bool>  _blip_isValid;
    std::vector<int>   _blip_TPC;
    std::vector<int>   _blip_NPlanes;
    std::vector<int>   _blip_MaxWireSpan;
    std::vector<float> _blip_Energy;
    std::vector<float> _blip_EnergyESTAR;
    std::vector<float> _blip_Time;
    std::vector<float> _blip_ProxTrkDist;
    std::vector<int>   _blip_ProxTrkID;
    std::vector<bool>  _blip_inCylinder;
    std::vector<float> _blip_X;
    std::vector<float> _blip_Y;
    std::vector<float> _blip_Z;
    std::vector<float> _blip_SigmaYZ;
    std::vector<float> _blip_dX;
    std::vector<float> _blip_dYZ;
    std::vector<float> _blip_Charge;
    std::vector<int>   _blip_LeadG4ID;
    std::vector<int>         _blip_pdg;
    std::vector<std::string> _blip_process;
    std::vector<float>       _blip_vx;
    std::vector<float>       _blip_vy;
    std::vector<float>       _blip_vz;
    std::vector<float>       _blip_E;
    std::vector<float>       _blip_mass;
    std::vector<int>         _blip_trkid;
};

BlipAnalysis::BlipAnalysis(const fhicl::ParameterSet &pset) {
    fBlipAlg = new blip::BlipRecoAlg(pset.get<fhicl::ParameterSet>("BlipAlg"));
}

void BlipAnalysis::configure(fhicl::ParameterSet const &pset) {
}

void BlipAnalysis::analyseEvent(art::Event const &e, bool fData) {
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();
    std::cout << "[BlipAnalysis::analyseEvent] Run: " << _run << ", SubRun: " << _sub << ", Event: " << _evt << std::endl;
}

void BlipAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) {
    fBlipAlg->RunBlipReco(e);
    std::vector<blip::Blip> blipVec = fBlipAlg->blips;
    std::cout << "number of blips = " << blipVec.size() << std::endl;
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    art::Handle< std::vector<simb::MCParticle> > mcparticleHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    std::vector<art::Ptr<simb::MCParticle> > mcparticle;
    if (e.getByLabel("generator", mctruthListHandle)) {
        art::fill_ptr_vector(mclist, mctruthListHandle);
    }
    if (e.getByLabel("largeant", mcparticleHandle)) {
        art::fill_ptr_vector(mcparticle, mcparticleHandle);
    }
    for (size_t i = 0; i < blipVec.size(); i++) {
        _blip_ID.push_back(blipVec[i].ID);
        _blip_isValid.push_back(blipVec[i].isValid);
        _blip_TPC.push_back(blipVec[i].TPC);
        _blip_NPlanes.push_back(blipVec[i].NPlanes);
        _blip_MaxWireSpan.push_back(blipVec[i].MaxWireSpan);
        _blip_Energy.push_back(blipVec[i].Energy);
        _blip_EnergyESTAR.push_back(blipVec[i].EnergyESTAR);
        _blip_Time.push_back(blipVec[i].Time);
        _blip_ProxTrkDist.push_back(blipVec[i].ProxTrkDist);
        _blip_ProxTrkID.push_back(blipVec[i].ProxTrkID);
        _blip_inCylinder.push_back(blipVec[i].inCylinder);
        _blip_X.push_back(blipVec[i].X());
        _blip_Y.push_back(blipVec[i].Y());
        _blip_Z.push_back(blipVec[i].Z());
        _blip_SigmaYZ.push_back(blipVec[i].SigmaYZ);
        _blip_dX.push_back(blipVec[i].dX);
        _blip_dYZ.push_back(blipVec[i].dYZ);
        _blip_Charge.push_back(blipVec[i].Charge);
        _blip_LeadG4ID.push_back(blipVec[i].truth.LeadG4ID);
        int b_pdg = -999;
        std::string b_pro = "null";
        float b_vx = -999;
        float b_vy = -999;
        float b_vz = -999;
        float b_E = -999;
        float b_mass = -999;
        int b_tid = -999;
        if (!fData) {
            b_pdg = blipVec[i].truth.LeadG4PDG;
            int blip_to_mcpar = -1;
            std::vector<int> _neutron_trkid;
            for (int i_mcp = 0; i_mcp < (int)mcparticle.size(); i_mcp++) {
                if (mcparticle[i_mcp]->TrackId() == blipVec[i].truth.LeadG4ID) {
                    blip_to_mcpar = i_mcp;
                }
            }
            if (blip_to_mcpar > -1) {
                for (int i_n = 0; i_n < (int)_neutron_trkid.size(); i_n++) {
                }
                b_pro = mcparticle[blip_to_mcpar]->Process();
                b_vx = mcparticle[blip_to_mcpar]->Vx();
                b_vy = mcparticle[blip_to_mcpar]->Vy();
                b_vz = mcparticle[blip_to_mcpar]->Vz();
                b_mass = mcparticle[blip_to_mcpar]->Mass();
                b_E = mcparticle[blip_to_mcpar]->E();
                b_tid = mcparticle[blip_to_mcpar]->TrackId();
            }
        }
        _blip_pdg.push_back(b_pdg);
        _blip_process.push_back(b_pro);
        _blip_vx.push_back(b_vx);
        _blip_vy.push_back(b_vy);
        _blip_vz.push_back(b_vz);
        _blip_E.push_back(b_E);
        _blip_mass.push_back(b_mass);
        _blip_trkid.push_back(b_tid);
    }
    std::cout << "finished analysing blip slice" << std::endl;
}

void BlipAnalysis::fillDefault() {
    _blip_ID.push_back(std::numeric_limits<int>::lowest());
    _blip_isValid.push_back(false);
    _blip_TPC.push_back(std::numeric_limits<int>::lowest());
    _blip_NPlanes.push_back(std::numeric_limits<int>::lowest());
    _blip_MaxWireSpan.push_back(std::numeric_limits<int>::lowest());
    _blip_Energy.push_back(std::numeric_limits<float>::lowest());
    _blip_EnergyESTAR.push_back(std::numeric_limits<float>::lowest());
    _blip_Time.push_back(std::numeric_limits<float>::lowest());
    _blip_ProxTrkDist.push_back(std::numeric_limits<float>::lowest());
    _blip_ProxTrkID.push_back(std::numeric_limits<int>::lowest());
    _blip_inCylinder.push_back(false);
    _blip_X.push_back(std::numeric_limits<float>::lowest());
    _blip_Y.push_back(std::numeric_limits<float>::lowest());
    _blip_Z.push_back(std::numeric_limits<float>::lowest());
    _blip_SigmaYZ.push_back(std::numeric_limits<float>::lowest());
    _blip_dX.push_back(std::numeric_limits<float>::lowest());
    _blip_dYZ.push_back(std::numeric_limits<float>::lowest());
    _blip_Charge.push_back(std::numeric_limits<float>::lowest());
    _blip_pdg.push_back(std::numeric_limits<int>::lowest());
    _blip_process.push_back("null");
    _blip_vx.push_back(std::numeric_limits<float>::lowest());
    _blip_vy.push_back(std::numeric_limits<float>::lowest());
    _blip_vz.push_back(std::numeric_limits<float>::lowest());
    _blip_E.push_back(std::numeric_limits<float>::lowest());
    _blip_mass.push_back(std::numeric_limits<float>::lowest());
    _blip_trkid.push_back(std::numeric_limits<int>::lowest());
}

void BlipAnalysis::setBranches(TTree *_tree) {
    _tree->Branch("blip_ID", "std::vector< int >", &_blip_ID);
    _tree->Branch("blip_isValid", "std::vector< bool >", &_blip_isValid);
    _tree->Branch("blip_TPC", "std::vector< int >", &_blip_TPC);
    _tree->Branch("blip_NPlanes", "std::vector< int >", &_blip_NPlanes);
    _tree->Branch("blip_MaxWireSpan", "std::vector< int >", &_blip_MaxWireSpan);
    _tree->Branch("blip_Energy", "std::vector< float >", &_blip_Energy);
    _tree->Branch("blip_EnergyESTAR", "std::vector< float >", &_blip_EnergyESTAR);
    _tree->Branch("blip_Time", "std::vector< float >", &_blip_Time);
    _tree->Branch("blip_ProxTrkDist", "std::vector< float >", &_blip_ProxTrkDist);
    _tree->Branch("blip_ProxTrkID", "std::vector< int >", &_blip_ProxTrkID);
    _tree->Branch("blip_inCylinder", "std::vector< bool >", &_blip_inCylinder);
    _tree->Branch("blip_X", "std::vector< float >", &_blip_X);
    _tree->Branch("blip_Y", "std::vector< float >", &_blip_Y);
    _tree->Branch("blip_Z", "std::vector< float >", &_blip_Z);
    _tree->Branch("blip_SigmaYZ", "std::vector< float >", &_blip_SigmaYZ);
    _tree->Branch("blip_dX", "std::vector< float >", &_blip_dX);
    _tree->Branch("blip_dYZ", "std::vector< float >", &_blip_dYZ);
    _tree->Branch("blip_Charge", "std::vector< float >", &_blip_Charge);
    _tree->Branch("blip_LeadG4ID", "std::vector< int >", &_blip_LeadG4ID);
    _tree->Branch("blip_pdg", "std::vector< int >", &_blip_pdg);
    _tree->Branch("blip_process", "std::vector< std::string >", &_blip_process);
    _tree->Branch("blip_vx", "std::vector< float >", &_blip_vx);
    _tree->Branch("blip_vy", "std::vector< float >", &_blip_vy);
    _tree->Branch("blip_vz", "std::vector< float >", &_blip_vz);
    _tree->Branch("blip_E", "std::vector< float >", &_blip_E);
    _tree->Branch("blip_mass", "std::vector< float >", &_blip_mass);
    _tree->Branch("blip_trkid", "std::vector< int >", &_blip_trkid);
}

void BlipAnalysis::resetTTree(TTree *_tree) {
    _blip_ID.clear();
    _blip_isValid.clear();
    _blip_TPC.clear();
    _blip_NPlanes.clear();
    _blip_MaxWireSpan.clear();
    _blip_Energy.clear();
    _blip_EnergyESTAR.clear();
    _blip_Time.clear();
    _blip_ProxTrkDist.clear();
    _blip_ProxTrkID.clear();
    _blip_inCylinder.clear();
    _blip_X.clear();
    _blip_Y.clear();
    _blip_Z.clear();
    _blip_SigmaYZ.clear();
    _blip_dX.clear();
    _blip_dYZ.clear();
    _blip_Charge.clear();
    _blip_LeadG4ID.clear();
    _blip_pdg.clear();
    _blip_process.clear();
    _blip_vx.clear();
    _blip_vy.clear();
    _blip_vz.clear();
    _blip_E.clear();
    _blip_mass.clear();
    _blip_trkid.clear();
}

DEFINE_ART_CLASS_TOOL(BlipAnalysis)
}

#endif