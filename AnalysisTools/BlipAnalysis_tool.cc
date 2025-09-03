#ifndef BLIP_ANALYSIS_CXX
#define BLIP_ANALYSIS_CXX

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
    std::vector<int>   _blip_id;
    std::vector<bool>  _blip_isvalid;
    std::vector<int>   _blip_tpc;
    std::vector<int>   _blip_nplanes;
    std::vector<int>   _blip_maxwirespan;
    std::vector<float> _blip_energy;
    std::vector<float> _blip_energyestar;
    std::vector<float> _blip_time;
    std::vector<float> _blip_proxtrkdist;
    std::vector<int>   _blip_proxtrkid;
    std::vector<bool>  _blip_incylinder;
    std::vector<float> _blip_x;
    std::vector<float> _blip_y;
    std::vector<float> _blip_z;
    std::vector<float> _blip_sigmayz;
    std::vector<float> _blip_dx;
    std::vector<float> _blip_dyz;
    std::vector<float> _blip_charge;
    std::vector<int>   _blip_leadg4id;
    std::vector<int>         _blip_pdg;
    std::vector<std::string> _blip_process;
    std::vector<float>       _blip_vx;
    std::vector<float>       _blip_vy;
    std::vector<float>       _blip_vz;
    std::vector<float>       _blip_e;
    std::vector<float>       _blip_mass;
    std::vector<int>         _blip_trkid;
};

BlipAnalysis::BlipAnalysis(const fhicl::ParameterSet &pset) {
    fBlipAlg = new blip::BlipRecoAlg(pset.get<fhicl::ParameterSet>("BlipAlg"));
}

void BlipAnalysis::configure(fhicl::ParameterSet const &pset) {}

void BlipAnalysis::analyseEvent(art::Event const &e, bool fData) {
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();
}

void BlipAnalysis::analyseSlice(art::Event const &e, std::vector<common::ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) {
    fBlipAlg->RunBlipReco(e);
    std::vector<blip::Blip> blipVec = fBlipAlg->blips;

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
        _blip_id.push_back(blipVec[i].ID);
        _blip_isvalid.push_back(blipVec[i].isValid);
        _blip_tpc.push_back(blipVec[i].TPC);
        _blip_nplanes.push_back(blipVec[i].NPlanes);
        _blip_maxwirespan.push_back(blipVec[i].MaxWireSpan);
        _blip_energy.push_back(blipVec[i].Energy);
        _blip_energyestar.push_back(blipVec[i].EnergyESTAR);
        _blip_time.push_back(blipVec[i].Time);
        _blip_proxtrkdist.push_back(blipVec[i].ProxTrkDist);
        _blip_proxtrkid.push_back(blipVec[i].ProxTrkID);
        _blip_incylinder.push_back(blipVec[i].inCylinder);
        _blip_x.push_back(blipVec[i].X());
        _blip_y.push_back(blipVec[i].Y());
        _blip_z.push_back(blipVec[i].Z());
        _blip_sigmayz.push_back(blipVec[i].SigmaYZ);
        _blip_dx.push_back(blipVec[i].dX);
        _blip_dyz.push_back(blipVec[i].dYZ);
        _blip_charge.push_back(blipVec[i].Charge);
        _blip_leadg4id.push_back(blipVec[i].truth.LeadG4ID);
        int b_pdg = -999;
        std::string b_pro = "null";
        float b_vx = -999;
        float b_vy = -999;
        float b_vz = -999;
        float b_e = -999;
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
                b_e = mcparticle[blip_to_mcpar]->E();
                b_tid = mcparticle[blip_to_mcpar]->TrackId();
            }
        }

        _blip_pdg.push_back(b_pdg);
        _blip_process.push_back(b_pro);
        _blip_vx.push_back(b_vx);
        _blip_vy.push_back(b_vy);
        _blip_vz.push_back(b_vz);
        _blip_e.push_back(b_e);
        _blip_mass.push_back(b_mass);
        _blip_trkid.push_back(b_tid);
    }
}

void BlipAnalysis::fillDefault() {
    _blip_id.push_back(std::numeric_limits<int>::lowest());
    _blip_isvalid.push_back(false);
    _blip_tpc.push_back(std::numeric_limits<int>::lowest());
    _blip_nplanes.push_back(std::numeric_limits<int>::lowest());
    _blip_maxwirespan.push_back(std::numeric_limits<int>::lowest());
    _blip_energy.push_back(std::numeric_limits<float>::lowest());
    _blip_energyestar.push_back(std::numeric_limits<float>::lowest());
    _blip_time.push_back(std::numeric_limits<float>::lowest());
    _blip_proxtrkdist.push_back(std::numeric_limits<float>::lowest());
    _blip_proxtrkid.push_back(std::numeric_limits<int>::lowest());
    _blip_incylinder.push_back(false);
    _blip_x.push_back(std::numeric_limits<float>::lowest());
    _blip_y.push_back(std::numeric_limits<float>::lowest());
    _blip_z.push_back(std::numeric_limits<float>::lowest());
    _blip_sigmayz.push_back(std::numeric_limits<float>::lowest());
    _blip_dx.push_back(std::numeric_limits<float>::lowest());
    _blip_dyz.push_back(std::numeric_limits<float>::lowest());
    _blip_charge.push_back(std::numeric_limits<float>::lowest());
    _blip_pdg.push_back(std::numeric_limits<int>::lowest());
    _blip_process.push_back("null");
    _blip_vx.push_back(std::numeric_limits<float>::lowest());
    _blip_vy.push_back(std::numeric_limits<float>::lowest());
    _blip_vz.push_back(std::numeric_limits<float>::lowest());
    _blip_e.push_back(std::numeric_limits<float>::lowest());
    _blip_mass.push_back(std::numeric_limits<float>::lowest());
    _blip_trkid.push_back(std::numeric_limits<int>::lowest());
}

void BlipAnalysis::setBranches(TTree *_tree) {
    _tree->Branch("blip_id", "std::vector< int >", &_blip_id);
    _tree->Branch("blip_isvalid", "std::vector< bool >", &_blip_isvalid);
    _tree->Branch("blip_tpc", "std::vector< int >", &_blip_tpc);
    _tree->Branch("blip_nplanes", "std::vector< int >", &_blip_nplanes);
    _tree->Branch("blip_maxwirespan", "std::vector< int >", &_blip_maxwirespan);
    _tree->Branch("blip_energy", "std::vector< float >", &_blip_energy);
    _tree->Branch("blip_energyestar", "std::vector< float >", &_blip_energyestar);
    _tree->Branch("blip_time", "std::vector< float >", &_blip_time);
    _tree->Branch("blip_proxtrkdist", "std::vector< float >", &_blip_proxtrkdist);
    _tree->Branch("blip_proxtrkid", "std::vector< int >", &_blip_proxtrkid);
    _tree->Branch("blip_incylinder", "std::vector< bool >", &_blip_incylinder);
    _tree->Branch("blip_x", "std::vector< float >", &_blip_x);
    _tree->Branch("blip_y", "std::vector< float >", &_blip_y);
    _tree->Branch("blip_z", "std::vector< float >", &_blip_z);
    _tree->Branch("blip_sigmayz", "std::vector< float >", &_blip_sigmayz);
    _tree->Branch("blip_dx", "std::vector< float >", &_blip_dx);
    _tree->Branch("blip_dyz", "std::vector< float >", &_blip_dyz);
    _tree->Branch("blip_charge", "std::vector< float >", &_blip_charge);
    _tree->Branch("blip_leadg4id", "std::vector< int >", &_blip_leadg4id);
    _tree->Branch("blip_pdg", "std::vector< int >", &_blip_pdg);
    _tree->Branch("blip_process", "std::vector< std::string >", &_blip_process);
    _tree->Branch("blip_vx", "std::vector< float >", &_blip_vx);
    _tree->Branch("blip_vy", "std::vector< float >", &_blip_vy);
    _tree->Branch("blip_vz", "std::vector< float >", &_blip_vz);
    _tree->Branch("blip_e", "std::vector< float >", &_blip_e);
    _tree->Branch("blip_mass", "std::vector< float >", &_blip_mass);
    _tree->Branch("blip_trkid", "std::vector< int >", &_blip_trkid);
}

void BlipAnalysis::resetTTree(TTree *_tree) {
    _blip_id.clear();
    _blip_isvalid.clear();
    _blip_tpc.clear();
    _blip_nplanes.clear();
    _blip_maxwirespan.clear();
    _blip_energy.clear();
    _blip_energyestar.clear();
    _blip_time.clear();
    _blip_proxtrkdist.clear();
    _blip_proxtrkid.clear();
    _blip_incylinder.clear();
    _blip_x.clear();
    _blip_y.clear();
    _blip_z.clear();
    _blip_sigmayz.clear();
    _blip_dx.clear();
    _blip_dyz.clear();
    _blip_charge.clear();
    _blip_leadg4id.clear();
    _blip_pdg.clear();
    _blip_process.clear();
    _blip_vx.clear();
    _blip_vy.clear();
    _blip_vz.clear();
    _blip_e.clear();
    _blip_mass.clear();
    _blip_trkid.clear();
}

DEFINE_ART_CLASS_TOOL(BlipAnalysis)

}

#endif
