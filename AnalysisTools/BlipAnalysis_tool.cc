#ifndef BLIP_ANALYSIS_CXX
#define BLIP_ANALYSIS_CXX

#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"
#include "ubreco/BlipReco/Alg/BlipRecoAlg.h"

#include "AnalysisToolBase.h"
#include "Common/BacktrackingUtilities.h"
#include "Common/EnergyCalibration.h"
#include "Common/GeometryUtils.h"
#include "Common/LLRCorrectionLookup.h"
#include "Common/LLRPIDCalculator.h"
#include "Common/LLRProtonMuonLookup.h"
#include "Common/ParticleIdentifierUtils.h"
#include "Common/ProxyTypes.h"
#include "Common/SpaceChargeCorrections.h"
#include "Common/TrackShowerScore.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <iostream>
#include <limits>

namespace analysis {

class BlipAnalysis : public AnalysisToolBase {
public:
    BlipAnalysis(const fhicl::ParameterSet &pset);
    ~BlipAnalysis(){};
    void configure(fhicl::ParameterSet const &pset);
    void analyseEvent(const art::Event &event, bool is_data) override;
    void analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) override;
    void SaveTruth(const art::Event &event);
    void fillDefault();
    void setBranches(TTree *_tree) override;
    void resetTTree(TTree *_tree) override;

private:
    int _run, _sub, _evt;
    blip::BlipRecoAlg* _blip_alg;
    std::vector<int>   _blip_id;
    std::vector<bool>  _blip_is_valid;
    std::vector<int>   _blip_tpc;
    std::vector<int>   _blip_n_planes;
    std::vector<int>   _blip_max_wire_span;
    std::vector<float> _blip_energy;
    std::vector<float> _blip_energy_estar;
    std::vector<float> _blip_time;
    std::vector<float> _blip_prox_trk_dist;
    std::vector<int>   _blip_prox_trk_id;
    std::vector<bool>  _blip_in_cylinder;
    std::vector<float> _blip_x;
    std::vector<float> _blip_y;
    std::vector<float> _blip_z;
    std::vector<float> _blip_sigma_yz;
    std::vector<float> _blip_dx;
    std::vector<float> _blip_dyz;
    std::vector<float> _blip_charge;
    std::vector<int>   _blip_lead_g4_id;
    std::vector<int>         _blip_pdg;
    std::vector<std::string> _blip_process;
    std::vector<float>       _blip_vx;
    std::vector<float>       _blip_vy;
    std::vector<float>       _blip_vz;
    std::vector<float>       _blip_e;
    std::vector<float>       _blip_mass;
    std::vector<int>         _blip_trk_id;
};

BlipAnalysis::BlipAnalysis(const fhicl::ParameterSet &pset) {
    _blip_alg = new blip::BlipRecoAlg(pset.get<fhicl::ParameterSet>("BlipAlg"));
}

void BlipAnalysis::configure(fhicl::ParameterSet const &pset) {}

void BlipAnalysis::analyseEvent(const art::Event &event, bool is_data) {
    _evt = event.event();
    _sub = event.subRun();
    _run = event.run();
}

void BlipAnalysis::analyseSlice(const art::Event &event, std::vector<common::ProxyPfpElem_t> &slice_pfp_vec, bool is_data, bool is_selected) {
    _blip_alg->RunBlipReco(event);
    std::vector<blip::Blip> blip_vec = _blip_alg->blips;

    art::Handle< std::vector<simb::MCTruth> > mc_truth_list_handle;
    art::Handle< std::vector<simb::MCParticle> > mc_particle_handle;
    std::vector<art::Ptr<simb::MCTruth> > mc_list;
    std::vector<art::Ptr<simb::MCParticle> > mc_particle;

    if (event.getByLabel("generator", mc_truth_list_handle)) {
        art::fill_ptr_vector(mc_list, mc_truth_list_handle);
    }

    if (event.getByLabel("largeant", mc_particle_handle)) {
        art::fill_ptr_vector(mc_particle, mc_particle_handle);
    }

    for (size_t i = 0; i < blip_vec.size(); i++) {
        _blip_id.push_back(blip_vec[i].ID);
        _blip_is_valid.push_back(blip_vec[i].isValid);
        _blip_tpc.push_back(blip_vec[i].TPC);
        _blip_n_planes.push_back(blip_vec[i].NPlanes);
        _blip_max_wire_span.push_back(blip_vec[i].MaxWireSpan);
        _blip_energy.push_back(blip_vec[i].Energy);
        _blip_energy_estar.push_back(blip_vec[i].EnergyESTAR);
        _blip_time.push_back(blip_vec[i].Time);
        _blip_prox_trk_dist.push_back(blip_vec[i].ProxTrkDist);
        _blip_prox_trk_id.push_back(blip_vec[i].ProxTrkID);
        _blip_in_cylinder.push_back(blip_vec[i].inCylinder);
        _blip_x.push_back(blip_vec[i].X());
        _blip_y.push_back(blip_vec[i].Y());
        _blip_z.push_back(blip_vec[i].Z());
        _blip_sigma_yz.push_back(blip_vec[i].SigmaYZ);
        _blip_dx.push_back(blip_vec[i].dX);
        _blip_dyz.push_back(blip_vec[i].dYZ);
        _blip_charge.push_back(blip_vec[i].Charge);
        _blip_lead_g4_id.push_back(blip_vec[i].truth.LeadG4ID);
        int b_pdg = 0;
        std::string b_process = "null";
        float b_vx = std::numeric_limits<float>::quiet_NaN();
        float b_vy = std::numeric_limits<float>::quiet_NaN();
        float b_vz = std::numeric_limits<float>::quiet_NaN();
        float b_e = std::numeric_limits<float>::quiet_NaN();
        float b_mass = std::numeric_limits<float>::quiet_NaN();
        int b_tid = -1;

        if (!is_data) {
            b_pdg = blip_vec[i].truth.LeadG4PDG;
            int blip_to_mc_par = -1;
            std::vector<int> neutron_trk_id;
            for (int i_mcp = 0; i_mcp < (int)mc_particle.size(); i_mcp++) {
                if (mc_particle[i_mcp]->TrackId() == blip_vec[i].truth.LeadG4ID) {
                    blip_to_mc_par = i_mcp;
                }
            }

            if (blip_to_mc_par > -1) {
                for (int i_n = 0; i_n < (int)neutron_trk_id.size(); i_n++) {
                }
                b_process = mc_particle[blip_to_mc_par]->Process();
                b_vx = mc_particle[blip_to_mc_par]->Vx();
                b_vy = mc_particle[blip_to_mc_par]->Vy();
                b_vz = mc_particle[blip_to_mc_par]->Vz();
                b_mass = mc_particle[blip_to_mc_par]->Mass();
                b_e = mc_particle[blip_to_mc_par]->E();
                b_tid = mc_particle[blip_to_mc_par]->TrackId();
            }
        }

        _blip_pdg.push_back(b_pdg);
        _blip_process.push_back(b_process);
        _blip_vx.push_back(b_vx);
        _blip_vy.push_back(b_vy);
        _blip_vz.push_back(b_vz);
        _blip_e.push_back(b_e);
        _blip_mass.push_back(b_mass);
        _blip_trk_id.push_back(b_tid);
    }
}

void BlipAnalysis::fillDefault() {
    _blip_id.push_back(-1);
    _blip_is_valid.push_back(false);
    _blip_tpc.push_back(-1);
    _blip_n_planes.push_back(-1);
    _blip_max_wire_span.push_back(-1);
    _blip_energy.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_energy_estar.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_time.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_prox_trk_dist.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_prox_trk_id.push_back(-1);
    _blip_in_cylinder.push_back(false);
    _blip_x.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_y.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_z.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_sigma_yz.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_dx.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_dyz.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_charge.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_pdg.push_back(0);
    _blip_process.push_back("null");
    _blip_vx.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_vy.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_vz.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_e.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_mass.push_back(std::numeric_limits<float>::quiet_NaN());
    _blip_trk_id.push_back(-1);
}

void BlipAnalysis::setBranches(TTree *_tree) {
    _tree->Branch("blip_id", "std::vector< int >", &_blip_id);
    _tree->Branch("blip_is_valid", "std::vector< bool >", &_blip_is_valid);
    _tree->Branch("blip_tpc", "std::vector< int >", &_blip_tpc);
    _tree->Branch("blip_n_planes", "std::vector< int >", &_blip_n_planes);
    _tree->Branch("blip_max_wire_span", "std::vector< int >", &_blip_max_wire_span);
    _tree->Branch("blip_energy", "std::vector< float >", &_blip_energy);
    _tree->Branch("blip_energy_estar", "std::vector< float >", &_blip_energy_estar);
    _tree->Branch("blip_time", "std::vector< float >", &_blip_time);
    _tree->Branch("blip_prox_trk_dist", "std::vector< float >", &_blip_prox_trk_dist);
    _tree->Branch("blip_prox_trk_id", "std::vector< int >", &_blip_prox_trk_id);
    _tree->Branch("blip_in_cylinder", "std::vector< bool >", &_blip_in_cylinder);
    _tree->Branch("blip_x", "std::vector< float >", &_blip_x);
    _tree->Branch("blip_y", "std::vector< float >", &_blip_y);
    _tree->Branch("blip_z", "std::vector< float >", &_blip_z);
    _tree->Branch("blip_sigma_yz", "std::vector< float >", &_blip_sigma_yz);
    _tree->Branch("blip_dx", "std::vector< float >", &_blip_dx);
    _tree->Branch("blip_dyz", "std::vector< float >", &_blip_dyz);
    _tree->Branch("blip_charge", "std::vector< float >", &_blip_charge);
    _tree->Branch("blip_lead_g4_id", "std::vector< int >", &_blip_lead_g4_id);
    _tree->Branch("blip_pdg", "std::vector< int >", &_blip_pdg);
    _tree->Branch("blip_process", "std::vector< std::string >", &_blip_process);
    _tree->Branch("blip_vx", "std::vector< float >", &_blip_vx);
    _tree->Branch("blip_vy", "std::vector< float >", &_blip_vy);
    _tree->Branch("blip_vz", "std::vector< float >", &_blip_vz);
    _tree->Branch("blip_e", "std::vector< float >", &_blip_e);
    _tree->Branch("blip_mass", "std::vector< float >", &_blip_mass);
    _tree->Branch("blip_trk_id", "std::vector< int >", &_blip_trk_id);
}

void BlipAnalysis::resetTTree(TTree *_tree) {
    _blip_id.clear();
    _blip_is_valid.clear();
    _blip_tpc.clear();
    _blip_n_planes.clear();
    _blip_max_wire_span.clear();
    _blip_energy.clear();
    _blip_energy_estar.clear();
    _blip_time.clear();
    _blip_prox_trk_dist.clear();
    _blip_prox_trk_id.clear();
    _blip_in_cylinder.clear();
    _blip_x.clear();
    _blip_y.clear();
    _blip_z.clear();
    _blip_sigma_yz.clear();
    _blip_dx.clear();
    _blip_dyz.clear();
    _blip_charge.clear();
    _blip_lead_g4_id.clear();
    _blip_pdg.clear();
    _blip_process.clear();
    _blip_vx.clear();
    _blip_vy.clear();
    _blip_vz.clear();
    _blip_e.clear();
    _blip_mass.clear();
    _blip_trk_id.clear();
}

DEFINE_ART_CLASS_TOOL(BlipAnalysis)

}

#endif
