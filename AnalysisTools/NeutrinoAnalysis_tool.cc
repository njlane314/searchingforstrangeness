#ifndef ANALYSIS_NEUTRINO_CXX
#define ANALYSIS_NEUTRINO_CXX

#include <iostream>
#include <map>
#include "AnalysisToolBase.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "TVector3.h"
#include "../CommonDefs/Backtracking.h"
#include "../CommonDefs/Scores.h"
#include "../CommonDefs/Corrections.h"
#include "../CommonDefs/Scatters.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/Pandora.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/MCBase/MCShower.h"

namespace analysis
{
    class NeutrinoAnalysis : public AnalysisToolBase {
    public:
        explicit NeutrinoAnalysis(const fhicl::ParameterSet& pset);
        ~NeutrinoAnalysis() = default;

        NeutrinoAnalysis(const NeutrinoAnalysis&) = delete;
        NeutrinoAnalysis(NeutrinoAnalysis&&) = delete;
        NeutrinoAnalysis& operator=(const NeutrinoAnalysis&) = delete;
        NeutrinoAnalysis& operator=(NeutrinoAnalysis&&) = delete;

        void configure(const fhicl::ParameterSet& pset) override;  
        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;
        void analyseEvent(const art::Event& event, bool is_data) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        art::InputTag _MCTproducer;
        art::InputTag _MCPproducer;
        
        struct Neutrino {
            int nu_pdg = -1;
            int ccnc = -1;
            int mode = -1;
            int interaction = -1;
            float nu_e = std::numeric_limits<float>::lowest();
            float nu_theta = std::numeric_limits<float>::lowest();
            float nu_pt = std::numeric_limits<float>::lowest();
            int target_nucleus = -1;
            int hit_nucleon = -1;
            float W = std::numeric_limits<float>::lowest();
            float X = std::numeric_limits<float>::lowest();
            float Y = std::numeric_limits<float>::lowest();
            float QSqr = std::numeric_limits<float>::lowest();
            float px = std::numeric_limits<float>::lowest();
            float py = std::numeric_limits<float>::lowest();
            float pz = std::numeric_limits<float>::lowest();
            float vtx_x = std::numeric_limits<float>::lowest();
            float vtx_y = std::numeric_limits<float>::lowest();
            float vtx_z = std::numeric_limits<float>::lowest();
            float vtx_u = std::numeric_limits<float>::lowest();
            float vtx_v = std::numeric_limits<float>::lowest();
            float vtx_w = std::numeric_limits<float>::lowest();
        };

        Neutrino _event_neutrino; 
        std::vector<int> _fs_pdg;
        std::vector<float> _fs_energy;
        std::vector<float> _fs_px;
        std::vector<float> _fs_py; 
        std::vector<float> _fs_pz; 

        void fillNeutrino(Neutrino& neutrino_data, const simb::MCNeutrino& neutrino, const simb::MCParticle& neutrino_particle);
    };

    NeutrinoAnalysis::NeutrinoAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void NeutrinoAnalysis::configure(const fhicl::ParameterSet& pset) {
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
        _MCPproducer = pset.get<art::InputTag>("MCPproducer", "largeant");
    }

    void NeutrinoAnalysis::analyseSlice(const art::Event& /*event*/, std::vector<common::ProxyPfpElem_t>& /*slice_pfp_v*/, bool /*is_data*/, bool /*selected*/) {}

    void NeutrinoAnalysis::analyseEvent(const art::Event& event, bool is_data) {
        if (is_data) 
            return;

        auto const& truth_handle = event.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);   
        if (!truth_handle->empty() && truth_handle->at(0).NeutrinoSet()) {
            const auto& truth = truth_handle->at(0);
            const auto& neutrino = truth.GetNeutrino();
            const auto& neutrino_particle = neutrino.Nu();
            this->fillNeutrino(_event_neutrino, neutrino, neutrino_particle);
        }

        auto const& mcp_handle = event.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);
        if (!mcp_handle.isValid()) 
            return;

        for (const auto& particle : *mcp_handle) {
            if (particle.Process() == "primary") {
                _fs_pdg.push_back(particle.PdgCode());
                _fs_energy.push_back(particle.E());
                _fs_px.push_back(particle.Px());
                _fs_py.push_back(particle.Py());
                _fs_pz.push_back(particle.Pz());
            }
        }
    }

    void NeutrinoAnalysis::fillNeutrino(Neutrino& neutrino_data, const simb::MCNeutrino& neutrino, const simb::MCParticle& neutrino_particle) {
        neutrino_data.nu_pdg = neutrino_particle.PdgCode();
        neutrino_data.nu_e = neutrino_particle.Trajectory().E(0);
        neutrino_data.nu_theta = neutrino.Theta();
        neutrino_data.nu_pt = neutrino.Pt();
        neutrino_data.ccnc = neutrino.CCNC();
        neutrino_data.mode = neutrino.Mode();
        neutrino_data.interaction = neutrino.InteractionType();
        neutrino_data.target_nucleus = neutrino.Target();
        neutrino_data.hit_nucleon = neutrino.HitNuc();
        neutrino_data.W = neutrino.W();
        neutrino_data.X = neutrino.X();
        neutrino_data.Y = neutrino.Y();
        neutrino_data.QSqr = neutrino.QSqr();

        neutrino_data.px = neutrino_particle.Px();
        neutrino_data.py = neutrino_particle.Py();
        neutrino_data.pz = neutrino_particle.Pz();

        float corrected_vertex[3];
        common::True2RecoMappingXYZ(neutrino_particle.T(), neutrino_particle.Vx(), neutrino_particle.Vy(), neutrino_particle.Vz(), corrected_vertex);
        neutrino_data.vtx_x = corrected_vertex[0];
        neutrino_data.vtx_y = corrected_vertex[1];
        neutrino_data.vtx_z = corrected_vertex[2];

        neutrino_data.vtx_u = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_U).Z();
        neutrino_data.vtx_v = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_V).Z();
        neutrino_data.vtx_w = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_W).Z();
    }

    void NeutrinoAnalysis::setBranches(TTree* tree) {
        tree->Branch("nu_pdg", &_event_neutrino.nu_pdg, "nu_pdg/I");
        tree->Branch("ccnc", &_event_neutrino.ccnc, "ccnc/I");
        tree->Branch("nu_e", &_event_neutrino.nu_e, "nu_e/F");
        tree->Branch("nu_theta", &_event_neutrino.nu_theta, "nu_theta/F");
        tree->Branch("nu_pt", &_event_neutrino.nu_pt, "nu_pt/F");
        tree->Branch("interaction", &_event_neutrino.interaction, "interaction/I");
        tree->Branch("target_nucleus", &_event_neutrino.target_nucleus, "target_nucleus/I");
        tree->Branch("hit_nucleon", &_event_neutrino.hit_nucleon, "hit_nucleon/I");
        tree->Branch("W", &_event_neutrino.W, "W/F");
        tree->Branch("X", &_event_neutrino.X, "X/F");
        tree->Branch("Y", &_event_neutrino.Y, "Y/F");
        tree->Branch("QSqr", &_event_neutrino.QSqr, "QSqr/F");
        tree->Branch("vtx_x", &_event_neutrino.vtx_x, "vtx_x/F");
        tree->Branch("vtx_y", &_event_neutrino.vtx_y, "vtx_y/F");
        tree->Branch("vtx_z", &_event_neutrino.vtx_z, "vtx_z/F");
        tree->Branch("vtx_u", &_event_neutrino.vtx_u, "vtx_u/F");
        tree->Branch("vtx_v", &_event_neutrino.vtx_v, "vtx_v/F");
        tree->Branch("vtx_w", &_event_neutrino.vtx_w, "vtx_w/F");
        tree->Branch("fs_pdg", &_fs_pdg);
        tree->Branch("fs_energy", &_fs_energy);
        tree->Branch("fs_px", &_fs_px);
        tree->Branch("fs_py", &_fs_py);
        tree->Branch("fs_pz", &_fs_pz);
    }

    void NeutrinoAnalysis::resetTTree(TTree* /*tree*/) {
        _event_neutrino = Neutrino();
        _fs_pdg.clear();
        _fs_energy.clear();
        _fs_px.clear();
        _fs_py.clear();
        _fs_pz.clear();
    }

    DEFINE_ART_CLASS_TOOL(NeutrinoAnalysis)
}

#endif