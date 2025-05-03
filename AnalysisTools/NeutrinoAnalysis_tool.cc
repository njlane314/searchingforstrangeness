#ifndef ANALYSIS_NEUTRINO_CXX
#define ANALYSIS_NEUTRINO_CXX

#include <iostream>
#include <map>
#include "AnalysisToolBase.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TVector3.h"
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Scatters.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/Pandora.h"
#include "larcore/Geometry/Geometry.h"

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
        
        int _nu_pdg = -1;
        int _ccnc = -1;
        int _mode = -1;
        int _interaction = -1;
        float _nu_e = std::numeric_limits<float>::lowest();
        float _nu_theta = std::numeric_limits<float>::lowest();
        float _nu_pt = std::numeric_limits<float>::lowest();
        int _target_nucleus = -1;
        int _hit_nucleon = -1;
        float _W = std::numeric_limits<float>::lowest();
        float _X = std::numeric_limits<float>::lowest();
        float _Y = std::numeric_limits<float>::lowest();
        float _QSqr = std::numeric_limits<float>::lowest();
        float _px = std::numeric_limits<float>::lowest();
        float _py = std::numeric_limits<float>::lowest();
        float _pz = std::numeric_limits<float>::lowest();
        float _vtx_x = std::numeric_limits<float>::lowest();
        float _vtx_y = std::numeric_limits<float>::lowest();
        float _vtx_z = std::numeric_limits<float>::lowest();
        float _vtx_u = std::numeric_limits<float>::lowest();
        float _vtx_v = std::numeric_limits<float>::lowest();
        float _vtx_w = std::numeric_limits<float>::lowest();
        bool _in_fidvol = false;

        float _fidvol_x_start;
        float _fidvol_x_end;
        float _fidvol_y_start;
        float _fidvol_y_end;
        float _fidvol_z_start;
        float _fidvol_z_end;

        bool isFiducial(const double x[3]) const;
    };

    NeutrinoAnalysis::NeutrinoAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void NeutrinoAnalysis::configure(const fhicl::ParameterSet& pset) {
        _MCTproducer = pset.get<art::InputTag>("MCTproducer", "generator");
        _fidvol_x_start = pset.get<float>("FidvolXstart", 10.0);
        _fidvol_x_end = pset.get<float>("FidvolXend", 10.0);
        _fidvol_y_start = pset.get<float>("FidvolYstart", 10.0);
        _fidvol_y_end = pset.get<float>("FidvolYend", 10.0);
        _fidvol_z_start = pset.get<float>("FidvolZstart", 10.0);
        _fidvol_z_end = pset.get<float>("FidvolZend", 50.0);
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
            
            _nu_pdg = neutrino_particle.PdgCode();
            _nu_e = neutrino_particle.Trajectory().E(0);
            _nu_theta = neutrino.Theta();
            _nu_pt = neutrino.Pt();
            _ccnc = neutrino.CCNC();
            _mode = neutrino.Mode();
            _interaction = neutrino.InteractionType();
            _target_nucleus = neutrino.Target();
            _hit_nucleon = neutrino.HitNuc();
            _W = neutrino.W();
            _X = neutrino.X();
            _Y = neutrino.Y();
            _QSqr = neutrino.QSqr();
            _px = neutrino_particle.Px();
            _py = neutrino_particle.Py();
            _pz = neutrino_particle.Pz();

            float corrected_vertex[3];
            common::True2RecoMappingXYZ(neutrino_particle.T(), neutrino_particle.Vx(), neutrino_particle.Vy(), neutrino_particle.Vz(), corrected_vertex);
            _vtx_x = corrected_vertex[0];
            _vtx_y = corrected_vertex[1];
            _vtx_z = corrected_vertex[2];

            _vtx_u = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_U).Z();
            _vtx_v = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_V).Z();
            _vtx_w = common::ProjectToWireView(corrected_vertex[0], corrected_vertex[1], corrected_vertex[2], common::TPC_VIEW_W).Z();

            double vtx[3] = {_vtx_x, _vtx_y, _vtx_z};
            _in_fidvol = isFiducial(vtx);
        }
    }

    bool NeutrinoAnalysis::isFiducial(const double x[3]) const {
        auto const& tpc = art::ServiceHandle<geo::Geometry>{}->TPC();
        std::vector<double> bnd = {0., 2. * tpc.HalfWidth(), -tpc.HalfHeight(), tpc.HalfHeight(), 0., tpc.Length()};
        bool is_x = x[0] > (bnd[0] + _fidvol_x_start) && x[0] < (bnd[1] - _fidvol_x_end);
        bool is_y = x[1] > (bnd[2] + _fidvol_y_start) && x[1] < (bnd[3] - _fidvol_y_end);
        bool is_z = x[2] > (bnd[4] + _fidvol_z_start) && x[2] < (bnd[5] - _fidvol_z_end);
        return is_x && is_y && is_z;
    }

    void NeutrinoAnalysis::setBranches(TTree* tree) {
        tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
        tree->Branch("nu_ccnc", &_ccnc, "nu_ccnc/I");
        tree->Branch("nu_mode", &_mode, "nu_mode/I");
        tree->Branch("nu_interaction", &_interaction, "nu_interaction/I");
        tree->Branch("nu_e", &_nu_e, "nu_e/F");
        tree->Branch("nu_theta", &_nu_theta, "nu_theta/F");
        tree->Branch("nu_pt", &_nu_pt, "nu_pt/F");
        tree->Branch("nu_target_nucleus", &_target_nucleus, "nu_target_nucleus/I");
        tree->Branch("nu_hit_nucleon", &_hit_nucleon, "nu_hit_nucleon/I");
        tree->Branch("nu_W", &_W, "nu_W/F");
        tree->Branch("nu_X", &_X, "nu_X/F");
        tree->Branch("nu_Y", &_Y, "nu_Y/F");
        tree->Branch("nu_QSqr", &_QSqr, "nu_QSqr/F");
        tree->Branch("nu_px", &_px, "nu_px/F");
        tree->Branch("nu_py", &_py, "nu_py/F");
        tree->Branch("nu_pz", &_pz, "nu_pz/F");
        tree->Branch("nu_vtx_x", &_vtx_x, "nu_vtx_x/F");
        tree->Branch("nu_vtx_y", &_vtx_y, "nu_vtx_y/F");
        tree->Branch("nu_vtx_z", &_vtx_z, "nu_vtx_z/F");
        tree->Branch("nu_vtx_u", &_vtx_u, "nu_vtx_u/F");
        tree->Branch("nu_vtx_v", &_vtx_v, "nu_vtx_v/F");
        tree->Branch("nu_vtx_w", &_vtx_w, "nu_vtx_w/F");
        tree->Branch("nu_in_fidvol", &_in_fidvol, "nu_in_fidvol/O");
    }

    void NeutrinoAnalysis::resetTTree(TTree* /*tree*/) {
        _nu_pdg = -1;
        _ccnc = -1;
        _mode = -1;
        _interaction = -1;
        _nu_e = std::numeric_limits<float>::lowest();
        _nu_theta = std::numeric_limits<float>::lowest();
        _nu_pt = std::numeric_limits<float>::lowest();
        _target_nucleus = -1;
        _hit_nucleon = -1;
        _W = std::numeric_limits<float>::lowest();
        _X = std::numeric_limits<float>::lowest();
        _Y = std::numeric_limits<float>::lowest();
        _QSqr = std::numeric_limits<float>::lowest();
        _px = std::numeric_limits<float>::lowest();
        _py = std::numeric_limits<float>::lowest();
        _pz = std::numeric_limits<float>::lowest();
        _vtx_x = std::numeric_limits<float>::lowest();
        _vtx_y = std::numeric_limits<float>::lowest();
        _vtx_z = std::numeric_limits<float>::lowest();
        _vtx_u = std::numeric_limits<float>::lowest();
        _vtx_v = std::numeric_limits<float>::lowest();
        _vtx_w = std::numeric_limits<float>::lowest();
        _in_fidvol = false;
    }

    DEFINE_ART_CLASS_TOOL(NeutrinoAnalysis)
}

#endif