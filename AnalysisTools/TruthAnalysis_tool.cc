#ifndef ANALYSIS_TRUTH_CXX
#define ANALYSIS_TRUTH_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TTree.h"
#include "larcore/Geometry/Geometry.h"
#include <limits>
#include "TVector3.h"
#include "../CommonDefs/BacktrackingFuncs.h"
#include "../CommonDefs/Geometry.h"
#include "../CommonDefs/SCECorrections.h"
#include "../CommonDefs/Pandora.h"

namespace analysis 
{
    class TruthAnalysis : public AnalysisToolBase {
    public:
        explicit TruthAnalysis(fhicl::ParameterSet const& p);
        virtual ~TruthAnalysis() = default;

        void configure(const fhicl::ParameterSet& pset) override;  

        void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) override;

        void analyseEvent(const art::Event& event, bool is_data) override;

        void setBranches(TTree* tree) override;
        
        void resetTTree(TTree* tree) override;

    private:
        std::string fMCTproducer;

        int _nu_pdg;
        int _ccnc;
        int _mode;
        int _interaction;
        float _nu_e;
        float _nu_theta;
        float _nu_pt;
        int _target_nucleus;
        int _hit_nucleon;
        float _W;
        float _X;
        float _Y;
        float _QSqr;
        float _px;
        float _py;
        float _pz;
        float _vtx_x;
        float _vtx_y;
        float _vtx_z;
        float _vtx_u;
        float _vtx_v;
        float _vtx_w;
        bool _in_fidvol;

        int _mcf_nmm;
        int _mcf_nmp;
        int _mcf_nem;
        int _mcf_nep;
        int _mcf_np0;
        int _mcf_npp;
        int _mcf_npm;
        int _mcf_nkp;
        int _mcf_nkm;
        int _mcf_nk0;
        int _mcf_npr;
        int _mcf_nne;
        int _mcf_nlambda;
        int _mcf_nsigma_p;
        int _mcf_nsigma_0;
        int _mcf_nsigma_m;
        int _mcf_nxi_0;
        int _mcf_nxi_m;
        int _mcf_nomega;

        float _fidvol_x_start;
        float _fidvol_x_end;
        float _fidvol_y_start;
        float _fidvol_y_end;
        float _fidvol_z_start;
        float _fidvol_z_end;

        bool isFiducial(const double x[3]) const;
    };

    TruthAnalysis::TruthAnalysis(fhicl::ParameterSet const& p) {
        this->configure(p);
    }

    void TruthAnalysis::configure(fhicl::ParameterSet const& p) {
        fMCTproducer = p.get<std::string>("MCTproducer", "generator");
        _fidvol_x_start = p.get<float>("fidvolXstart", 10.0);
        _fidvol_x_end = p.get<float>("fidvolXend", 10.0);
        _fidvol_y_start = p.get<float>("fidvolYstart", 10.0);
        _fidvol_y_end = p.get<float>("fidvolYend", 10.0);
        _fidvol_z_start = p.get<float>("fidvolZstart", 10.0);
        _fidvol_z_end = p.get<float>("fidvolZend", 50.0);
    }

    void TruthAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("nu_pdg", &_nu_pdg, "nu_pdg/I");
        _tree->Branch("nu_ccnc", &_ccnc, "nu_ccnc/I");
        _tree->Branch("nu_mode", &_mode, "nu_mode/I");
        _tree->Branch("nu_interaction", &_interaction, "nu_interaction/I");
        _tree->Branch("nu_e", &_nu_e, "nu_e/F");
        _tree->Branch("nu_theta", &_nu_theta, "nu_theta/F");
        _tree->Branch("nu_pt", &_nu_pt, "nu_pt/F");
        _tree->Branch("nu_target_nucleus", &_target_nucleus, "nu_target_nucleus/I");
        _tree->Branch("nu_hit_nucleon", &_hit_nucleon, "nu_hit_nucleon/I");
        _tree->Branch("nu_W", &_W, "nu_W/F");
        _tree->Branch("nu_X", &_X, "nu_X/F");
        _tree->Branch("nu_Y", &_Y, "nu_Y/F");
        _tree->Branch("nu_QSqr", &_QSqr, "nu_QSqr/F");
        _tree->Branch("nu_px", &_px, "nu_px/F");
        _tree->Branch("nu_py", &_py, "nu_py/F");
        _tree->Branch("nu_pz", &_pz, "nu_pz/F");
        _tree->Branch("nu_vtx_x", &_vtx_x, "nu_vtx_x/F");
        _tree->Branch("nu_vtx_y", &_vtx_y, "nu_vtx_y/F");
        _tree->Branch("nu_vtx_z", &_vtx_z, "nu_vtx_z/F");
        _tree->Branch("nu_vtx_u", &_vtx_u, "nu_vtx_u/F");
        _tree->Branch("nu_vtx_v", &_vtx_v, "nu_vtx_v/F");
        _tree->Branch("nu_vtx_w", &_vtx_w, "nu_vtx_w/F");
        _tree->Branch("nu_in_fidvol", &_in_fidvol, "nu_in_fidvol/O");

        _tree->Branch("mcf_nmm", &_mcf_nmm, "mcf_nmm/I");
        _tree->Branch("mcf_nmp", &_mcf_nmp, "mcf_nmp/I");
        _tree->Branch("mcf_nem", &_mcf_nem, "mcf_nem/I");
        _tree->Branch("mcf_nep", &_mcf_nep, "mcf_nep/I");
        _tree->Branch("mcf_np0", &_mcf_np0, "mcf_np0/I");
        _tree->Branch("mcf_npp", &_mcf_npp, "mcf_npp/I");
        _tree->Branch("mcf_npm", &_mcf_npm, "mcf_npm/I");
        _tree->Branch("mcf_nkp", &_mcf_nkp, "mcf_nkp/I");
        _tree->Branch("mcf_nkm", &_mcf_nkm, "mcf_nkm/I");
        _tree->Branch("mcf_nk0", &_mcf_nk0, "mcf_nk0/I");
        _tree->Branch("mcf_npr", &_mcf_npr, "mcf_npr/I");
        _tree->Branch("mcf_nne", &_mcf_nne, "mcf_nne/I");
        _tree->Branch("mcf_nlambda", &_mcf_nlambda, "mcf_nlambda/I");
        _tree->Branch("mcf_nsigma_p", &_mcf_nsigma_p, "mcf_nsigma_p/I");
        _tree->Branch("mcf_nsigma_0", &_mcf_nsigma_0, "mcf_nsigma_0/I");
        _tree->Branch("mcf_nsigma_m", &_mcf_nsigma_m, "mcf_nsigma_m/I");
        _tree->Branch("mcf_nxi_0", &_mcf_nxi_0, "mcf_nxi_0/I");
        _tree->Branch("mcf_nxi_m", &_mcf_nxi_m, "mcf_nxi_m/I");
        _tree->Branch("mcf_nomega", &_mcf_nomega, "mcf_nomega/I");
    }

    void TruthAnalysis::resetTTree(TTree* tree) {
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

        _mcf_nmm = 0;
        _mcf_nmp = 0;
        _mcf_nem = 0;
        _mcf_nep = 0;
        _mcf_np0 = 0;
        _mcf_npp = 0;
        _mcf_npm = 0;
        _mcf_nkp = 0;
        _mcf_nkm = 0;
        _mcf_nk0 = 0;
        _mcf_npr = 0;
        _mcf_nne = 0;
        _mcf_nlambda = 0;
        _mcf_nsigma_p = 0;
        _mcf_nsigma_0 = 0;
        _mcf_nsigma_m = 0;
        _mcf_nxi_0 = 0;
        _mcf_nxi_m = 0;
        _mcf_nomega = 0;
    }

    void TruthAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, 
                                    bool selected) {}

    void TruthAnalysis::analyseEvent(art::Event const& e, bool is_data) {
        if (is_data) 
            return;

        art::Handle<std::vector<simb::MCTruth>> mct_h;
        if (!e.getByLabel(fMCTproducer, mct_h) || mct_h->empty()) 
            return;

        auto mct = mct_h->at(0);
        
        if (mct.NeutrinoSet()) {
            const auto& neutrino = mct.GetNeutrino();
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
            _in_fidvol = this->isFiducial(vtx);
        } else {
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

        size_t npart = mct.NParticles();

        for (size_t i = 0; i < npart; i++) {
            auto const& part = mct.GetParticle(i);
            if (part.StatusCode() != 1) continue;
            int pdg = part.PdgCode();
            int abs_pdg = std::abs(pdg);
            if (pdg == 13) _mcf_nmm += 1;        // mu-
            else if (pdg == -13) _mcf_nmp += 1;  // mu+
            else if (pdg == 11) _mcf_nem += 1;   // e-
            else if (pdg == -11) _mcf_nep += 1;  // e+
            else if (abs_pdg == 111) _mcf_np0 += 1; // pi0
            else if (pdg == 211) _mcf_npp += 1;  // pi+
            else if (pdg == -211) _mcf_npm += 1; // pi-
            else if (pdg == 321) _mcf_nkp += 1;  // K+
            else if (pdg == -321) _mcf_nkm += 1; // K-
            else if (abs_pdg == 311) _mcf_nk0 += 1; // K0
            else if (pdg == 2212) _mcf_npr += 1; // proton
            else if (pdg == 2112) _mcf_nne += 1; // neutron
            else if (abs_pdg == 3122) _mcf_nlambda += 1; // Lambda
            else if (abs_pdg == 3222) _mcf_nsigma_p += 1; // Sigma+
            else if (abs_pdg == 3212) _mcf_nsigma_0 += 1; // Sigma0
            else if (abs_pdg == 3112) _mcf_nsigma_m += 1; // Sigma-
            else if (abs_pdg == 3322) _mcf_nxi_0 += 1; // Xi0
            else if (abs_pdg == 3312) _mcf_nxi_m += 1; // Xi-
            else if (abs_pdg == 3334) _mcf_nomega += 1; // Omega-
        }
    }

    bool TruthAnalysis::isFiducial(const double x[3]) const {
        auto const& tpc = art::ServiceHandle<geo::Geometry>{}->TPC();
        std::vector<double> bnd = {0., 2. * tpc.HalfWidth(), -tpc.HalfHeight(), tpc.HalfHeight(), 0., tpc.Length()};
        bool is_x = x[0] > (bnd[0] + _fidvol_x_start) && x[0] < (bnd[1] - _fidvol_x_end);
        bool is_y = x[1] > (bnd[2] + _fidvol_y_start) && x[1] < (bnd[3] - _fidvol_y_end);
        bool is_z = x[2] > (bnd[4] + _fidvol_z_start) && x[2] < (bnd[5] - _fidvol_z_end);
        return is_x && is_y && is_z;
    }

    DEFINE_ART_CLASS_TOOL(TruthAnalysis)
}

#endif