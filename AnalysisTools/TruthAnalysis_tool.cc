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
        bool isFiducial(const double x[3]) const;
        bool isActive(const double x[3]) const;

        std::string fMCTproducer;

        bool _is_mc;

        float _mcf_nu_e;
        float _mcf_lep_e;
        int _mcf_nu_pdg;
        int _mcf_ccnc;
        int _mcf_mode;
        int _mcf_interaction;
        int _mcf_target_nucleus;
        int _mcf_hit_nucleon;
        float _mcf_W;
        float _mcf_X;
        float _mcf_Y;
        float _mcf_QSqr;
        float _mcf_px;
        float _mcf_py;
        float _mcf_pz;
        float _mcf_theta;
        float _mcf_pt;

        int _mcf_actvol;
        int _mcf_fidvol; 

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

        float m_fidvolXstart;
        float m_fidvolXend;
        float m_fidvolYstart;
        float m_fidvolYend;
        float m_fidvolZstart;
        float m_fidvolZend;
    };

    TruthAnalysis::TruthAnalysis(fhicl::ParameterSet const& p) {
        this->configure(p);
    }

    void TruthAnalysis::configure(fhicl::ParameterSet const& p) {
        fMCTproducer = p.get<std::string>("MCTproducer", "generator");
        m_fidvolXstart = p.get<float>("fidvolXstart", 10.0);
        m_fidvolXend = p.get<float>("fidvolXend", 10.0);
        m_fidvolYstart = p.get<float>("fidvolYstart", 10.0);
        m_fidvolYend = p.get<float>("fidvolYend", 10.0);
        m_fidvolZstart = p.get<float>("fidvolZstart", 10.0);
        m_fidvolZend = p.get<float>("fidvolZend", 10.0);
    }

    void TruthAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("is_mc", &_is_mc, "is_mc/B");

        _tree->Branch("mcf_nu_e", &_mcf_nu_e, "mcf_nu_e/F");
        _tree->Branch("mcf_lep_e", &_mcf_lep_e, "mcf_lep_e/F");
        _tree->Branch("mcf_nu_pdg", &_mcf_nu_pdg, "mcf_nu_pdg/I");
        _tree->Branch("mcf_ccnc", &_mcf_ccnc, "mcf_ccnc/I");
        _tree->Branch("mcf_mode", &_mcf_mode, "mcf_mode/I");
        _tree->Branch("mcf_interaction", &_mcf_interaction, "mcf_interaction/I");
        _tree->Branch("mcf_target_nucleus", &_mcf_target_nucleus, "mcf_target_nucleus/I");
        _tree->Branch("mcf_hit_nucleon", &_mcf_hit_nucleon, "mcf_hit_nucleon/I");
        _tree->Branch("mcf_W", &_mcf_W, "mcf_W/F");
        _tree->Branch("mcf_X", &_mcf_X, "mcf_X/F");
        _tree->Branch("mcf_Y", &_mcf_Y, "mcf_Y/F");
        _tree->Branch("mcf_QSqr", &_mcf_QSqr, "mcf_QSqr/F");
        _tree->Branch("mcf_px", &_mcf_px, "mcf_px/F");
        _tree->Branch("mcf_py", &_mcf_py, "mcf_py/F");
        _tree->Branch("mcf_pz", &_mcf_pz, "mcf_pz/F");
        _tree->Branch("mcf_theta", &_mcf_theta, "mcf_theta/F");
        _tree->Branch("mcf_pt", &_mcf_pt, "mcf_pt/F");

        _tree->Branch("mcf_actvol", &_mcf_actvol, "mcf_actvol/I");
        _tree->Branch("mcf_fidvol", &_mcf_fidvol, "mcf_fidvol/I"); 

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
        _is_mc = true;

        _mcf_nu_e = std::numeric_limits<float>::lowest();
        _mcf_lep_e = std::numeric_limits<float>::lowest();
        _mcf_nu_pdg = std::numeric_limits<int>::lowest();
        _mcf_ccnc = std::numeric_limits<int>::lowest();
        _mcf_mode = std::numeric_limits<int>::lowest();
        _mcf_interaction = std::numeric_limits<int>::lowest();
        _mcf_target_nucleus = std::numeric_limits<int>::lowest();
        _mcf_hit_nucleon = std::numeric_limits<int>::lowest();
        _mcf_W = std::numeric_limits<float>::lowest();
        _mcf_X = std::numeric_limits<float>::lowest();
        _mcf_Y = std::numeric_limits<float>::lowest();
        _mcf_QSqr = std::numeric_limits<float>::lowest();
        _mcf_px = std::numeric_limits<float>::lowest();
        _mcf_py = std::numeric_limits<float>::lowest();
        _mcf_pz = std::numeric_limits<float>::lowest();
        _mcf_theta = std::numeric_limits<float>::lowest();
        _mcf_pt = std::numeric_limits<float>::lowest();

        _mcf_actvol = std::numeric_limits<int>::lowest();
        _mcf_fidvol = std::numeric_limits<int>::lowest(); 
       
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

    bool TruthAnalysis::isFiducial(const double x[3]) const {
        auto const& tpc = art::ServiceHandle<geo::Geometry>()->TPC();
        geo::BoxBoundedGeo theTpcGeo = tpc.ActiveBoundingBox();
        std::vector<double> bnd = {theTpcGeo.MinX(), theTpcGeo.MaxX(), theTpcGeo.MinY(), theTpcGeo.MaxY(), theTpcGeo.MinZ(), theTpcGeo.MaxZ()};
        bool is_x = x[0] > (bnd[0] + m_fidvolXstart) && x[0] < (bnd[1] - m_fidvolXend);
        bool is_y = x[1] > (bnd[2] + m_fidvolYstart) && x[1] < (bnd[3] - m_fidvolYend);
        bool is_z = x[2] > (bnd[4] + m_fidvolZstart) && x[2] < (bnd[5] - m_fidvolZend);
        return is_x && is_y && is_z;
    }

    bool TruthAnalysis::isActive(const double x[3]) const {
        auto const& tpc = art::ServiceHandle<geo::Geometry>()->TPC();
        std::vector<double> bnd = {
            0., 2. * tpc.HalfWidth(), -tpc.HalfHeight(), tpc.HalfHeight(),
            0., tpc.Length()
        };
        bool is_x = x[0] > bnd[0] && x[0] < bnd[1];
        bool is_y = x[1] > bnd[2] && x[1] < bnd[3];
        bool is_z = x[2] > bnd[4] && x[2] < bnd[5];
        return is_x && is_y && is_z;
    }

    void TruthAnalysis::analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool is_data, bool selected) {}

    void TruthAnalysis::analyseEvent(art::Event const& e, bool is_data) {
        if (is_data) {
            _is_mc = false;
            return;
        }

        art::Handle<std::vector<simb::MCTruth>> mct_h;
        if (!e.getByLabel(fMCTproducer, mct_h) || mct_h->empty()) {
            _is_mc = false;
            return;
        }

        auto mct = mct_h->at(0);
        auto neutrino = mct.GetNeutrino();
        auto nu = neutrino.Nu();

        _mcf_nu_e = nu.Trajectory().E(0);
        _mcf_lep_e = neutrino.Lepton().E();
        _mcf_nu_pdg = nu.PdgCode();
        _mcf_ccnc = neutrino.CCNC();
        _mcf_mode = neutrino.Mode();
        _mcf_interaction = neutrino.InteractionType();
        _mcf_target_nucleus = neutrino.Target();
        _mcf_hit_nucleon = neutrino.HitNuc();
        _mcf_W = neutrino.W();
        _mcf_X = neutrino.X();
        _mcf_Y = neutrino.Y();
        _mcf_QSqr = neutrino.QSqr();
        _mcf_px = nu.Px();
        _mcf_py = nu.Py();
        _mcf_pz = nu.Pz();
        _mcf_theta = neutrino.Theta();
        _mcf_pt = neutrino.Pt();

        double vertex[3] = {nu.Vx(), nu.Vy(), nu.Vz()};
        _mcf_actvol = this->isActive(vertex);
        _mcf_fidvol = this->isFiducial(vertex); 

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

    DEFINE_ART_CLASS_TOOL(TruthAnalysis)
}

#endif