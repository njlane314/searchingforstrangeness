#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TTree.h"
#include "larcore/Geometry/Geometry.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include <limits>

namespace analysis 
{
    class TruthAnalysis : public art::EDFilter {
    public:
        explicit TruthAnalysis(fhicl::ParameterSet const& p);
        virtual ~TruthAnalysis() = default;

        void configure(fhicl::ParameterSet const& p);

        void analyseEvent(art::Event const &e, bool fData) override;
        void analyseSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

        void setBranches(TTree* _tree);
        void resetTTree();

    private:
        bool isFiducial(const double x[3]) const;
        bool isActive(const double x[3]) const;

        enum class EventCategory {
            k_data = 0,
            k_external = 1,                               
            k_nc = 2,                                   
            k_nu_e_cc = 3,                
            k_nu_mu_cc_with_strange = 4,  
            k_nu_mu_cc_with_protons = 5, 
            k_nu_mu_cc_with_pions = 6, 
            k_nu_mu_cc_with_protons_pions = 7, 
            k_nu_mu_cc_other = 8,         
            k_other = 9                   
        };

        std::string fMCTproducer;

        float _mcf_nu_e;
        float _mcf_lep_e;
        int _mcf_actvol;
        int _mcf_fidvol; 
        int _event_category;

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
        _tree->Branch("event_category", &_event_category, "event_category/I");
        _tree->Branch("mcf_nu_e", &_mcf_nu_e, "mcf_nu_e/F");
        _tree->Branch("mcf_lep_e", &_mcf_lep_e, "mcf_lep_e/F");
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

    void TruthAnalysis::resetTTree() {
        _mcf_nu_e = std::numeric_limits<float>::lowest();
        _mcf_lep_e = std::numeric_limits<float>::lowest();
        _mcf_actvol = std::numeric_limits<int>::lowest();
        _mcf_fidvol = std::numeric_limits<int>::lowest(); 
        _event_category = std::numeric_limits<int>::lowest();
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

    void TruthAnalysis::analyseEvent(art::Event const& e, bool fData) {
        if (fData) {
            _event_category = static_cast<int>(EventCategory::k_data);
            return;
        }

        auto const& mct_h = e.getValidHandle<std::vector<simb::MCTruth>>(fMCTproducer);
        if (mct_h->empty()) {
            _event_category = static_cast<int>(EventCategory::k_external);
            return;
        }

        auto mct = mct_h->at(0);
        auto neutrino = mct.GetNeutrino();
        auto nu = neutrino.Nu();

        _mcf_nu_e = nu.Trajectory().E(0);
        _mcf_lep_e = neutrino.Lepton().E();

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

        if (!_mcf_actvol) {
            _event_category = static_cast<int>(EventCategory::k_external);
        } else {
            if (neutrino.CCNC() == 1) { // Neutral current
                _event_category = static_cast<int>(EventCategory::k_nc);
            } else { // Charged current
                int lepton_pdg = std::abs(neutrino.Lepton().PdgCode());
                if (lepton_pdg == 11) { // Electron-neutrino CC
                    _event_category = static_cast<int>(EventCategory::k_nu_e_cc);
                } else if (lepton_pdg == 13) { // Muon-neutrino CC
                    bool has_strange_hadron = (_mcf_nkp > 0 || _mcf_nkm > 0 || _mcf_nk0 > 0 ||
                                            _mcf_nlambda > 0 || _mcf_nsigma_p > 0 || _mcf_nsigma_0 > 0 ||
                                            _mcf_nsigma_m > 0 || _mcf_nxi_0 > 0 || _mcf_nxi_m > 0 ||
                                            _mcf_nomega > 0);
                    bool has_protons = (_mcf_npr > 0); 
                    bool has_pions = (_mcf_npp > 0 || _mcf_npm > 0);
                    bool has_pi0 = (_mcf_np0 > 0);
                    if (has_strange_hadron) {
                        _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_strange);
                    } else if (has_protons && !has_pions && !has_pi0) {
                        _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_protons); 
                    } else if (!has_protons && has_pions && !has_pi0) {
                        _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_pions); 
                    } else if (has_protons && has_pions && !has_pi0) {
                        _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_with_protons_pions);
                    } else {
                        _event_category = static_cast<int>(EventCategory::k_nu_mu_cc_other);
                    }
                } else {
                    _event_category = static_cast<int>(EventCategory::k_other);
                }
            }
        }
    }

    DEFINE_ART_MODULE(TruthAnalysis)
}