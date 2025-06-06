#ifndef ANALYSIS_TRUTH_CXX
#define ANALYSIS_TRUTH_CXX

#include "AnalysisToolBase.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "TTree.h"
#include <limits>
#include <vector>
#include <string>
#include "TVector3.h"

#include "../CommonDefs/BacktrackingFuncs.h" 
#include "../CommonDefs/Geometry.h"         
#include "../CommonDefs/SCECorrections.h"   

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
    };

    TruthAnalysis::TruthAnalysis(fhicl::ParameterSet const& p) {
        this->configure(p);
    }

    void TruthAnalysis::configure(fhicl::ParameterSet const& p) {
        fMCTproducer = p.get<std::string>("MCTproducer", "generator");
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
                                    bool selected) {
    }

    void TruthAnalysis::analyseEvent(art::Event const& e, bool is_data) {
        if (is_data) {
            _nu_pdg = -1; _ccnc = -1; _mode = -1; _interaction = -1;
            _nu_e = std::numeric_limits<float>::lowest(); _nu_theta = std::numeric_limits<float>::lowest(); _nu_pt = std::numeric_limits<float>::lowest();
            _target_nucleus = -1; _hit_nucleon = -1;
            _W = std::numeric_limits<float>::lowest(); _X = std::numeric_limits<float>::lowest(); _Y = std::numeric_limits<float>::lowest(); _QSqr = std::numeric_limits<float>::lowest();
            _px = std::numeric_limits<float>::lowest(); _py = std::numeric_limits<float>::lowest(); _pz = std::numeric_limits<float>::lowest();
            _mcf_nmm = 0; _mcf_nmp = 0; _mcf_nem = 0; _mcf_nep = 0; _mcf_np0 = 0; _mcf_npp = 0; _mcf_npm = 0;
            _mcf_nkp = 0; _mcf_nkm = 0; _mcf_nk0 = 0; _mcf_npr = 0; _mcf_nne = 0; _mcf_nlambda = 0;
            _mcf_nsigma_p = 0; _mcf_nsigma_0 = 0; _mcf_nsigma_m = 0; _mcf_nxi_0 = 0; _mcf_nxi_m = 0; _mcf_nomega = 0;
            return;
        }

        art::Handle<std::vector<simb::MCTruth>> mct_h;
        if (!e.getByLabel(fMCTproducer, mct_h) || mct_h->empty()) {
            _nu_pdg = -1; _ccnc = -1; _mode = -1; _interaction = -1;
            _nu_e = std::numeric_limits<float>::lowest(); _nu_theta = std::numeric_limits<float>::lowest(); _nu_pt = std::numeric_limits<float>::lowest();
            _target_nucleus = -1; _hit_nucleon = -1;
            _W = std::numeric_limits<float>::lowest(); _X = std::numeric_limits<float>::lowest(); _Y = std::numeric_limits<float>::lowest(); _QSqr = std::numeric_limits<float>::lowest();
            _px = std::numeric_limits<float>::lowest(); _py = std::numeric_limits<float>::lowest(); _pz = std::numeric_limits<float>::lowest();
            _mcf_nmm = 0; _mcf_nmp = 0; _mcf_nem = 0; _mcf_nep = 0; _mcf_np0 = 0; _mcf_npp = 0; _mcf_npm = 0;
            _mcf_nkp = 0; _mcf_nkm = 0; _mcf_nk0 = 0; _mcf_npr = 0; _mcf_nne = 0; _mcf_nlambda = 0;
            _mcf_nsigma_p = 0; _mcf_nsigma_0 = 0; _mcf_nsigma_m = 0; _mcf_nxi_0 = 0; _mcf_nxi_m = 0; _mcf_nomega = 0;
            return;
        }

        const simb::MCTruth& mct = mct_h->at(0);
        
        if (mct.NeutrinoSet()) {
            const simb::MCNeutrino& neutrino = mct.GetNeutrino();
            const simb::MCParticle& neutrino_particle = neutrino.Nu();
            _nu_pdg = neutrino_particle.PdgCode();
            _nu_e = neutrino_particle.E(); 
            _nu_theta = neutrino_particle.Momentum().Theta();
            _nu_pt = neutrino_particle.Pt();
            _ccnc = neutrino.CCNC();
            std::cout << _mode << std::endl;
            _mode = neutrino.Mode();
            _interaction = neutrino.InteractionType();
            std::cout << _interaction << std::endl;
            _target_nucleus = neutrino.Target();
            _hit_nucleon = neutrino.HitNuc();
            _W = neutrino.W();
            _X = neutrino.X();
            _Y = neutrino.Y();
            _QSqr = neutrino.QSqr();
            _px = neutrino_particle.Px();
            _py = neutrino_particle.Py();
            _pz = neutrino_particle.Pz();
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
        }

        _mcf_nmm = 0; _mcf_nmp = 0; _mcf_nem = 0; _mcf_nep = 0; _mcf_np0 = 0;
        _mcf_npp = 0; _mcf_npm = 0; _mcf_nkp = 0; _mcf_nkm = 0; _mcf_nk0 = 0;
        _mcf_npr = 0; _mcf_nne = 0; _mcf_nlambda = 0; _mcf_nsigma_p = 0;
        _mcf_nsigma_0 = 0; _mcf_nsigma_m = 0; _mcf_nxi_0 = 0; _mcf_nxi_m = 0; _mcf_nomega = 0;

        for (int i = 0; i < mct.NParticles(); ++i) {
            const simb::MCParticle& part = mct.GetParticle(i);
            if (part.StatusCode() != 1) continue; 
            
            int pdg = part.PdgCode();
            int abs_pdg = std::abs(pdg);

            if (pdg == 13) _mcf_nmm++;        
            else if (pdg == -13) _mcf_nmp++;  
            else if (pdg == 11) _mcf_nem++;   
            else if (pdg == -11) _mcf_nep++;  
            else if (pdg == 111) _mcf_np0++; 
            else if (pdg == 211) _mcf_npp++;  
            else if (pdg == -211) _mcf_npm++; 
            else if (pdg == 321) _mcf_nkp++;  
            else if (pdg == -321) _mcf_nkm++; 
            else if (pdg == 130 || pdg == 310 || pdg == 311) _mcf_nk0++; 
            else if (pdg == 2212) _mcf_npr++; 
            else if (pdg == 2112) _mcf_nne++; 
            else if (abs_pdg == 3122) _mcf_nlambda++; 
            else if (pdg == 3222) _mcf_nsigma_p++; 
            else if (pdg == 3212) _mcf_nsigma_0++; 
            else if (pdg == 3112) _mcf_nsigma_m++; 
            else if (abs_pdg == 3322) _mcf_nxi_0++; 
            else if (abs_pdg == 3312) _mcf_nxi_m++; 
            else if (abs_pdg == 3334) _mcf_nomega++;
        }
    }

    DEFINE_ART_CLASS_TOOL(TruthAnalysis)
}

#endif
