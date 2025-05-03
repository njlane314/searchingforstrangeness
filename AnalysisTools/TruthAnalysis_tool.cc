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
        std::string fMCTproducer;

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