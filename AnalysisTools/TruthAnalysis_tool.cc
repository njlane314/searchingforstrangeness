#ifndef ANALYSIS_MCFILTER_CXX
#define ANALYSIS_MCFILTER_CXX

#include <iostream>
#include "AnalysisToolBase.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "larcore/Geometry/Geometry.h"

namespace analysis
{
class TruthAnalysis : public AnalysisToolBase
{

public:
  TruthAnalysis(const fhicl::ParameterSet &pset);

  ~TruthAnalysis(){};

  void configure(fhicl::ParameterSet const &pset);

  void analyzeEvent(art::Event const &e, bool fData) override;

  void analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected) override;

  void SaveTruth(art::Event const &e);

  void setBranches(TTree *_tree) override;

  void resetTTree(TTree *_tree) override;

private:

  art::InputTag fMCTproducer;
  art::InputTag fMCRproducer;

  float _mcf_nu_e;
  float _mcf_lep_e;
  int _mcf_actvol;
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
  float _mcf_mcshr_elec_etot;
  int _mcf_pass_ccpi0;
  int _mcf_pass_ncpi0;
  int _mcf_pass_ccnopi;
  int _mcf_pass_ncnopi;
  int _mcf_pass_cccpi;
  int _mcf_pass_nccpi;
};

TruthAnalysis::TruthAnalysis(const fhicl::ParameterSet &p)
{
  fMCTproducer = p.get<art::InputTag>("MCTproducer");
  fMCRproducer = p.get<art::InputTag>("MCRproducer");
}

void TruthAnalysis::configure(fhicl::ParameterSet const &p) {}

void TruthAnalysis::analyzeEvent(art::Event const &e, bool fData)
{
  if (fData) return;

  auto const &mct_h = e.getValidHandle<std::vector<simb::MCTruth> >(fMCTproducer);

  auto mct = mct_h->at(0);
  auto neutrino = mct.GetNeutrino();
  auto nu = neutrino.Nu();

  _mcf_nu_e = nu.Trajectory().E(0);
  _mcf_lep_e = neutrino.Lepton().E();

  art::ServiceHandle<geo::Geometry const> geo;
  auto const& tpcActiveBox = geo->TPC().ActiveBoundingBox();
  _mcf_actvol = tpcActiveBox.ContainsPosition(geo::Point_t(nu.Vx(),nu.Vy(),nu.Vz()));

  size_t npart = mct.NParticles();
  for (size_t i = 0; i < npart; i++) {
    auto const &part = mct.GetParticle(i);
    if (part.StatusCode() != 1) continue;
    if (part.PdgCode() == 13) _mcf_nmm += 1;
    if (part.PdgCode() == -13) _mcf_nmp += 1;
    if (part.PdgCode() == 11) _mcf_nem += 1;
    if (part.PdgCode() == -11) _mcf_nep += 1;
    if (part.PdgCode() == 111) _mcf_np0 += 1;
    if (part.PdgCode() == 211) _mcf_npp += 1;
    if (part.PdgCode() == -211) _mcf_npm += 1;
    if (part.PdgCode() == 321) _mcf_nkp += 1;
    if (part.PdgCode() == -321) _mcf_nkm += 1;
    if (part.PdgCode() == 311) _mcf_nk0 += 1;
    if (part.PdgCode() == 2212) _mcf_npr += 1;
    if (part.PdgCode() == 2112) _mcf_nne += 1;
  }

  float maxElecMCShwEMeV = 0.;
  const std::vector<sim::MCShower> &inputMCShower = *(e.getValidHandle<std::vector<sim::MCShower> >(fMCRproducer));
  for (auto& mcs : inputMCShower) {
    if (std::abs(mcs.PdgCode())==11) {
      if (mcs.Start().E()>maxElecMCShwEMeV) {
	_mcf_mcshr_elec_etot = mcs.Start().E();
	maxElecMCShwEMeV = _mcf_mcshr_elec_etot;
      }
    }
  }

  _mcf_pass_ccpi0 = 0;
  if (_mcf_actvol==1 && _mcf_nmm==1 && _mcf_nem==0 && _mcf_nep==0 && _mcf_np0==1) _mcf_pass_ccpi0 = 1;

  _mcf_pass_ncpi0 = 0;
  if (_mcf_actvol==1 && _mcf_nmm==0 && _mcf_nmp==0 && _mcf_nem==0 && _mcf_nep==0 && _mcf_np0==1) _mcf_pass_ncpi0 = 1;

  _mcf_pass_ccnopi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==1 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && _mcf_npp==0 && _mcf_npm==0 &&
      (((_mcf_lep_e-0.105)>0.02 && _mcf_lep_e<0.3) || _mcf_mcshr_elec_etot>15)) _mcf_pass_ccnopi = 1;

  _mcf_pass_ncnopi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==0 && _mcf_nmp==0 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && _mcf_npp==0 && _mcf_npm==0 && _mcf_nu_e>0.9) _mcf_pass_ncnopi = 1;

  _mcf_pass_cccpi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==1 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && (_mcf_npp==1 || _mcf_npm==1) &&
      (((_mcf_lep_e-0.105)>0.02 && _mcf_lep_e<0.4) || _mcf_mcshr_elec_etot>35)) _mcf_pass_cccpi = 1;

  _mcf_pass_nccpi = 0;
  if (_mcf_actvol==1 && _mcf_nmm==0 && _mcf_nmp==0 && _mcf_nem==0 && _mcf_nep==0 &&
      _mcf_np0==0 && (_mcf_npp==1 || _mcf_npm==1)) _mcf_pass_nccpi = 1;
}

void TruthAnalysis::analyzeSlice(art::Event const &e, std::vector<ProxyPfpElem_t> &slice_pfp_v, bool fData, bool selected){}

void TruthAnalysis::setBranches(TTree *_tree)
{
  _tree->Branch("mcf_nu_e", &_mcf_nu_e, "mcf_nu_e/F");
  _tree->Branch("mcf_lep_e", &_mcf_lep_e, "mcf_lep_e/F");
  _tree->Branch("mcf_actvol", &_mcf_actvol, "mcf_actvol/I");
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
  _tree->Branch("mcf_mcshr_elec_etot", &_mcf_mcshr_elec_etot, "mcf_mcshr_elec_etot/F");
  _tree->Branch("mcf_pass_ccpi0", &_mcf_pass_ccpi0, "mcf_pass_ccpi0/I");
  _tree->Branch("mcf_pass_ncpi0", &_mcf_pass_ncpi0, "mcf_pass_ncpi0/I");
  _tree->Branch("mcf_pass_ccnopi", &_mcf_pass_ccnopi, "mcf_pass_ccnopi/I");
  _tree->Branch("mcf_pass_ncnopi", &_mcf_pass_ncnopi, "mcf_pass_ncnopi/I");
  _tree->Branch("mcf_pass_cccpi", &_mcf_pass_cccpi, "mcf_pass_cccpi/I");
  _tree->Branch("mcf_pass_nccpi", &_mcf_pass_nccpi, "mcf_pass_nccpi/I");
}

void TruthAnalysis::resetTTree(TTree *_tree)
{
  _mcf_nu_e = -1.;
  _mcf_lep_e = -1;
  _mcf_actvol = -1;
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
  _mcf_mcshr_elec_etot = -1.;
  _mcf_pass_ccpi0 = -1;
  _mcf_pass_ncpi0 = -1;
  _mcf_pass_ccnopi = -1;
  _mcf_pass_ncnopi = -1;
  _mcf_pass_cccpi = -1;
  _mcf_pass_nccpi = -1;
}

DEFINE_ART_CLASS_TOOL(TruthAnalysis)
} // namespace analysis

#endif