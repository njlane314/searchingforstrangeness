#ifndef METAANALYSIS_TOOL_CXX
#define METAANALYSIS_TOOL_CXX

#include "AnalysisTools/AnalysisToolBase.h"

#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpFlash.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"

#include "larsim/EventWeight/Base/MCEventWeight.h"

#include "TTree.h"

#include <cctype>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace analysis {

class MetaAnalysis_tool : public AnalysisToolBase {
public:
  explicit MetaAnalysis_tool(const fhicl::ParameterSet& p) { configure(p); }
  ~MetaAnalysis_tool() override = default;

  void configure(const fhicl::ParameterSet&) override {}
  void analyseEvent(const art::Event& event, bool is_data) override;
  void analyseSlice(const art::Event&, std::vector<common::ProxyPfpElem_t>&, bool, bool) override {}
  void setBranches(TTree* tree) override;
  void resetTTree(TTree* tree) override;

private:
  bool _is_data = false;
  int _run = -1;
  int _sub = -1;
  int _evt = -1;
  uint64_t _event_time_value = 0ULL;

  std::vector<std::string> _process_names;
  std::vector<std::string> _process_release_versions;
  std::vector<std::string> _process_pset_ids;

  std::string _generator_module_label;
  std::string _fhicl_pset_id_generator;
  std::string _genie_version;
  std::string _genie_tune;
  std::string _genie_knobs_fcl;

  std::string _flux_module_label;
  std::string _fhicl_pset_id_flux;
  std::string _flux_tag;
  std::string _flux_release;
  std::string _flux_hist_version;
  float _horn_current = std::numeric_limits<float>::quiet_NaN();
  std::string _beam_mode;
  float _beam_energy_nominal_GeV = std::numeric_limits<float>::quiet_NaN();

  std::string _g4_module_label;
  std::string _fhicl_pset_id_g4;
  std::string _geant4_release;
  std::string _physics_list;

  std::string _reco_pfp_module_label;
  std::string _fhicl_pset_id_reco;
  std::string _pandora_version;
  std::string _pandora_settings_tag;

  std::vector<std::string> _swtrigger_algo_names;
  std::vector<int> _swtrigger_algo_passed;
  float _opfilter_pe_beam = std::numeric_limits<float>::quiet_NaN();
  float _opfilter_pe_veto = std::numeric_limits<float>::quiet_NaN();

  std::vector<std::string> _mc_weight_names;

  int _mc_nu_pdg = 0;
  int _mc_ccnc = -1;
  int _mc_mode = -1;

  std::vector<std::string> _hit_producers;
  std::vector<int> _hit_counts;

  std::vector<std::string> _slice_producers;
  std::vector<int> _slice_counts;

  std::vector<std::string> _pfp_producers;
  std::vector<int> _pfp_counts;

  std::vector<std::string> _track_producers;
  std::vector<int> _track_counts;

  std::vector<std::string> _shower_producers;
  std::vector<int> _shower_counts;

  std::vector<std::string> _opflash_producers;
  std::vector<int> _opflash_counts;

  std::vector<std::string> _crt_producers;
  std::vector<int> _crt_counts;

  template <class VecT>
  void collectCounts(const art::Event& e,
                     std::vector<std::string>& producers,
                     std::vector<int>& counts) {
    producers.clear();
    counts.clear();
    std::vector<art::Handle<VecT>> handles;
    e.getManyByType(handles);
    producers.reserve(handles.size());
    counts.reserve(handles.size());
    for (auto const& h : handles) {
      auto const* prov = h.provenance();
      std::string label = prov ? prov->productDescription().moduleLabel() : std::string{};
      producers.push_back(label);
      int n = static_cast<int>(h->size());
      counts.push_back(n);
    }
  }

  static std::string compactFcl(const fhicl::ParameterSet& ps) {
    std::string s = ps.to_string();
    std::string out;
    out.reserve(s.size());
    bool ws = false;
    for (char c : s) {
      if (c == '\n' || c == '\t' || c == '\r') {
        if (!ws) {
          out.push_back(' ');
          ws = true;
        }
        continue;
      }
      if (std::isspace(static_cast<unsigned char>(c))) {
        if (!ws) {
          out.push_back(' ');
          ws = true;
        }
      } else {
        out.push_back(c);
        ws = false;
      }
    }
    return out;
  }

  static bool get_if_present_str(const fhicl::ParameterSet& ps,
                                 std::initializer_list<const char*> keys,
                                 std::string& out) {
    for (auto k : keys) {
      if (ps.has_key(k)) {
        try {
          out = ps.get<std::string>(k);
          return true;
        } catch (...) {
        }
      }
    }
    return false;
  }

  static bool get_if_present_f(const fhicl::ParameterSet& ps,
                               std::initializer_list<const char*> keys,
                               float& out) {
    for (auto k : keys) {
      if (ps.has_key(k)) {
        try {
          out = ps.get<float>(k);
          return true;
        } catch (...) {
        }
      }
    }
    return false;
  }

  template <class CollT>
  bool fillModulePSetInfo_(const art::Handle<CollT>& h,
                           std::string& out_label,
                           std::string& out_pset_id,
                           fhicl::ParameterSet& out_pset,
                           std::string* out_release = nullptr) {
    auto const* prov = h.provenance();
    if (!prov) return false;
    auto const& desc = prov->productDescription();
    out_label = desc.moduleLabel();
    out_pset_id = desc.parameterSetID().to_string();
    if (out_release) *out_release = desc.processConfiguration().releaseVersion();

    try {
      fhicl::ParameterSet ptmp;
      if (fhicl::ParameterSetRegistry::get(fhicl::ParameterSetID{out_pset_id}, ptmp)) {
        out_pset = std::move(ptmp);
        return true;
      }
    } catch (...) {
    }
    return false;
  }
};

void MetaAnalysis_tool::setBranches(TTree* t)
{
  t->Branch("is_data", &_is_data, "is_data/O");
  t->Branch("run", &_run, "run/I");
  t->Branch("sub", &_sub, "sub/I");
  t->Branch("evt", &_evt, "evt/I");
  t->Branch("event_time_value", &_event_time_value, "event_time_value/l");

  t->Branch("process_names", "std::vector<std::string>", &_process_names);
  t->Branch("process_release_versions", "std::vector<std::string>", &_process_release_versions);
  t->Branch("process_pset_ids", "std::vector<std::string>", &_process_pset_ids);

  t->Branch("generator_module_label", &_generator_module_label);
  t->Branch("fhicl_pset_id_generator", &_fhicl_pset_id_generator);
  t->Branch("genie_version", &_genie_version);
  t->Branch("genie_tune", &_genie_tune);
  t->Branch("genie_knobs_fcl", &_genie_knobs_fcl);

  t->Branch("flux_module_label", &_flux_module_label);
  t->Branch("fhicl_pset_id_flux", &_fhicl_pset_id_flux);
  t->Branch("flux_tag", &_flux_tag);
  t->Branch("flux_release", &_flux_release);
  t->Branch("flux_hist_version", &_flux_hist_version);
  t->Branch("horn_current", &_horn_current, "horn_current/F");
  t->Branch("beam_mode", &_beam_mode);
  t->Branch("beam_energy_nominal_GeV", &_beam_energy_nominal_GeV, "beam_energy_nominal_GeV/F");

  t->Branch("g4_module_label", &_g4_module_label);
  t->Branch("fhicl_pset_id_g4", &_fhicl_pset_id_g4);
  t->Branch("geant4_release", &_geant4_release);
  t->Branch("physics_list", &_physics_list);

  t->Branch("reco_pfp_module_label", &_reco_pfp_module_label);
  t->Branch("fhicl_pset_id_reco", &_fhicl_pset_id_reco);
  t->Branch("pandora_version", &_pandora_version);
  t->Branch("pandora_settings_tag", &_pandora_settings_tag);

  t->Branch("swtrigger_algo_names", "std::vector<std::string>", &_swtrigger_algo_names);
  t->Branch("swtrigger_algo_passed", "std::vector<int>", &_swtrigger_algo_passed);
  t->Branch("opfilter_pe_beam", &_opfilter_pe_beam, "opfilter_pe_beam/F");
  t->Branch("opfilter_pe_veto", &_opfilter_pe_veto, "opfilter_pe_veto/F");

  t->Branch("mc_weight_names", "std::vector<std::string>", &_mc_weight_names);

  t->Branch("mc_nu_pdg", &_mc_nu_pdg, "mc_nu_pdg/I");
  t->Branch("mc_ccnc", &_mc_ccnc, "mc_ccnc/I");
  t->Branch("mc_mode", &_mc_mode, "mc_mode/I");

  t->Branch("hit_producers", "std::vector<std::string>", &_hit_producers);
  t->Branch("hit_counts", "std::vector<int>", &_hit_counts);
  t->Branch("slice_producers", "std::vector<std::string>", &_slice_producers);
  t->Branch("slice_counts", "std::vector<int>", &_slice_counts);
  t->Branch("pfp_producers", "std::vector<std::string>", &_pfp_producers);
  t->Branch("pfp_counts", "std::vector<int>", &_pfp_counts);
  t->Branch("track_producers", "std::vector<std::string>", &_track_producers);
  t->Branch("track_counts", "std::vector<int>", &_track_counts);
  t->Branch("shower_producers", "std::vector<std::string>", &_shower_producers);
  t->Branch("shower_counts", "std::vector<int>", &_shower_counts);
  t->Branch("opflash_producers", "std::vector<std::string>", &_opflash_producers);
  t->Branch("opflash_counts", "std::vector<int>", &_opflash_counts);
  t->Branch("crt_producers", "std::vector<std::string>", &_crt_producers);
  t->Branch("crt_counts", "std::vector<int>", &_crt_counts);
}

void MetaAnalysis_tool::resetTTree(TTree*)
{
  _is_data = false;
  _run = _sub = _evt = -1;
  _event_time_value = 0ULL;

  _process_names.clear();
  _process_release_versions.clear();
  _process_pset_ids.clear();

  _generator_module_label.clear();
  _fhicl_pset_id_generator.clear();
  _genie_version.clear();
  _genie_tune.clear();
  _genie_knobs_fcl.clear();

  _flux_module_label.clear();
  _fhicl_pset_id_flux.clear();
  _flux_tag.clear();
  _flux_release.clear();
  _flux_hist_version.clear();
  _horn_current = std::numeric_limits<float>::quiet_NaN();
  _beam_mode.clear();
  _beam_energy_nominal_GeV = std::numeric_limits<float>::quiet_NaN();

  _g4_module_label.clear();
  _fhicl_pset_id_g4.clear();
  _geant4_release.clear();
  _physics_list.clear();

  _reco_pfp_module_label.clear();
  _fhicl_pset_id_reco.clear();
  _pandora_version.clear();
  _pandora_settings_tag.clear();

  _swtrigger_algo_names.clear();
  _swtrigger_algo_passed.clear();
  _opfilter_pe_beam = std::numeric_limits<float>::quiet_NaN();
  _opfilter_pe_veto = std::numeric_limits<float>::quiet_NaN();

  _mc_weight_names.clear();

  _mc_nu_pdg = 0;
  _mc_ccnc = -1;
  _mc_mode = -1;

  _hit_producers.clear();
  _hit_counts.clear();
  _slice_producers.clear();
  _slice_counts.clear();
  _pfp_producers.clear();
  _pfp_counts.clear();
  _track_producers.clear();
  _track_counts.clear();
  _shower_producers.clear();
  _shower_counts.clear();
  _opflash_producers.clear();
  _opflash_counts.clear();
  _crt_producers.clear();
  _crt_counts.clear();
}

void MetaAnalysis_tool::analyseEvent(const art::Event& e, bool is_data)
{
  _is_data = is_data;
  _run = static_cast<int>(e.run());
  _sub = static_cast<int>(e.subRun());
  _evt = static_cast<int>(e.event());
  try {
    _event_time_value = e.time().value();
  } catch (...) {
    _event_time_value = 0ULL;
  }

  {
    _process_names.clear();
    _process_release_versions.clear();
    _process_pset_ids.clear();
    auto const& ph = e.processHistory();
    _process_names.reserve(ph.size());
    _process_release_versions.reserve(ph.size());
    _process_pset_ids.reserve(ph.size());
    for (auto const& pc : ph) {
      _process_names.emplace_back(pc.processName());
      _process_release_versions.emplace_back(pc.releaseVersion());
      _process_pset_ids.emplace_back(pc.parameterSetID().to_string());
    }
  }

  {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> gens;
    e.getManyByType(gens);
    for (auto const& h : gens) {
      fhicl::ParameterSet gps;
      std::string pset_id;
      std::string rel;
      if (!fillModulePSetInfo_(h, _generator_module_label, _fhicl_pset_id_generator, gps, &rel)) continue;

      get_if_present_str(gps, {"GenieVersion", "GENIEVersion", "genie_version"}, _genie_version);
      get_if_present_str(gps, {"TuneName", "tune_name", "GenieTune", "genie_tune", "tune"}, _genie_tune);
      if (gps.has_key("genie")) {
        try {
          auto g = gps.get<fhicl::ParameterSet>("genie");
          get_if_present_str(g, {"version", "GenieVersion"}, _genie_version);
          get_if_present_str(g, {"tune", "TuneName", "tune_name"}, _genie_tune);
          if (g.has_key("reweight")) {
            auto rw = g.get<fhicl::ParameterSet>("reweight");
            _genie_knobs_fcl = compactFcl(rw);
          }
        } catch (...) {
        }
      }

      if (!_mc_nu_pdg && !h->empty()) {
        for (auto const& mct : *h) {
          if (!mct.NeutrinoSet()) continue;
          auto const& nu = mct.GetNeutrino();
          _mc_nu_pdg = nu.Nu().PdgCode();
          _mc_ccnc = nu.CCNC();
          _mc_mode = nu.Mode();
          break;
        }
      }

      break;
    }
  }

  {
    std::vector<art::Handle<std::vector<simb::MCFlux>>> fluxes;
    e.getManyByType(fluxes);
    for (auto const& h : fluxes) {
      fhicl::ParameterSet fps;
      std::string pset_id;
      std::string rel;
      if (!fillModulePSetInfo_(h, _flux_module_label, _fhicl_pset_id_flux, fps, &rel)) continue;

      get_if_present_str(fps, {"flux_tag", "FluxTag"}, _flux_tag);
      get_if_present_str(fps, {"flux_release", "FluxRelease"}, _flux_release);
      get_if_present_str(fps, {"flux_hist_version", "FluxHistVersion"}, _flux_hist_version);
      get_if_present_str(fps, {"BeamMode", "beam_mode", "nu_mode"}, _beam_mode);
      get_if_present_f(fps, {"HornCurrent", "horn_current_kA", "horncurrent"}, _horn_current);
      get_if_present_f(fps, {"BeamEnergyGeV", "beam_energy_GeV", "ProtonMomentumGeV"}, _beam_energy_nominal_GeV);
      break;
    }
  }

  {
    std::vector<art::Handle<std::vector<simb::MCParticle>>> g4s;
    e.getManyByType(g4s);
    for (auto const& h : g4s) {
      fhicl::ParameterSet g4ps;
      std::string pset_id;
      std::string rel;
      if (!fillModulePSetInfo_(h, _g4_module_label, _fhicl_pset_id_g4, g4ps, &rel)) continue;
      _geant4_release = rel;
      get_if_present_str(g4ps, {"PhysicsList", "physics_list", "G4PhysList", "Geant4PhysList"}, _physics_list);
      break;
    }
  }

  {
    std::vector<art::Handle<std::vector<recob::PFParticle>>> pfps;
    e.getManyByType(pfps);
    for (auto const& h : pfps) {
      fhicl::ParameterSet rps;
      std::string pset_id;
      std::string rel;
      if (!fillModulePSetInfo_(h, _reco_pfp_module_label, _fhicl_pset_id_reco, rps, &rel)) continue;
      get_if_present_str(rps, {"PandoraVersion", "pandora_version", "Version"}, _pandora_version);
      get_if_present_str(rps, {"PandoraSettings", "PandoraSettingsTag", "ConfigFile", "settings_file"}, _pandora_settings_tag);
      break;
    }
  }

  {
    _swtrigger_algo_names.clear();
    _swtrigger_algo_passed.clear();
    std::vector<art::Handle<raw::ubdaqSoftwareTriggerData>> trigHandles;
    e.getManyByType(trigHandles);
    for (auto const& h : trigHandles) {
      auto names = h->getListOfAlgorithms();
      _swtrigger_algo_names.assign(names.begin(), names.end());
      _swtrigger_algo_passed.reserve(_swtrigger_algo_names.size());
      for (auto const& nm : _swtrigger_algo_names) {
        int pass = 0;
        try {
          pass = h->passedAlgo(nm) ? 1 : 0;
        } catch (...) {
          pass = 0;
        }
        _swtrigger_algo_passed.push_back(pass);
      }
      break;
    }
  }

  {
    std::vector<art::Handle<uboone::UbooneOpticalFilter>> oph;
    e.getManyByType(oph);
    for (auto const& h : oph) {
      _opfilter_pe_beam = h->PE_Beam();
      _opfilter_pe_veto = h->PE_Veto();
      break;
    }
  }

  {
    _mc_weight_names.clear();
    std::set<std::string> names;
    std::vector<art::Handle<std::vector<evwgh::MCEventWeight>>> wgh;
    e.getManyByType(wgh);
    for (auto const& h : wgh) {
      for (auto const& w : *h) {
        try {
          for (auto const& kv : w.fWeight) names.insert(kv.first);
        } catch (...) {
        }
      }
      if (!names.empty()) break;
    }
    _mc_weight_names.assign(names.begin(), names.end());
  }

  collectCounts<std::vector<recob::Hit>>(e, _hit_producers, _hit_counts);
  collectCounts<std::vector<recob::Slice>>(e, _slice_producers, _slice_counts);
  collectCounts<std::vector<recob::PFParticle>>(e, _pfp_producers, _pfp_counts);
  collectCounts<std::vector<recob::Track>>(e, _track_producers, _track_counts);
  collectCounts<std::vector<recob::Shower>>(e, _shower_producers, _shower_counts);
  collectCounts<std::vector<recob::OpFlash>>(e, _opflash_producers, _opflash_counts);
  collectCounts<std::vector<crt::CRTHit>>(e, _crt_producers, _crt_counts);
}

}

DEFINE_ART_CLASS_TOOL(analysis::MetaAnalysis_tool)

#endif
