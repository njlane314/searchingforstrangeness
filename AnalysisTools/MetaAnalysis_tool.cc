#ifndef METAANALYSIS_TOOL_CXX
#define METAANALYSIS_TOOL_CXX

#include "AnalysisTools/AnalysisToolBase.h"

#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h" // larpandoraobj::PFParticleMetadata

#include "TTree.h"

#include <limits>
#include <set>
#include <string>
#include <type_traits>
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
  // ---- generator (from MCTruth::MCGeneratorInfo) ----
  std::string _gen_name;
  std::string _gen_version;
  std::string _gen_tune;

  // ---- GEANT4 / largeant provenance (from MCParticle product) ----
  std::string _g4_module_label;
  std::string _g4_process_name;
  std::string _g4_release_version;
  std::string _g4_pset_id;

  // ---- PANDORA / reco provenance (from PFParticle product) ----
  std::string _pandora_module_label;
  std::string _pandora_process_name;
  std::string _pandora_release_version;
  std::string _pandora_pset_id;

  // Union of metadata keys from PFParticleMetadata objects found in the event
  std::vector<std::string> _pandora_metadata_keys;

  // ---------- helpers for MCGeneratorInfo (guard for LArSoft variations) ----------
  template <typename T, typename = void>
  struct has_GeneratorInfo : std::false_type {};
  template <typename T>
  struct has_GeneratorInfo<T, std::void_t<decltype(std::declval<const T&>().GeneratorInfo())>> : std::true_type {};

  template <typename GI, typename = void>
  struct has_GI_GeneratorName : std::false_type {};
  template <typename GI>
  struct has_GI_GeneratorName<GI, std::void_t<decltype(std::declval<const GI&>().GeneratorName())>> : std::true_type {};

  template <typename GI, typename = void>
  struct has_GI_Generator : std::false_type {};
  template <typename GI>
  struct has_GI_Generator<GI, std::void_t<decltype(std::declval<const GI&>().Generator())>> : std::true_type {};

  template <typename GI, typename = void>
  struct has_GI_Version : std::false_type {};
  template <typename GI>
  struct has_GI_Version<GI, std::void_t<decltype(std::declval<const GI&>().Version())>> : std::true_type {};

  template <typename GI, typename = void>
  struct has_GI_GeneratorVersion : std::false_type {};
  template <typename GI>
  struct has_GI_GeneratorVersion<GI, std::void_t<decltype(std::declval<const GI&>().GeneratorVersion())>> : std::true_type {};

  template <typename GI, typename = void>
  struct has_GI_Tune : std::false_type {};
  template <typename GI>
  struct has_GI_Tune<GI, std::void_t<decltype(std::declval<const GI&>().Tune())>> : std::true_type {};

  template <typename GI, typename = void>
  struct has_GI_TuneName : std::false_type {};
  template <typename GI>
  struct has_GI_TuneName<GI, std::void_t<decltype(std::declval<const GI&>().TuneName())>> : std::true_type {};

  template <typename GI>
  void extractGI_(const GI& gi, std::string& name, std::string& ver, std::string& tune) {
    if constexpr (has_GI_GeneratorName<GI>::value) name = gi.GeneratorName();
    else if constexpr (has_GI_Generator<GI>::value) name = std::to_string(static_cast<int>(gi.Generator()));
    if constexpr (has_GI_GeneratorVersion<GI>::value) ver = std::to_string(static_cast<int>(gi.GeneratorVersion()));
    else if constexpr (has_GI_Version<GI>::value)     ver = std::to_string(static_cast<int>(gi.Version()));
    if constexpr (has_GI_TuneName<GI>::value)         tune = gi.TuneName();
    else if constexpr (has_GI_Tune<GI>::value)        tune = std::to_string(static_cast<int>(gi.Tune()));
  }

  template <class CollT>
  static bool fillProvInfo_(const art::Handle<CollT>& h,
                            std::string& module_label,
                            std::string& process_name,
                            std::string& release_version,
                            std::string& pset_id)
  {
    auto const* prov = h.provenance();
    if (!prov) return false;
    module_label   = prov->branchDescription().moduleLabel();
    process_name   = prov->processName();
    release_version= prov->processConfiguration().releaseVersion();
    pset_id        = prov->parameterSetID().to_string(); // ID only; not reading FHiCL
    return true;
  }
};

// ---------------- setBranches: only new fields ----------------
void MetaAnalysis_tool::setBranches(TTree* t)
{
  // Generator info (from MCTruth::MCGeneratorInfo)
  t->Branch("gen_name",    &_gen_name);
  t->Branch("gen_version", &_gen_version);
  t->Branch("gen_tune",    &_gen_tune);

  // GEANT4 provenance (largeant / MCParticle producer)
  t->Branch("g4_module_label",   &_g4_module_label);
  t->Branch("g4_process_name",   &_g4_process_name);
  t->Branch("g4_release_version",&_g4_release_version);
  t->Branch("g4_pset_id",        &_g4_pset_id);

  // PANDORA provenance (PFParticle producer)
  t->Branch("pandora_module_label",   &_pandora_module_label);
  t->Branch("pandora_process_name",   &_pandora_process_name);
  t->Branch("pandora_release_version",&_pandora_release_version);
  t->Branch("pandora_pset_id",        &_pandora_pset_id);

  // Pandora object-level metadata keys present in the file
  t->Branch("pandora_metadata_keys", "std::vector<std::string>", &_pandora_metadata_keys);
}

// ---------------- reset ----------------
void MetaAnalysis_tool::resetTTree(TTree*)
{
  _gen_name.clear(); _gen_version.clear(); _gen_tune.clear();

  _g4_module_label.clear();
  _g4_process_name.clear();
  _g4_release_version.clear();
  _g4_pset_id.clear();

  _pandora_module_label.clear();
  _pandora_process_name.clear();
  _pandora_release_version.clear();
  _pandora_pset_id.clear();

  _pandora_metadata_keys.clear();
}

// ---------------- analyseEvent ----------------
void MetaAnalysis_tool::analyseEvent(const art::Event& e, bool /*is_data*/)
{
  // --- generator info from MCTruth::MCGeneratorInfo (first available) ---
  {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mct_handles;
    e.getManyByType(mct_handles);
    for (auto const& hmct : mct_handles) {
      if (!hmct.isValid()) continue;
      for (auto const& mct : *hmct) {
        if constexpr (has_GeneratorInfo<simb::MCTruth>::value) {
          try {
            auto const& gi = mct.GeneratorInfo();
            extractGI_(gi, _gen_name, _gen_version, _gen_tune);
            if (!_gen_name.empty() || !_gen_version.empty() || !_gen_tune.empty())
              break;
          } catch (...) {}
        }
      }
      if (!_gen_name.empty() || !_gen_version.empty() || !_gen_tune.empty()) break;
    }
  }

  // --- GEANT4 provenance (from first MCParticle collection) ---
  {
    std::vector<art::Handle<std::vector<simb::MCParticle>>> g4_handles;
    e.getManyByType(g4_handles);
    for (auto const& hg4 : g4_handles) {
      if (!hg4.isValid()) continue;
      if (fillProvInfo_(hg4, _g4_module_label, _g4_process_name, _g4_release_version, _g4_pset_id)) break;
    }
  }

  // --- PANDORA provenance (from first PFParticle collection) ---
  {
    std::vector<art::Handle<std::vector<recob::PFParticle>>> pfp_handles;
    e.getManyByType(pfp_handles);
    for (auto const& hpfp : pfp_handles) {
      if (!hpfp.isValid()) continue;
      if (fillProvInfo_(hpfp, _pandora_module_label, _pandora_process_name, _pandora_release_version, _pandora_pset_id)) break;
    }
  }

  // --- Pandora object-level metadata keys (union across products) ---
  {
    _pandora_metadata_keys.clear();
    std::set<std::string> allkeys;
    std::vector<art::Handle<std::vector<larpandoraobj::PFParticleMetadata>>> md_handles;
    e.getManyByType(md_handles);
    for (auto const& hmd : md_handles) {
      if (!hmd.isValid()) continue;
      for (auto const& md : *hmd) {
        auto const& m = md.GetPropertiesMap();
        for (auto const& kv : m) allkeys.insert(kv.first);
      }
    }
    _pandora_metadata_keys.assign(allkeys.begin(), allkeys.end());
  }
}

} // namespace analysis

DEFINE_ART_CLASS_TOOL(analysis::MetaAnalysis_tool)

#endif // METAANALYSIS_TOOL_CXX
