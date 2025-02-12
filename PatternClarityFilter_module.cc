#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"
#include "CommonFunctions/Visualisation.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/VertexToolBase.h"

#include "ClarityTools/ClarityToolBase.h"

#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include "TDatabasePDG.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <chrono>

class PatternClarityFilter : public art::EDFilter 
{
public:
    explicit PatternClarityFilter(fhicl::ParameterSet const &pset);

    PatternClarityFilter(PatternClarityFilter const &) = delete;
    PatternClarityFilter(PatternClarityFilter &&) = delete;
    PatternClarityFilter &operator=(PatternClarityFilter const &) = delete;
    PatternClarityFilter &operator=(PatternClarityFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag;

    const geo::GeometryCore* _geo;

    std::string _bad_channel_file;
    std::vector<bool> _bad_channel_mask;

    double _patt_hit_comp_thresh;
    int _patt_hit_thresh;
    double _sig_hit_comp_thresh; 
    int _chan_act_reg;
    double _hit_exclus_thresh;
    double _sig_exclus_thresh;

    calo::CalorimetryAlg* _calo_alg;
    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;
    std::vector<std::unique_ptr<::claritytools::ClarityToolBase>> _clarityToolsVec;
    int _targetDetectorPlane;
    bool _quickVisualise;

    bool filterPatternCompleteness(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filterSignatureIntegrity(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc);
    bool filterHitExclusivity(art::Event &e, signature::Pattern& patt, const std::vector<art::Ptr<recob::Hit>> mc_hits, const std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>& mcp_bkth_assoc); 
};

PatternClarityFilter::PatternClarityFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _patt_hit_comp_thresh{pset.get<double>("PatternHitCompletenessThreshold", 0.5)}
    , _patt_hit_thresh{pset.get<int>("PatternHitThreshold", 100)}
    , _sig_hit_comp_thresh{pset.get<double>("SignatureHitCompletenessThreshold", 0.1)}
    , _chan_act_reg{pset.get<int>("ChannelActiveRegion", 3)}
    , _hit_exclus_thresh{pset.get<double>("HitExclusivityThreshold", 0.5)}
    , _sig_exclus_thresh{pset.get<double>("SignatureExclusivityThreshold", 0.8)}
    , _targetDetectorPlane{pset.get<int>("TargetDetectorPlane", 2)}
    , _quickVisualise{pset.get<bool>("QuickVisualise", true)}
{
    _calo_alg = new calo::CalorimetryAlg(pset.get<fhicl::ParameterSet>("CaloAlg"));

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    };

    const fhicl::ParameterSet &claritytool_psets = pset.get<fhicl::ParameterSet>("ClarityTools");
    for (auto const &tool_pset_label : claritytool_psets.get_pset_names())
    {
      auto const tool_pset = claritytool_psets.get<fhicl::ParameterSet>(tool_pset_label);
      _clarityToolsVec.push_back(art::make_tool<::claritytools::ClarityToolBase>(tool_pset));
    };

}

bool PatternClarityFilter::filter(art::Event &e) 
{
    signature::Pattern patt;
    for (auto &signatureTool : _signatureToolsVec) {
        signature::Signature signature;
        if (!signatureTool->constructSignature(e, signature))
            return false;

        patt.push_back(signature);
    }

    for (auto &clarityTool : _clarityToolsVec){
      std::vector<bool> filter_result =  clarityTool->filter(e, patt, static_cast<common::PandoraView>(_targetDetectorPlane));
      if(std::find(filter_result.begin(),filter_result.end(),false) != filter_result.end()) return false;
    }

    if (_quickVisualise)
    {
        std::string filename = "event_" + std::to_string(e.run()) + "_" + std::to_string(e.subRun()) + "_" + std::to_string(e.event());
        common::visualiseTrueEvent(e, _MCPproducer, _HitProducer, _BacktrackTag, filename);
        common::visualiseSignature(e, _MCPproducer, _HitProducer, _BacktrackTag, patt, filename);
    }

    return true; 
}

DEFINE_ART_MODULE(PatternClarityFilter)
