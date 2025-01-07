#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/MCBase/MCShower.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "CommonFunctions/Types.h"
#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Visualisation.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"

#include "SignatureTools/SignatureToolBase.h"
#include "SignatureTools/VertexToolBase.h"

#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <tuple>

class VisualiseEventFilter : public art::EDFilter
{
public:
    explicit VisualiseEventFilter(fhicl::ParameterSet const &pset);

    VisualiseEventFilter(VisualiseEventFilter const &) = delete;
    VisualiseEventFilter(VisualiseEventFilter &&) = delete;
    VisualiseEventFilter &operator=(VisualiseEventFilter const &) = delete;
    VisualiseEventFilter &operator=(VisualiseEventFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag;

    std::string _mode;
    std::vector<std::tuple<int, int, int>> _target_events;

    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;
};

VisualiseEventFilter::VisualiseEventFilter(fhicl::ParameterSet const &pset)
    : EDFilter{pset}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _mode{pset.get<std::string>("Mode", "nominal")}
{
    if (pset.has_key("TargetEvents")) {
        for (auto const &entry : pset.get<std::vector<std::vector<int>>>("TargetEvents")) {
            if (entry.size() == 3) 
                _target_events.emplace_back(entry[0], entry[1], entry[2]);
        }
    } 
    else if (pset.has_key("TargetEventsFile")) {
        std::ifstream file(pset.get<std::string>("TargetEventsFile"));
        int run, subrun, event;
        while (file >> run >> subrun >> event) 
            _target_events.emplace_back(run, subrun, event);
    }

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    };
}

bool VisualiseEventFilter::filter(art::Event &e)
{
    if (_target_events.empty()) 
        return false;

    std::cout << _mode << std::endl;

    if (_mode == "target") {
        auto current_event = std::make_tuple(e.run(), e.subRun(), e.event());
        if (std::find(_target_events.begin(), _target_events.end(), current_event) == _target_events.end()) 
            return false;
    }

    signature::Pattern pattern;
    for (auto& signatureTool : _signatureToolsVec) {
        signature::Signature signature;
        if (!signatureTool->constructSignature(e, signature))
            return false;

        pattern.push_back(signature);
    }

    std::string filename = "event_" + std::to_string(e.run()) + "_" + std::to_string(e.subRun()) + "_" + std::to_string(e.event());
    common::visualiseSignature(e, _MCPproducer, _HitProducer, _BacktrackTag, pattern, filename);

    return true;
}

DEFINE_ART_MODULE(VisualiseEventFilter)
