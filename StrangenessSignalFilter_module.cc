#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "AnalysisTools/AnalysisToolBase.h"

class StrangenessSignalFilter : public art::EDFilter
{
public:
    explicit StrangenessSignalFilter(fhicl::ParameterSet const &p);

    StrangenessSignalFilter(StrangenessSignalFilter const &) = delete;
    StrangenessSignalFilter(StrangenessSignalFilter &&) = delete;
    StrangenessSignalFilter &operator=(StrangenessSignalFilter const &) = delete;
    StrangenessSignalFilter &operator=(StrangenessSignalFilter &&) = delete;

    bool filter(art::Event &e) override;

private:
    bool fFilter;
    std::unique_ptr<::analysis::AnalysisToolBase> _analysisTool;
};

StrangenessSignalFilter::StrangenessSignalFilter(fhicl::ParameterSet const &p)
    : EDFilter{p}, fFilter(p.get<bool>("Filter", false))
{
    auto const tool_psets = p.get<fhicl::ParameterSet>("SignalDefinition");
    auto const sig_ana_pset = tool_psets.get<fhicl::ParameterSet>("signal");
    _analysisTool = art::make_tool<::analysis::AnalysisToolBase>(sig_ana_pset);
}

bool StrangenessSignalFilter::filter(art::Event &e)
{
    if (_analysisTool) 
    {
        std::cout << "Analysing signal filter..." << std::endl;
        return _analysisTool->passesAnalysis(e, false);
    }
    else 
    {
        throw cet::exception("StrangenessSignalFilter") << "SignalAnalysis tool configuration failed.";
        
        return false;
    }
}

DEFINE_ART_MODULE(StrangenessSignalFilter)
