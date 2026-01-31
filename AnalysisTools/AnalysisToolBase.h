#ifndef ANALYSISTOOLBASE_H
#define ANALYSISTOOLBASE_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"

#include "Common/ProxyTypes.h"

#include "TTree.h"

#include <limits>

namespace analysis 
{
    class AnalysisToolBase {
    public:
        virtual ~AnalysisToolBase() noexcept = default;

        virtual void configure(const fhicl::ParameterSet&) = 0;

        virtual void analyseEvent(const art::Event& event, bool is_data) = 0;

        virtual void analyseSlice(const art::Event& event, std::vector<common::ProxyPfpElem_t>& slice_pfp_vec, bool is_data, bool is_selected) = 0;

        virtual void setBranches(TTree* _tree) = 0;
        
        virtual void resetTTree(TTree* _tree) = 0;
    };
}

#endif // ANALYSISTOOLBASE_H
