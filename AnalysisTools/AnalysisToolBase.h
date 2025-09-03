#ifndef ANALYSISTOOLBASE_H
#define ANALYSISTOOLBASE_H

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "../Common/Types.h"
#include "TTree.h"
#include <limits>

namespace analysis 
{
    class AnalysisToolBase {
    public:
        virtual ~AnalysisToolBase() noexcept = default;

        virtual void configure(const fhicl::ParameterSet&) = 0;

        virtual void analyseEvent(art::Event const& e, bool _is_data) = 0;

        virtual void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) = 0;

        virtual void setBranches(TTree* _tree) = 0;
        
        virtual void resetTTree(TTree* _tree) = 0;
    };
} 

#endif