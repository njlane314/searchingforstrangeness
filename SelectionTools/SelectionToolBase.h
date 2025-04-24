#ifndef SELECTION_TOOLBASE_H
#define SELECTION_TOOLBASE_H

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "../CommonDefs/Types.h"
#include "../CommonDefs/Image.h"

#include "TTree.h"
#include <limits>

namespace selection 
{
    class SelectionToolBase {
    public:
        virtual ~SelectionToolBase() noexcept = default;

        void configure(const fhicl::ParameterSet&){};

        virtual bool selectEvent(art::Event const& e, 
                                const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, 
                                const std::vector<image::Image>& calo_images, 
                                const std::vector<image::Image>& reco_images, 
                                const std::vector<image::Image>& label_images) = 0;

        virtual void setBranches(TTree* _tree) = 0;

        virtual void resetTTree(TTree* _tree) = 0;
        
        void SetData(bool isdata) { _isdata = isdata; }

    protected:
        bool _isdata;
    };  
} 

#endif