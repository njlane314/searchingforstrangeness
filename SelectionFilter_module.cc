#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"

#include "SelectionTools/SelectionToolBase.h"
#include "AnalysisTools/AnalysisToolBase.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

#include "CommonFunctions/Geometry.h"
#include "CommonFunctions/Corrections.h"

class SelectionFilter;

class SelectionFilter : public art::EDFilter
{
    
public:
    explicit SelectionFilter(fhicl::ParameterSet const &p);

    SelectionFilter(SelectionFilter const &) = delete;
    SelectionFilter(SelectionFilter &&) = delete;
    SelectionFilter &operator=(SelectionFilter const &) = delete;
    SelectionFilter &operator=(SelectionFilter &&) = delete;

    bool filter(art::Event &e) override;
    bool endSubRun(art::SubRun &subrun) override;

    using ProxyPfpColl_t = common::ProxyPfpColl_t;
    using ProxyPfpElem_t = common::ProxyPfpElem_t;

private:
    art::InputTag _PFPproducer;
    art::InputTag _CLSproducer; 
    art::InputTag _SLCproducer; 
    art::InputTag _HITproducer; 
    art::InputTag _SHRproducer; 
    art::InputTag _VTXproducer; 
    art::InputTag _TRKproducer; 
    art::InputTag _PCAproducer; 
    art::InputTag _MCTproducer;
    bool _is_data, _is_fake_data;
    bool _filter;

    TTree *_tree;
    int _run, _sub, _evt;
    int _selected;

    TTree *_subrun_tree;
    int _run_sr; 
    int _sub_sr;
    float _pot; 

    std::map<unsigned int, unsigned int> _pfpmap;

    //std::unique_ptr<::selection::SelectionToolBase> _selectionTool;
    std::vector<std::unique_ptr<::analysis::AnalysisToolBase>> _analysisToolsVec;

    void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);
    void AddDaughters(const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v);
    void ResetTTree();
};

SelectionFilter::SelectionFilter(fhicl::ParameterSet const &p)
    : EDFilter{p}
{
    _PFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
    _SHRproducer = p.get<art::InputTag>("SHRproducer", "pandora");
    _HITproducer = p.get<art::InputTag>("HITproducer", "gaushit");
    _CLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
    _SLCproducer = p.get<art::InputTag>("SLCproducer", "pandora");
    _VTXproducer = p.get<art::InputTag>("VTXproducer", "pandora");
    _PCAproducer = p.get<art::InputTag>("PCAproducer", "pandora");
    _TRKproducer = p.get<art::InputTag>("TRKproducer", "pandora");
    _MCTproducer = p.get<art::InputTag>("MCTproducer", "generator");
    _is_data = p.get<bool>("IsData", false);
    _is_fake_data = p.get<bool>("IsFakeData",false);
    _filter = p.get<bool>("Filter", false);

    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("SelectionFilter", "Selection TTree");
    _tree->Branch("selected", &_selected, "selected/I");
    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("sub", &_sub, "sub/I");
    _tree->Branch("evt", &_evt, "evt/I");

    _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
    _subrun_tree->Branch("run", &_run_sr, "run/I");
    _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

    if ((!_is_data) || (_is_fake_data))
        _subrun_tree->Branch("pot", &_pot, "pot/F");

    /*const fhicl::ParameterSet &selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");
    _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);
    _selectionTool->setBranches(_tree);
    _selectionTool->SetData(_is_data);*/

    auto const tool_psets = p.get<fhicl::ParameterSet>("AnalysisTools");
    for (auto const &tool_pset_labels : tool_psets.get_pset_names()) {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_labels);
        _analysisToolsVec.push_back(art::make_tool<::analysis::AnalysisToolBase>(tool_pset));
    }

    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->setBranches(_tree);
}

void SelectionFilter::ResetTTree() {
    _selected = 0;
    _run = std::numeric_limits<int>::lowest();
    _sub = std::numeric_limits<int>::lowest();
    _evt = std::numeric_limits<int>::lowest();

    //_selectionTool->resetTTree(_tree);
    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->resetTTree(_tree);
}

bool SelectionFilter::filter(art::Event &e) {
    this->ResetTTree();
    std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl;
    
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();
    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, _PFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(_PFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(_CLSproducer),
                                                        proxy::withAssociated<recob::Slice>(_SLCproducer),
                                                        proxy::withAssociated<recob::Track>(_TRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(_VTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(_PCAproducer),
                                                        proxy::withAssociated<recob::Shower>(_SHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(_PFPproducer));

    this->BuildPFPMap(pfp_proxy);
    for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
        _analysisToolsVec[i]->analyseEvent(e, _is_data); 
    }

    bool keep_event = false;
    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy) {
        if (pfp_pxy->IsPrimary() == false)
            continue;

        auto nu_pdg = fabs(pfp_pxy->PdgCode());
        if ((nu_pdg == 12) || (nu_pdg == 14)) {
            std::vector<ProxyPfpElem_t> slice_pfp_v;
            this->AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);
            
            bool selected = true;
            /*bool selected = _selectionTool->selectEvent(e, slice_pfp_v);
            if (selected) {
                keep_event = true;
                _selected = 1;
            }*/

            for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
                _analysisToolsVec[i]->analyseSlice(e, slice_pfp_v, _is_data, selected);
            }
        }
    }

    _tree->Fill();
    if (_filter == true)
        return keep_event;

    return true;
}

void SelectionFilter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col) {
    _pfpmap.clear();
    unsigned int p = 0;
    for (const auto &pfp_pxy : pfp_pxy_col) {
        _pfpmap[pfp_pxy->Self()] = p;
        p++;
    }

    return;
} 

void SelectionFilter::AddDaughters(const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v) {
    auto daughters = pfp_pxy->Daughters();
    slice_v.push_back(pfp_pxy);

    for (auto const &daughterid : daughters) {
        if (_pfpmap.find(daughterid) == _pfpmap.end())
            continue;
        auto pfp_pxy2 = pfp_pxy_col.begin();
        for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
            ++pfp_pxy2;

        this->AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
    } 

    return;
} 

bool SelectionFilter::endSubRun(art::SubRun &subrun) {
    if ( (!_is_data) || (_is_fake_data) ) {
        art::Handle<sumdata::POTSummary> potSummaryHandle;
        _pot = subrun.getByLabel(_MCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();

    return true;
}

DEFINE_ART_MODULE(SelectionFilter)