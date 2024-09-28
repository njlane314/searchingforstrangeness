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
    art::InputTag _CLSproducer; // cluster associated to PFP
    art::InputTag _SLCproducer; // slice associated to PFP
    art::InputTag _HITproducer; // hit associated to cluster
    art::InputTag _SHRproducer; // shower associated to PFP
    art::InputTag _VTXproducer; // vertex associated to PFP
    art::InputTag _TRKproducer; // track associated to PFP
    art::InputTag _PCAproducer; // PCAxis associated to PFP
    art::InputTag _MCTproducer;
    bool _is_data, _is_fake_data;
    bool _filter;
    std::string _bdt_branch;
    float _bdt_cut;

    TTree *_tree;
    int _run, _sub, _evt;
    int _selected;

    TTree *_subrun_tree;
    int _run_sr; // The run number
    int _sub_sr; // The subRun number
    float _pot;  // The total amount of POT for the current sub run

    std::map<unsigned int, unsigned int> _pfpmap;

    std::unique_ptr<::selection::SelectionToolBase> _selectionTool;
    std::vector<std::unique_ptr<::analysis::AnalysisToolBase>> _analysisToolsVec;

    void BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);

    template <typename T>
    void printPFParticleMetadata(const ProxyPfpElem_t &pfp_pxy,
                                const T &pfParticleMetadataList);

    void AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                        const ProxyPfpColl_t &pfp_pxy_col,
                        std::vector<ProxyPfpElem_t> &slice_v);

    void AddProximateParticles(const ProxyPfpElem_t &nu_pfp,
                                    const ProxyPfpColl_t &pfp_pxy_col,
                                    std::vector<ProxyPfpElem_t> &slice_v);

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
    _bdt_branch = p.get<std::string>("BDT_branch", "");
    _bdt_cut = p.get<float>("BDT_cut", -1);

    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("SelectionFilter", "Selection TTree");
    _tree->Branch("selected", &_selected, "selected/I");
    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("sub", &_sub, "sub/I");
    _tree->Branch("evt", &_evt, "evt/I");

    _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
    _subrun_tree->Branch("run", &_run_sr, "run/I");
    _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

    if ( (!_is_data) || (_is_fake_data) )
        _subrun_tree->Branch("pot", &_pot, "pot/F");

    const fhicl::ParameterSet &selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");
    _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);

    _selectionTool->setBranches(_tree);
    _selectionTool->SetData(_is_data);

    auto const tool_psets = p.get<fhicl::ParameterSet>("AnalysisTools");
    for (auto const &tool_pset_labels : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_labels);
        _analysisToolsVec.push_back(art::make_tool<::analysis::AnalysisToolBase>(tool_pset));
    }

    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->setBranches(_tree);
}

bool SelectionFilter::filter(art::Event &e)
{
    ResetTTree();

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

    BuildPFPMap(pfp_proxy);

    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
    {
        _analysisToolsVec[i]->analyzeEvent(e, _is_data); 
    }

    bool keepEvent = false;

    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy)
    {
        const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

        if (pfp_pxy->IsPrimary() == false)
            continue;

        auto PDG = fabs(pfp_pxy->PdgCode());

        if ((PDG == 12) || (PDG == 14))
        {
            printPFParticleMetadata(pfp_pxy, pfParticleMetadataList);

            std::vector<ProxyPfpElem_t> slice_pfp_v;
            AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);

            std::vector<art::Ptr<recob::Track>> sliceTracks;
            std::vector<art::Ptr<recob::Shower>> sliceShowers;

            for (auto pfp : slice_pfp_v)
            {
                auto const &ass_trk_v = pfp.get<recob::Track>();
                if (ass_trk_v.size() == 1)
                    sliceTracks.push_back(ass_trk_v.at(0));
               
                auto const &ass_shr_v = pfp.get<recob::Shower>();
                if (ass_shr_v.size() == 1)
                    sliceShowers.push_back(ass_shr_v.at(0));
            } 

            bool selected = _selectionTool->selectEvent(e, slice_pfp_v);

            if (selected)
            {
                keepEvent = true;
                _selected = 1;
            }

            for (size_t i = 0; i < _analysisToolsVec.size(); i++) 
            {
                _analysisToolsVec[i]->analyzeSlice(e, slice_pfp_v, _is_data, selected);
            }
        } // if a neutrino PFParticle
    } // for all PFParticles

    _tree->Fill();

    if (_bdt_branch != "" && _bdt_cut > 0 && _bdt_cut < 1) {
        float* bdtscore = (float*) _tree->GetBranch(_bdt_branch.c_str())->GetAddress();
        std::cout << "bdtscore=" << *bdtscore << std::endl;
        keepEvent = keepEvent && ( (*bdtscore) < _bdt_cut );
    }

    if (_filter == true)
        return keepEvent;

    return true;
}

template <typename T>
void SelectionFilter::printPFParticleMetadata(const ProxyPfpElem_t &pfp_pxy,
                                                      const T &pfParticleMetadataList)
{
    if (pfParticleMetadataList.size() != 0)
    {

        for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j)
        {

            const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
            auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
            if (!pfParticlePropertiesMap.empty())
            {
                std::cout << " Found PFParticle " << pfp_pxy->Self() << " with: " << std::endl;
                for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it)
                {
                    std::cout << "  - " << it->first << " = " << it->second << std::endl;
                }
            }
        }
    } 

    return;
}

void SelectionFilter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
{
    _pfpmap.clear();

    unsigned int p = 0;
    for (const auto &pfp_pxy : pfp_pxy_col)
    {
        _pfpmap[pfp_pxy->Self()] = p;
        p++;
    }

    return;
} 

void SelectionFilter::AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                                           const ProxyPfpColl_t &pfp_pxy_col,
                                           std::vector<ProxyPfpElem_t> &slice_v)
{
    auto daughters = pfp_pxy->Daughters();

    slice_v.push_back(pfp_pxy);

    std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;

    for (auto const &daughterid : daughters)
    {

        if (_pfpmap.find(daughterid) == _pfpmap.end())
            continue;

        auto pfp_pxy2 = pfp_pxy_col.begin();
        for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
            ++pfp_pxy2;

        AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);

    } // for all daughters

    return;
} 
void SelectionFilter::ResetTTree()
{
    _selected = 0;
    _run = std::numeric_limits<int>::lowest();
    _sub = std::numeric_limits<int>::lowest();
    _evt = std::numeric_limits<int>::lowest();

    _selectionTool->resetTTree(_tree);
    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->resetTTree(_tree);
}

bool SelectionFilter::endSubRun(art::SubRun &subrun)
{
    if ( (!_is_data) || (_is_fake_data) )
    {
        art::Handle<sumdata::POTSummary> potSummaryHandle;
        _pot = subrun.getByLabel(_MCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
        std::cout << "Storing POT info!" << std::endl;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();

    return true;
}

DEFINE_ART_MODULE(SelectionFilter)