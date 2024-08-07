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

class NeutrinoSelectionFilter;

class NeutrinoSelectionFilter : public art::EDFilter
{
public:
    explicit NeutrinoSelectionFilter(fhicl::ParameterSet const &p);

    NeutrinoSelectionFilter(NeutrinoSelectionFilter const &) = delete;
    NeutrinoSelectionFilter(NeutrinoSelectionFilter &&) = delete;
    NeutrinoSelectionFilter &operator=(NeutrinoSelectionFilter const &) = delete;
    NeutrinoSelectionFilter &operator=(NeutrinoSelectionFilter &&) = delete;

    bool filter(art::Event &e) override;

    bool endSubRun(art::SubRun &subrun) override;

    using ProxyPfpColl_t = common::ProxyPfpColl_t;
    using ProxyPfpElem_t = common::ProxyPfpElem_t;

private:
    art::InputTag fPFPproducer;
    art::InputTag fCLSproducer; // cluster associated to PFP
    art::InputTag fSLCproducer; // slice associated to PFP
    art::InputTag fHITproducer; // hit associated to cluster
    art::InputTag fSHRproducer; // shower associated to PFP
    art::InputTag fVTXproducer; // vertex associated to PFP
    art::InputTag fPCAproducer;
    art::InputTag fTRKproducer; // track associated to PFP
    art::InputTag fMCTproducer;
    bool fVerbose;
    bool fData, fFakeData;

    // TTree
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

    void ResetTTree();
};

NeutrinoSelectionFilter::NeutrinoSelectionFilter(fhicl::ParameterSet const &p)
    : EDFilter{p} 
{ 
    fPFPproducer = p.get<art::InputTag>("PFPproducer");
    fSHRproducer = p.get<art::InputTag>("SHRproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fVTXproducer = p.get<art::InputTag>("VTXproducer");
    fPCAproducer = p.get<art::InputTag>("PCAproducer");
    fTRKproducer = p.get<art::InputTag>("TRKproducer");
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fVerbose = p.get<bool>("Verbose");
    fData = p.get<bool>("IsData");
    fFakeData = p.get<bool>("IsFakeData",false);

    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("NeutrinoSelectionFilter", "Neutrino Selection TTree");
    _tree->Branch("selected", &_selected, "selected/I");
    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("sub", &_sub, "sub/I");
    _tree->Branch("evt", &_evt, "evt/I");

    _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
    _subrun_tree->Branch("run", &_run_sr, "run/I");
    _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

    if ( (!fData) || (fFakeData) )
        _subrun_tree->Branch("pot", &_pot, "pot/F");

    const fhicl::ParameterSet &selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");
    _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);

    _selectionTool->setBranches(_tree);
    _selectionTool->SetData(fData);

    auto const tool_psets = p.get<fhicl::ParameterSet>("AnalysisTools");
    for (auto const &tool_pset_labels : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_labels);
        _analysisToolsVec.push_back(art::make_tool<::analysis::AnalysisToolBase>(tool_pset));
    }

    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->setBranches(_tree);
}

bool NeutrinoSelectionFilter::filter(art::Event &e)
{
    ResetTTree();

    if (fVerbose)
    {
        std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl;
    }
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();

    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(fCLSproducer),
                                                        proxy::withAssociated<recob::Slice>(fSLCproducer),
                                                        proxy::withAssociated<recob::Track>(fTRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(fVTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(fPCAproducer),
                                                        proxy::withAssociated<recob::Shower>(fSHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(fPFPproducer));

    BuildPFPMap(pfp_proxy);

    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
    {
        _analysisToolsVec[i]->analyzeEvent(e, fData);
    }

    //bool keepEvent = false;

    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy)
    {
        //const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

        if (pfp_pxy->IsPrimary() == false)
            continue;

        auto PDG = fabs(pfp_pxy->PdgCode());

        if ((PDG == 12) || (PDG == 14) || (PDG == 16))
        {
            std::vector<ProxyPfpElem_t> slice_pfp_v;
            AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);

            if (fVerbose)
            {
                std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl;
            }

            bool selected = _selectionTool->selectEvent(e, slice_pfp_v);
            if (fVerbose && selected)
            {
                std::cout << "Slice was selected!" << std::endl;
            }

            if (selected)
            {
                //keepEvent = true;
                _selected = 1;
            }

            for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
                _analysisToolsVec[i]->analyzeSlice(e, slice_pfp_v, fData, selected);
            }
        }
    } 

    _tree->Fill();

    return true;
}

void NeutrinoSelectionFilter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
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

void NeutrinoSelectionFilter::AddDaughters(const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v)
{
    auto daughters = pfp_pxy->Daughters();

    slice_v.push_back(pfp_pxy);

    if (fVerbose)
        std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;

    for (auto const &daughterid : daughters)
    {
        if (_pfpmap.find(daughterid) == _pfpmap.end())
        {
            continue;
        }

        auto pfp_pxy2 = pfp_pxy_col.begin();
        for (size_t j = 0; j < _pfpmap.at(daughterid); ++j) 
            ++pfp_pxy2;

        AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
    } 

    return;
} 

void NeutrinoSelectionFilter::ResetTTree()
{
    _selected = 0;
    _run = std::numeric_limits<int>::lowest();
    _sub = std::numeric_limits<int>::lowest();
    _evt = std::numeric_limits<int>::lowest();

    _selectionTool->resetTTree(_tree);
    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->resetTTree(_tree);
}

bool NeutrinoSelectionFilter::endSubRun(art::SubRun &subrun)
{
    if ( (!fData) || (fFakeData) )
    {
        art::Handle<sumdata::POTSummary> potSummaryHandle;
        _pot = subrun.getByLabel(fMCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
        std::cout << "[NeutrinoSelectionFilter::endSubRun] Storing POT info!" << std::endl;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();
    return true;
}

DEFINE_ART_MODULE(NeutrinoSelectionFilter)