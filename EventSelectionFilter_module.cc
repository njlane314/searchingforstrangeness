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
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

class EventSelectionFilter;

class EventSelectionFilter : public art::EDFilter {
public:
    explicit EventSelectionFilter(fhicl::ParameterSet const &p);

    EventSelectionFilter(EventSelectionFilter const &) = delete;
    EventSelectionFilter(EventSelectionFilter &&) = delete;
    EventSelectionFilter &operator=(EventSelectionFilter const &) = delete;
    EventSelectionFilter &operator=(EventSelectionFilter &&) = delete;

    bool filter(art::Event &e) override;

    bool endSubRun(art::SubRun &subrun) override;

    using ProxyPfpColl_t = selection::ProxyPfpColl_t;
    using ProxyPfpElem_t = selection::ProxyPfpElem_t;

private:
    art::InputTag fPFPproducer;
    art::InputTag fCLSproducer; // cluster associated to PFP
    art::InputTag fSLCproducer; // slice associated to PFP
    art::InputTag fHITproducer; // hit associated to cluster
    art::InputTag fSHRproducer; // shower associated to PFP
    art::InputTag fVTXproducer; // vertex associated to PFP
    art::InputTag fTRKproducer; // track associated to PFP
    art::InputTag fPCAproducer; // PCAxis associated to PFP
    art::InputTag fMCTproducer;
    bool _verbose;
    bool _data, _fake_data;
    bool _filter;

    TTree *_tree;
    int _run, _sub, _evt;
    int _selected;

    TTree *_subrun_tree;
    int _run_sr; // The run number
    int _sub_sr; // The subRun number
    float _pot;  // The total amount of POT for the current sub run

    // a map linking the PFP Self() attribute used for hierarchy building to the PFP index in the event record
    std::map<unsigned int, unsigned int> _pfpmap;

    // selection tool
    std::unique_ptr<::selection::SelectionToolBase> _selectionTool;

    // analysis tool
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

EventSelectionFilter::EventSelectionFilter(fhicl::ParameterSet const &p)
    : EDFilter{p} {
    fPFPproducer = p.get<art::InputTag>("PFPproducer");
    fSHRproducer = p.get<art::InputTag>("SHRproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fVTXproducer = p.get<art::InputTag>("VTXproducer");
    fPCAproducer = p.get<art::InputTag>("PCAproducer");
    fTRKproducer = p.get<art::InputTag>("TRKproducer");
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    _verbose = p.get<bool>("Verbose");
    _data = p.get<bool>("IsData");
    _fake_data = p.get<bool>("IsFakeData",false);
    _filter = p.get<bool>("Filter", false);

    art::ServiceHandle<art::TFileService> tfs;
    _tree = tfs->make<TTree>("EventSelectionFilter", "Neutrino Selection TTree");
    _tree->Branch("selected", &_selected, "selected/I");
    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("sub", &_sub, "sub/I");
    _tree->Branch("evt", &_evt, "evt/I");

    _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
    _subrun_tree->Branch("run", &_run_sr, "run/I");
    _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

    if ( (!_data) || (_fake_data) )
        _subrun_tree->Branch("pot", &_pot, "pot/F");

    const fhicl::ParameterSet &selection_pset = p.get<fhicl::ParameterSet>("SelectionTool");
    _selectionTool = art::make_tool<::selection::SelectionToolBase>(selection_pset);

    _selectionTool->setBranches(_tree);
    _selectionTool->SetData(_data);

    auto const tool_psets = p.get<fhicl::ParameterSet>("AnalysisTools");
    for (auto const &tool_pset_labels : tool_psets.get_pset_names()) {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_labels);
        _analysisToolsVec.push_back(art::make_tool<::analysis::AnalysisToolBase>(tool_pset));
    }

    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->setBranches(_tree);
}

bool EventSelectionFilter::filter(art::Event &e) {
    this->ResetTTree();

    if (_verbose) {
      std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl;
    }
    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();

    selection::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
                              proxy::withAssociated<larpandoraobj::PFParticleMetadata>(fPFPproducer),
                              proxy::withAssociated<recob::Cluster>(fCLSproducer),
                              proxy::withAssociated<recob::Slice>(fSLCproducer),
                              proxy::withAssociated<recob::Track>(fTRKproducer),
                              proxy::withAssociated<recob::Vertex>(fVTXproducer),
                              proxy::withAssociated<recob::PCAxis>(fPCAproducer),
                              proxy::withAssociated<recob::Shower>(fSHRproducer),
                              proxy::withAssociated<recob::SpacePoint>(fPFPproducer));

    this->BuildPFPMap(pfp_proxy);

    for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
      _analysisToolsVec[i]->analyzeEvent(e, _data); 
    }

    bool keepEvent = false;

    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy) {
      const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

      if (pfp_pxy->IsPrimary() == false)
        continue;

      auto PDG = fabs(pfp_pxy->PdgCode());

      if ((PDG == 12) || (PDG == 14)) {
            if (_verbose)
                this->printPFParticleMetadata(pfp_pxy, pfParticleMetadataList);

            std::vector<ProxyPfpElem_t> slice_pfp_v;
            this->AddDaughters(pfp_pxy, pfp_proxy, slice_pfp_v);

            if (_verbose)
                std::cout << "This slice has " << slice_pfp_v.size() << " daughter PFParticles" << std::endl;

            std::vector<art::Ptr<recob::Track>> sliceTracks;
            std::vector<art::Ptr<recob::Shower>> sliceShowers;

            for (auto pfp : slice_pfp_v) {
                auto const &ass_trk_v = pfp.get<recob::Track>();
                if (ass_trk_v.size() == 1)
                    sliceTracks.push_back(ass_trk_v.at(0));
                auto const &ass_shr_v = pfp.get<recob::Shower>();
                if (ass_shr_v.size() == 1)
                    sliceShowers.push_back(ass_shr_v.at(0));
            } 

            bool selected = _selectionTool->selectEvent(e, slice_pfp_v);
            if (selected) {
                keepEvent = true;
                _selected = 1;
            }

            for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
                _analysisToolsVec[i]->analyzeSlice(e, slice_pfp_v, _data, selected);
            }
        } // if a neutrino PFParticle
    } // for all PFParticles

    _tree->Fill();

    if (_filter == true)
        return keepEvent;

    return true;
}

template <typename T>
void EventSelectionFilter::printPFParticleMetadata(const ProxyPfpElem_t &pfp_pxy,
                                                      const T &pfParticleMetadataList) {

    if (pfParticleMetadataList.size() != 0) {
        for (unsigned int j = 0; j < pfParticleMetadataList.size(); ++j) {
            const art::Ptr<larpandoraobj::PFParticleMetadata> &pfParticleMetadata(pfParticleMetadataList.at(j));
            auto pfParticlePropertiesMap = pfParticleMetadata->GetPropertiesMap();
            if (!pfParticlePropertiesMap.empty()) {
                if (_verbose)
                std::cout << " Found PFParticle " << pfp_pxy->Self() << " with: " << std::endl;
                for (std::map<std::string, float>::const_iterator it = pfParticlePropertiesMap.begin(); it != pfParticlePropertiesMap.end(); ++it) {
                    if (_verbose)
                        std::cout << "  - " << it->first << " = " << it->second << std::endl;
                }
            }
        }
    } // if PFP metadata exists
    return;
}

void EventSelectionFilter::BuildPFPMap(const ProxyPfpColl_t &pfp_pxy_col) {
    _pfpmap.clear();

    unsigned int p = 0;
    for (const auto &pfp_pxy : pfp_pxy_col) {
        _pfpmap[pfp_pxy->Self()] = p;
        p++;
    }

    return;
} // BuildPFPMap

void EventSelectionFilter::AddDaughters(const ProxyPfpElem_t &pfp_pxy,
                                           const ProxyPfpColl_t &pfp_pxy_col,
                                           std::vector<ProxyPfpElem_t> &slice_v) {

    auto daughters = pfp_pxy->Daughters();

    slice_v.push_back(pfp_pxy);

    if (_verbose)
        std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;

    for (auto const &daughterid : daughters)
    {
        if (_pfpmap.find(daughterid) == _pfpmap.end())
            continue;

        auto pfp_pxy2 = pfp_pxy_col.begin();
        for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
        ++pfp_pxy2;

        this->AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
    } // for all daughters

    return;
} // AddDaughters

void EventSelectionFilter::ResetTTree() {
    _selected = 0;
    _run = std::numeric_limits<int>::lowest();
    _sub = std::numeric_limits<int>::lowest();
    _evt = std::numeric_limits<int>::lowest();

    _selectionTool->resetTTree(_tree);
    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->resetTTree(_tree);
}

bool EventSelectionFilter::endSubRun(art::SubRun &subrun) {
    if ((!_data) || (_fake_data)) {
        art::Handle<sumdata::POTSummary> potSummaryHandle;
        _pot = subrun.getByLabel(fMCTproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();
    return true;
}

DEFINE_ART_MODULE(EventSelectionFilter)