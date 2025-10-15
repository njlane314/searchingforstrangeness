#ifdef ClassDef
#undef ClassDef
#endif

// (file header elided for brevity)

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

// for counting gates via swtrigger and (optionally) optical filter
#include "ubobj/Trigger/ubdaqSoftwareTriggerData.h"
#include "ubobj/Optical/UbooneOpticalFilter.h"
#include "Rtypes.h" // Long64_t
#include "AnalysisTools/AnalysisToolBase.h"
#include "Common/ProxyTypes.h"
#include "SelectionTools/SelectionToolBase.h"

#include "TTree.h"
#include "TVector3.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

class EventSelectionFilter : public art::EDFilter {
public:
    explicit EventSelectionFilter(fhicl::ParameterSet const &p);
    EventSelectionFilter(EventSelectionFilter const &) = delete;
    EventSelectionFilter(EventSelectionFilter &&) = delete;
    EventSelectionFilter &operator=(EventSelectionFilter const &) = delete;
    EventSelectionFilter &operator=(EventSelectionFilter &&) = delete;
    void beginSubRun(art::SubRun &subrun) override;
    bool filter(art::Event &e) override;
    bool endSubRun(art::SubRun &subrun) override;
    using ProxyPfpColl_t = common::ProxyPfpColl_t;
    using ProxyPfpElem_t = common::ProxyPfpElem_t;

private:
    art::InputTag fPFPproducer;
    art::InputTag fCLSproducer;
    art::InputTag fSLCproducer;
    art::InputTag fHITproducer;
    art::InputTag fSHRproducer;
    art::InputTag fVTXproducer;
    art::InputTag fTRKproducer;
    art::InputTag fPCAproducer;
    art::InputTag fMCPproducer;
    art::InputTag fWIREproducer;
    art::InputTag fMCTproducer;

    bool _verbose;
    bool _data;
    bool _fake_data;
    bool _numi;
    bool _filter;

    TTree *_tree;
    int _run;
    int _sub;
    int _evt;
    int _selected;

    TTree *_subrun_tree;
    int _run_sr;
    int _sub_sr;
    double _pot; // now always written, for both data & MC

    std::map<unsigned int, unsigned int> _pfpmap;

    std::unique_ptr<::selection::SelectionToolBase> _selectionTool;
    std::vector<std::unique_ptr<::analysis::AnalysisToolBase>> _analysisToolsVec;

    // ---- Normalization / counting config ----
    std::string _swtrig_proc;                   // process name for swtrigger
    std::vector<std::string> _beam_gate_algos;  // OR of these defines "beam gate"
    std::vector<std::string> _ext_gate_algos;   // OR of these defines "EXT gate"
    bool _count_with_opfilter;                  // also require optical filter when counting?
    std::string _opfilter_proc;                 // process name for opfiltercommon
    float _op_pe_beam_min;                      // PE_Beam >= this
    float _op_pe_veto_max;                      // PE_Veto <= this
    art::InputTag _pot_prod_data;               // SubRun POT producer (data)
    art::InputTag _pot_prod_mc;                 // SubRun POT producer (MC)

    // ---- Per-subrun counters / derived ----
    Long64_t _n_beam_gates_sr = 0;              // beam gate count in this subrun
    Long64_t _n_ext_gates_sr  = 0;              // EXT gate count in this subrun
    double _pot_per_gate_sr   = -1.0;           // POT / beam_gates (if both known)
    double _ext_equiv_pot_sr  = -1.0;           // POT_per_gate * n_ext_gates

    void BuildPFPMap(const common::ProxyPfpColl_t &pfp_pxy_col);
    template <typename T>
    void printPFParticleMetadata(const common::ProxyPfpElem_t &pfp_pxy, const T &pfParticleMetadataList);
    void AddDaughters(const common::ProxyPfpElem_t &pfp_pxy, const common::ProxyPfpColl_t &pfp_pxy_col, std::vector<common::ProxyPfpElem_t> &slice_v);
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
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fMCTproducer = p.get<art::InputTag>("MCTproducer");
    fWIREproducer = p.get<art::InputTag>("WIREproducer");

    _verbose = p.get<bool>("Verbose", false);
    _data = p.get<bool>("IsData", false);
    _fake_data = p.get<bool>("IsFakeData", false);
    _numi = p.get<bool>("IsNuMI", false);
    _filter = p.get<bool>("Filter", false);

    // --- new normalization config (with sensible defaults) ---
    _swtrig_proc         = p.get<std::string>("BeamSWTrigProcName", "DataOverlayOptical");
    _beam_gate_algos     = p.get<std::vector<std::string>>("BeamGateAlgoNames", {});
    _ext_gate_algos      = p.get<std::vector<std::string>>("ExtGateAlgoNames", {});
    _count_with_opfilter = p.get<bool>("CountWithOpFilter", false);
    _opfilter_proc       = p.get<std::string>("OpFilterProcName", "DataStage1Optical");
    _op_pe_beam_min      = p.get<float>("OpFilterPEBeamMin", 0.0);
    _op_pe_veto_max      = p.get<float>("OpFilterPEVetoMax", 1e9);
    _pot_prod_data       = art::InputTag(p.get<std::string>("POTProducerData", "beamdata"));
    _pot_prod_mc         = art::InputTag(p.get<std::string>("POTProducerMC",  "generator"));

    // per-subrun counters (explicit init)
    _n_beam_gates_sr = 0;
    _n_ext_gates_sr  = 0;
    _pot_per_gate_sr = -1.0;
    _ext_equiv_pot_sr= -1.0;

    art::ServiceHandle<art::TFileService> tfs;

    _tree = tfs->make<TTree>("EventSelectionFilter", "Neutrino Selection TTree");
    _tree->Branch("selected", &_selected, "selected/I");
    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("sub", &_sub, "sub/I");
    _tree->Branch("evt", &_evt, "evt/I");
    _tree->Branch("is_data", &_data, "is_data/O");
    _tree->Branch("is_fake_data", &_fake_data, "is_fake_data/O");
    _tree->Branch("is_numi", &_numi, "is_numi/O");

    _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
    _subrun_tree->Branch("run", &_run_sr, "run/I");
    _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");
    // Always persist POT (data or MC) plus gate counters and derived quantities
    _subrun_tree->Branch("pot", &_pot, "pot/D");
    _subrun_tree->Branch("n_beam_gates", &_n_beam_gates_sr, "n_beam_gates/L");
    _subrun_tree->Branch("n_ext_gates",  &_n_ext_gates_sr,  "n_ext_gates/L");
    _subrun_tree->Branch("pot_per_gate", &_pot_per_gate_sr, "pot_per_gate/D");
    _subrun_tree->Branch("ext_equiv_pot", &_ext_equiv_pot_sr, "ext_equiv_pot/D");

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

void EventSelectionFilter::beginSubRun(art::SubRun &)
{
    // reset per-subrun counters/derived values
    _n_beam_gates_sr = 0;
    _n_ext_gates_sr  = 0;
    _pot_per_gate_sr = -1.0;
    _ext_equiv_pot_sr= -1.0;
}

bool EventSelectionFilter::filter(art::Event &e) {
    this->ResetTTree();
    if (_verbose) 
        std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl;

    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();

    // -------- Count gates (software trigger; optionally optical filter) --------
    bool pass_beam_gate = false;
    bool pass_ext_gate  = false;
    try {
        art::InputTag trigTag("swtrigger", "", _swtrig_proc);
        const auto& trigH = e.getValidHandle<raw::ubdaqSoftwareTriggerData>(trigTag);
        auto has_any = [&](const std::vector<std::string>& wanted)->bool{
            bool ok = false;
            for (auto const& name : wanted) ok = ok || trigH->passedAlgo(name);
            return ok;
        };
        pass_beam_gate = has_any(_beam_gate_algos);
        pass_ext_gate  = has_any(_ext_gate_algos);
    } catch (...) {
        // missing swtrigger product => leave as false
    }
    bool op_ok = true;
    if (_count_with_opfilter) {
        try {
            art::Handle<uboone::UbooneOpticalFilter> opH;
            art::InputTag opTag("opfiltercommon", "", _opfilter_proc);
            e.getByLabel(opTag, opH);
            if (opH.isValid()) {
                const float pe_beam = opH->PE_Beam();
                const float pe_veto = opH->PE_Veto();
                op_ok = (pe_beam >= _op_pe_beam_min) && (pe_veto <= _op_pe_veto_max);
            }
        } catch (...) {
            // if required but missing, keep op_ok=true to avoid accidental zeroing
            op_ok = true;
        }
    }
    if (pass_beam_gate && op_ok) ++_n_beam_gates_sr;
    if (pass_ext_gate  && op_ok) ++_n_ext_gates_sr;
    // --------------------------------------------------------------------------

    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, fPFPproducer,
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
        _analysisToolsVec[i]->analyseEvent(e, _data);
    }

    bool keepEvent = false;
    for (const common::ProxyPfpElem_t &pfp_pxy : pfp_proxy) {
        if (!pfp_pxy->IsPrimary() || (fabs(pfp_pxy->PdgCode()) != 12 && fabs(pfp_pxy->PdgCode()) != 14)) {
            continue;
        }

        std::vector<common::ProxyPfpElem_t> neutrino_slice;
        this->AddDaughters(pfp_pxy, pfp_proxy, neutrino_slice);

        if (!neutrino_slice.empty()) {
            bool selected = _selectionTool->selectEvent(e, neutrino_slice);
            if (selected) {
                keepEvent = true;
                _selected++;
            }

            for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
                _analysisToolsVec[i]->analyseSlice(e, neutrino_slice, _data, selected);
            }
        }
    }

    _tree->Fill();

    if (_filter)
        return keepEvent;

    return true;
}

void EventSelectionFilter::BuildPFPMap(const common::ProxyPfpColl_t &pfp_pxy_col) {
    _pfpmap.clear();
    unsigned int p = 0;
    for (const auto &pfp_pxy : pfp_pxy_col) {
        _pfpmap[pfp_pxy->Self()] = p;
        p++;
    }
}

template <typename T>
void EventSelectionFilter::printPFParticleMetadata(const common::ProxyPfpElem_t &pfp_pxy, const T &pfParticleMetadataList) {
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
    }
}

void EventSelectionFilter::AddDaughters(const common::ProxyPfpElem_t &pfp_pxy, const common::ProxyPfpColl_t &pfp_pxy_col, std::vector<common::ProxyPfpElem_t> &slice_v) {
    auto daughters = pfp_pxy->Daughters();
    slice_v.push_back(pfp_pxy);
    if (_verbose)
        std::cout << "\t PFP w/ PdgCode " << pfp_pxy->PdgCode() << " has " << daughters.size() << " daughters" << std::endl;
    for (auto const &daughterid : daughters) {
        if (_pfpmap.find(daughterid) == _pfpmap.end())
            continue;
        auto pfp_pxy2 = pfp_pxy_col.begin();
        for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
            ++pfp_pxy2;
        this->AddDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
    }
}

void EventSelectionFilter::ResetTTree() {
    _selected = 0;
    _run = -1;
    _sub = -1;
    _evt = -1;
    _selectionTool->resetTTree(_tree);
    for (size_t i = 0; i < _analysisToolsVec.size(); i++)
        _analysisToolsVec[i]->resetTTree(_tree);
}

bool EventSelectionFilter::endSubRun(art::SubRun &subrun) {
    // Try to read data POT first; if absent, fall back to MC POT
    _pot = -1.0;
    {
        art::Handle<sumdata::POTSummary> h;
        if (subrun.getByLabel(_pot_prod_data, h) && h.isValid()) {
            _pot = static_cast<double>(h->totpot);
        } else if (subrun.getByLabel(_pot_prod_mc, h) && h.isValid()) {
            _pot = static_cast<double>(h->totpot);
        }
    }

    // Derived per-subrun quantities
    _pot_per_gate_sr = (_n_beam_gates_sr > 0 && _pot > 0.0)
                         ? (_pot / static_cast<double>(_n_beam_gates_sr))
                         : -1.0;
    if (_pot_per_gate_sr > 0.0 && _n_ext_gates_sr > 0) {
        _ext_equiv_pot_sr = _pot_per_gate_sr * static_cast<double>(_n_ext_gates_sr);
    } else {
        _ext_equiv_pot_sr = -1.0;
    }
    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();
    return true;
}

DEFINE_ART_MODULE(EventSelectionFilter)
