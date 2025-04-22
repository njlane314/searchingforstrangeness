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
#include "CommonDefs/Image.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "CommonDefs/Pandora.h"
#include "CommonDefs/Types.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "CommonDefs/ReconstructionLabels.h"
#include "CommonDefs/PrimaryLabels.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>

class EventSelectionFilter : public art::EDFilter {
public:
    explicit EventSelectionFilter(fhicl::ParameterSet const &p);
    EventSelectionFilter(EventSelectionFilter const &) = delete;
    EventSelectionFilter(EventSelectionFilter &&) = delete;
    EventSelectionFilter &operator=(EventSelectionFilter const &) = delete;
    EventSelectionFilter &operator=(EventSelectionFilter &&) = delete;
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
    art::InputTag fDeadChannelTag;
    std::string fBackTrackerLabel;  
    double fGammaThreshold;         
    double fHadronThreshold;        

    bool _verbose;
    bool _data;
    bool _fake_data;
    bool _filter;

    TTree *_tree;
    int _run;
    int _sub;
    int _evt;
    int _selected;

    TTree *_subrun_tree;
    int _run_sr;
    int _sub_sr;
    float _pot;

    image::ImageProperties image_properties_;
    float _adc_image_threshold;

    std::map<unsigned int, unsigned int> _pfpmap;

    std::unique_ptr<::selection::SelectionToolBase> _selectionTool;
    std::vector<std::unique_ptr<::analysis::AnalysisToolBase>> _analysisToolsVec;

    const geo::GeometryCore* _geo;
    detinfo::DetectorPropertiesData _detp;

    float _drift_step;
    float _wire_pitch_u;
    float _wire_pitch_v;
    float _wire_pitch_w;
    int _image_width;
    int _image_height;

    void BuildPFPMap(const common::ProxyPfpColl_t &pfp_pxy_col);
    template <typename T>
    void printPFParticleMetadata(const common::ProxyPfpElem_t &pfp_pxy, const T &pfParticleMetadataList);
    void AddDaughters(const common::ProxyPfpElem_t &pfp_pxy, const common::ProxyPfpColl_t &pfp_pxy_col, std::vector<common::ProxyPfpElem_t> &slice_v);
    void ResetTTree();

    std::vector<common::ProxyPfpElem_t> collectNeutrinoSlice(const common::ProxyPfpColl_t& pfp_proxy);

    std::vector<art::Ptr<recob::Hit>> collectNeutrinoHits(const art::Event& e, 
                                                        const std::vector<common::ProxyPfpElem_t>& neutrino_slice);

    std::pair<double, double> calculateCentroid(const art::Event& e, 
                                                common::PandoraView view, 
                                                const std::vector<art::Ptr<recob::Hit>>& hits);

    void constructImages(const art::Event& e,
                        const std::vector<image::ImageProperties>& properties,
                        std::vector<image::Image>& calo_images,
                        std::vector<image::Image>& reco_images,
                        std::vector<image::Image>& label_images);
};

EventSelectionFilter::EventSelectionFilter(fhicl::ParameterSet const &p)
    : EDFilter{p},
      _detp(art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(
          art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob())) {
    // Retrieve parameters from the parameter set
    fPFPproducer = p.get<art::InputTag>("PFPproducer");
    fSHRproducer = p.get<art::InputTag>("SHRproducer");
    fHITproducer = p.get<art::InputTag>("HITproducer");
    fCLSproducer = p.get<art::InputTag>("CLSproducer");
    fSLCproducer = p.get<art::InputTag>("SLCproducer");
    fVTXproducer = p.get<art::InputTag>("VTXproducer");
    fPCAproducer = p.get<art::InputTag>("PCAproducer");
    fTRKproducer = p.get<art::InputTag>("TRKproducer");
    fMCPproducer = p.get<art::InputTag>("MCPproducer");
    fWIREproducer = p.get<art::InputTag>("WIREproducer");
    fDeadChannelTag = p.get<art::InputTag>("DeadChannelTag", "nfbadchannel:badchannels:OverlayDetsim");
    fBackTrackerLabel = p.get<std::string>("BackTrackerLabel", "gaushit");
    fGammaThreshold = p.get<double>("GammaThreshold", 0.1);         
    fHadronThreshold = p.get<double>("HadronThreshold", 0.1);

    _verbose = p.get<bool>("Verbose");
    _data = p.get<bool>("IsData");
    _fake_data = p.get<bool>("IsFakeData", false);
    _filter = p.get<bool>("Filter", false);
    _adc_image_threshold = p.get<float>("ADCthreshold", 10.0);

    _image_width = p.get<int>("ImageWidth", 512);
    _image_height = p.get<int>("ImageHeight", 512);

    _geo = art::ServiceHandle<geo::Geometry>()->provider();

    auto const& channelMap = art::ServiceHandle<geo::WireReadout>()->Get();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

    double tick_period = clockData.TPCClock().TickPeriod();
    double drift_velocity = _detp.DriftVelocity();
    _drift_step = tick_period * drift_velocity * 1e1;

    _wire_pitch_u = channelMap.Plane(geo::PlaneID(0, 0, geo::kU)).WirePitch();
    _wire_pitch_v = channelMap.Plane(geo::PlaneID(0, 0, geo::kV)).WirePitch();
    _wire_pitch_w = channelMap.Plane(geo::PlaneID(0, 0, geo::kW)).WirePitch();

    art::ServiceHandle<art::TFileService> tfs;

    _tree = tfs->make<TTree>("EventSelectionFilter", "Neutrino Selection TTree");
    _tree->Branch("selected", &_selected, "selected/I");
    _tree->Branch("run", &_run, "run/I");
    _tree->Branch("sub", &_sub, "sub/I");
    _tree->Branch("evt", &_evt, "evt/I");

    _subrun_tree = tfs->make<TTree>("SubRun", "SubRun TTree");
    _subrun_tree->Branch("run", &_run_sr, "run/I");
    _subrun_tree->Branch("subRun", &_sub_sr, "subRun/I");

    if ((!_data) || (_fake_data))
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
    if (_verbose) 
        std::cout << "new event : [run,event] : [" << e.run() << ", " << e.event() << "]" << std::endl;

    _evt = e.event();
    _sub = e.subRun();
    _run = e.run();

    std::vector<art::Ptr<recob::Wire>> wire_vec;
    if (auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer))
        art::fill_ptr_vector(wire_vec, wireHandle);

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
    auto neutrino_slice = this->collectNeutrinoSlice(pfp_proxy);
    if (!neutrino_slice.empty()) {
        auto neutrino_hits = this->collectNeutrinoHits(e, neutrino_slice);
        
        std::vector<image::ImageProperties> properties;
        auto [centroid_wire_u, centroid_drift_u] = this->calculateCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateCentroid(e, common::TPC_VIEW_W, neutrino_hits);

        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_height, _image_width, _wire_pitch_u, _drift_step, geo::kU);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_height, _image_width, _wire_pitch_v, _drift_step, geo::kV);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_height, _image_width, _wire_pitch_w, _drift_step, geo::kW);

        std::vector<image::Image> calo_images, reco_images, label_images;
        this->constructImages(e, properties, calo_images, reco_images, label_images);

        bool selected = _selectionTool->selectEvent(e, neutrino_slice, calo_images, reco_images, label_images);
        if (selected) {
            keepEvent = true;
            _selected = 1;
        }

        for (size_t i = 0; i < _analysisToolsVec.size(); i++) {
            _analysisToolsVec[i]->analyseSlice(e, neutrino_slice, _data, selected);
        }
    }

    _tree->Fill();
    if (_filter)
        return keepEvent;

    return true;
}

void EventSelectionFilter::constructImages(const art::Event& e,
                                            const std::vector<image::ImageProperties>& properties,
                                            std::vector<image::Image>& calo_images,
                                            std::vector<image::Image>& reco_images,
                                            std::vector<image::Image>& label_images) {
    calo_images.clear();
    reco_images.clear();
    label_images.clear();

    for (const auto& prop : properties) {
        image::Image calo_image(prop);
        calo_image.clear(static_cast<float>(0.0));
        calo_images.push_back(std::move(calo_image));

        image::Image reco_image(prop);
        reco_image.clear(static_cast<float>(reco_labels::ReconstructionLabel::empty));
        reco_images.push_back(std::move(reco_image));

        image::Image label_image(prop);
        label_image.clear(static_cast<float>(truth_labels::PrimaryLabel::empty));
        label_images.push_back(std::move(label_image));
    }

    auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
    auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
    auto mcpHandle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCPproducer);

    art::FindManyP<recob::Hit> wire_hit_assoc(wireHandle, e, fHITproducer);
    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hitHandle, e, fBackTrackerLabel);

    std::vector<truth_labels::PrimaryLabel> primary_labels = truth_labels::classifyParticles(e, fMCPproducer, fGammaThreshold, fHadronThreshold);
    std::vector<reco_labels::ReconstructionLabel> reco_labels = reco_labels::classifyParticles(e, fMCPproducer, fGammaThreshold, fHadronThreshold);

    std::map<int, size_t> trackid_to_index;
    for (size_t i = 0; i < mcpHandle->size(); ++i) {
        trackid_to_index[mcpHandle->at(i).TrackId()] = i;
    }

    auto const det_props = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
    auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

    for (size_t wire_idx = 0; wire_idx < wireHandle->size(); ++wire_idx) {
        const auto& wire = wireHandle->at(wire_idx);
        auto ch_id = wire.Channel();
        //if (ch_id < _bad_channel_mask.size() && _bad_channel_mask[ch_id]) continue;

        std::vector<geo::WireID> wire_ids = wireReadout.ChannelToWire(ch_id);
        if (wire_ids.empty()) continue;
        geo::View_t view = wireReadout.Plane(wire_ids.front().planeID()).View();
        size_t view_idx = view - geo::kU;

        geo::Point_t center = wireReadout.Wire(wire_ids.front()).GetCenter();
        TVector3 wire_center(center.X(), center.Y(), center.Z());
        double wire_coord = (view == geo::kW) ? wire_center.Z() :
                            (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                                (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));

        const auto& hits = wire_hit_assoc.at(wire_idx);

        for (const auto& range : wire.SignalROI().get_ranges()) {
            const auto& adcs = range.data();
            int start_tick = range.begin_index();
            for (size_t idx = 0; idx < adcs.size(); ++idx) {
                int tick = start_tick + idx;
                double x = det_props.ConvertTicksToX(tick, wire_ids.front().planeID());
                size_t row = properties[view_idx].row(x);
                size_t col = properties[view_idx].col(wire_coord);
                if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;

                reco_labels::ReconstructionLabel reco_label = reco_labels::ReconstructionLabel::cosmic;
                truth_labels::PrimaryLabel primary_label = truth_labels::PrimaryLabel::cosmic;

                for (const auto& hit : hits) {
                    if (tick >= hit->StartTick() && tick < hit->EndTick()) {
                        auto bkth_data = mcp_bkth_assoc.data(hit.key());
                        if (!bkth_data.empty()) {
                            for (size_t i = 0; i < bkth_data.size(); ++i) {
                                if (bkth_data[i]->isMaxIDE == 1) {
                                    int track_id = mcp_bkth_assoc.at(hit.key())[i]->TrackId();
                                    auto it = trackid_to_index.find(track_id);
                                    if (it != trackid_to_index.end()) {
                                        size_t particle_idx = it->second;
                                        if (particle_idx < primary_labels.size()) {
                                            primary_label = primary_labels[particle_idx];
                                        }
                                        if (particle_idx < reco_labels.size()) {
                                            reco_label = reco_labels[particle_idx];
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                        break;
                    }
                }

                if (adcs[idx] > _adc_image_threshold) {
                    calo_images[view_idx].set(row, col, adcs[idx]);
                    reco_images[view_idx].set(row, col, static_cast<float>(reco_label), false);
                    label_images[view_idx].set(row, col, static_cast<float>(primary_label), false);
                }
            }
        }
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

void EventSelectionFilter::BuildPFPMap(const common::ProxyPfpColl_t &pfp_pxy_col) {
    _pfpmap.clear();
    unsigned int p = 0;
    for (const auto &pfp_pxy : pfp_pxy_col) {
        _pfpmap[pfp_pxy->Self()] = p;
        p++;
    }
}

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
        _pot = subrun.getByLabel(fMCPproducer, potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    }

    _run_sr = subrun.run();
    _sub_sr = subrun.subRun();
    _subrun_tree->Fill();

    return true;
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

std::vector<common::ProxyPfpElem_t> EventSelectionFilter::collectNeutrinoSlice(const common::ProxyPfpColl_t& pfp_proxy) {
    std::vector<common::ProxyPfpElem_t> neutrino_slice;
    for (const common::ProxyPfpElem_t& pfp_pxy : pfp_proxy) {
        if (pfp_pxy->IsPrimary() && (fabs(pfp_pxy->PdgCode()) == 12 || fabs(pfp_pxy->PdgCode()) == 14)) {
            this->AddDaughters(pfp_pxy, pfp_proxy, neutrino_slice);
        }
    }

    return neutrino_slice;
}

std::vector<art::Ptr<recob::Hit>> EventSelectionFilter::collectNeutrinoHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& neutrino_slice) {
    std::vector<art::Ptr<recob::Hit>> neutrino_hits;
    auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(e, fCLSproducer, proxy::withAssociated<recob::Hit>(fCLSproducer));
    
    for (const auto& pfp : neutrino_slice) {
        if (pfp->IsPrimary()) continue;
        for (auto ass_clus : pfp.get<recob::Cluster>()) {
            auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
            neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
        }
    }

    return neutrino_hits;
}

std::pair<double, double> EventSelectionFilter::calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) {
    double sum_charge = 0.0;
    double sum_wire = 0.0;
    double sum_drift = 0.0;

    for (const auto& hit : hits) {
        if (common::GetPandoraView(hit) != view) continue;
        double charge = hit->Integral();
        TVector3 hit_pos = common::GetPandoraHitPosition(e, hit, view);
        sum_charge += charge;
        sum_wire += hit_pos.Z() * charge;
        sum_drift += hit_pos.X() * charge;
    }

    if (sum_charge == 0.0) return {0.0, 0.0};

    return {sum_wire / sum_charge, sum_drift / sum_charge};
}

DEFINE_ART_MODULE(EventSelectionFilter)