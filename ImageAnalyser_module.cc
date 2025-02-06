#ifdef ClassDef
#undef ClassDef
#endif
#include <torch/torch.h>
#include <torch/script.h>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"
#include "CommonFunctions/Geometry.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"
#include "EventClassifier.h"

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TTree.h"

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <algorithm>

#include "ImageProcessor.h"

class ImageAnalyser : public art::EDAnalyzer
{
public:
    explicit ImageAnalyser(fhicl::ParameterSet const& pset);
    
    ImageAnalyser(ImageAnalyser const&) = delete;
    ImageAnalyser(ImageAnalyser&&) = delete;
    ImageAnalyser& operator=(ImageAnalyser const&) = delete;
    ImageAnalyser& operator=(ImageAnalyser&&) = delete;

    void analyze(art::Event const& e) override;
    void beginSubRun(art::SubRun const& sbr) override;

    void beginJob() override;
    void endJob() override;

    using ProxyPfpColl_t = common::ProxyPfpColl_t;
    using ProxyPfpElem_t = common::ProxyPfpElem_t;

private:
    fhicl::ParameterSet _pset;

    std::string _training_output_file;
    TFile* _root_file;
    TTree* _job_tree;
    TTree* _image_tree;
    std::unique_ptr<image::ImageTrainingHandler> _image_handler;

    int _image_width, _image_height;
    float _drift_step;
    float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;

    int _kernel_size;

    std::string _bad_channel_file; 
    std::vector<bool> _bad_channel_mask;

    const geo::GeometryCore* _geo;
    const detinfo::DetectorProperties* _detp;

    std::map<int, unsigned int> _pfpmap;

    double _job_pot;

    art::InputTag _WREproducer, _POTlabel, _HITproducer, _MCPproducer, _MCTproducer, _BKTproducer, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer;

    void initialiseBadChannelMask();

    void filterBadChannels(std::vector<art::Ptr<recob::Wire>>& wires);

    void produceTrainingSample(art::Event const& e);
    void addDaughters(const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v);
    void buildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);

    bool constructSignatures(const art::Event& e, signature::Pattern& pattern);
    std::vector<common::ProxyPfpElem_t> collectNeutrinoSlice(const ProxyPfpColl_t& pfp_proxy);
    std::vector<art::Ptr<recob::Hit>> collectNeutrinoHits(const art::Event& e, const std::vector<ProxyPfpElem_t>& neutrino_slice);
    std::pair<double, double> calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
};

ImageAnalyser::ImageAnalyser(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset}
    , _pset(pset)
    , _training_output_file{pset.get<std::string>("TrainingOutputFile", "training_output.root")}
    , _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")} 
    , _image_width{pset.get<int>("ImageWidth", 512)}
    , _image_height{pset.get<int>("ImageHeight", 512)}
    , _kernel_size{pset.get<int>("KernelSize", 5)}
    , _WREproducer{pset.get<art::InputTag>("WREproducer", "butcher")}
    , _HITproducer{pset.get<art::InputTag>("HITpoducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BKTproducer{pset.get<art::InputTag>("BKTproducer", "gaushitTruthMatch")}
    , _PFPproducer{pset.get<art::InputTag>("PFPproducer", "pandora")}
    , _CLSproducer{pset.get<art::InputTag>("CLSproducer", "pandora")}
    , _SHRproducer{pset.get<art::InputTag>("SHRproducer", "pandora")}
    , _SLCproducer{pset.get<art::InputTag>("SLCproducer", "pandora")}
    , _VTXproducer{pset.get<art::InputTag>("VTXproducer", "pandora")}
    , _PCAproducer{pset.get<art::InputTag>("PCAproducer", "pandora")}
    , _TRKproducer{pset.get<art::InputTag>("TRKproducer", "pandora")}
    , _POTlabel{pset.get<art::InputTag>("POTlabel", "generator")}
{
    _geo = art::ServiceHandle<geo::Geometry>()->provider();

    _drift_step = 0.5;
    _wire_pitch_u = _geo->WirePitch(geo::kU);                 // U plane
    _wire_pitch_v = _geo->WirePitch(geo::kV);                 // V plane
    _wire_pitch_w = _geo->WirePitch(geo::kW);                 // W plane

    size_t n_channels = _geo->Nchannels();
    _bad_channel_mask.resize(n_channels, false);
    this->initialiseBadChannelMask();
}

void ImageAnalyser::analyze(const art::Event& e) 
{   
    _image_handler->reset();
    this->produceTrainingSample(e);
}

void ImageAnalyser::beginJob() 
{
    _root_file = new TFile(_training_output_file.c_str(), "RECREATE");
    _job_tree = new TTree("JobTree", "Tree containing job-level information");
    _image_tree = new TTree("ImageTree", "Tree containing training images");

    _job_pot = 0.0;
    _job_tree->Branch("job_pot", &_job_pot);

    _image_handler = std::make_unique<image::ImageTrainingHandler>(_image_tree, *_geo);
}

void ImageAnalyser::beginSubRun(art::SubRun const& sbr) 
{  
    if (const auto potHandle = sbr.getValidHandle<sumdata::POTSummary>(_POTlabel))
        _job_pot += potHandle->totgoodpot;  
}

void ImageAnalyser::endJob() 
{
    if (_root_file) {
        _root_file->cd();    
        _job_tree->Write();     
        _image_tree->Write();     
        _root_file->Close();      

        delete _root_file;        
        _root_file = nullptr;     
    }
}

void ImageAnalyser::produceTrainingSample(const art::Event& e) 
{
    std::vector<art::Ptr<recob::Wire>> wire_vec;
    if (auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(_WREproducer)) 
        art::fill_ptr_vector(wire_vec, wireHandle);
    else 
        return;

    std::vector<art::Ptr<recob::Hit>> hit_vec;
    if (auto hitHandle = e.getValidHandle<std::vector<recob::Hit>>(_HITproducer))
        art::fill_ptr_vector(hit_vec, hitHandle);
    else 
        return;

    this->filterBadChannels(wire_vec);

    auto pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, _PFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(_PFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(_CLSproducer),
                                                        proxy::withAssociated<recob::Slice>(_SLCproducer),
                                                        proxy::withAssociated<recob::Track>(_TRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(_VTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(_PCAproducer),
                                                        proxy::withAssociated<recob::Shower>(_SHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(_PFPproducer));

    this->buildPFPMap(pfp_proxy);

    auto neutrino_slice = this->collectNeutrinoSlice(pfp_proxy);
    if (neutrino_slice.empty())
        return;

    auto neutrino_hits = this->collectNeutrinoHits(e, neutrino_slice);
    if (neutrino_hits.size() < 200) 
        return;

    std::vector<image::ImageProperties> properties;

    auto [centroid_wire_u, centroid_drift_u] = this->calculateCentroid(e, common::TPC_VIEW_U, neutrino_hits);
    auto [centroid_wire_v, centroid_drift_v] = this->calculateCentroid(e, common::TPC_VIEW_V, neutrino_hits);
    auto [centroid_wire_w, centroid_drift_w] = this->calculateCentroid(e, common::TPC_VIEW_W, neutrino_hits);

    properties.emplace_back(centroid_wire_u, centroid_drift_u,
        _image_height, _image_width, _wire_pitch_u, _drift_step, geo::kU
    );
    properties.emplace_back(centroid_wire_v, centroid_drift_v,
        _image_height, _image_width, _wire_pitch_v, _drift_step, geo::kV
    );
    properties.emplace_back(centroid_wire_w, centroid_drift_w,
        _image_height, _image_width, _wire_pitch_w, _drift_step, geo::kW
    );

    signature::EventClassifier event_classifier(_pset);
    signature::EventType event_type = event_classifier.classifyEvent(e);
    signature::Pattern pattern = event_classifier.getPattern(e);

    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hit_vec, e, _BKTproducer);

    _image_handler->add(e, event_type, wire_vec, hit_vec, mcp_bkth_assoc, pattern, _kernel_size, properties);
}

void ImageAnalyser::initialiseBadChannelMask()
{
    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) 
            throw cet::exception("ImageAnalyser") << "Bad channel file not found: " << _bad_channel_file;

        std::ifstream inFile(fullname, std::ios::in);
        std::string line;
        while (std::getline(inFile, line)) {
            if (line.find("#") != std::string::npos) continue;
            std::istringstream ss(line);
            int ch1, ch2;
            ss >> ch1;
            if (!(ss >> ch2)) ch2 = ch1;
            for (int i = ch1; i <= ch2; ++i) {
                _bad_channel_mask[i] = true;
            }
        }
    }
}

void ImageAnalyser::filterBadChannels(std::vector<art::Ptr<recob::Wire>>& wires)
{
    wires.erase(std::remove_if(wires.begin(), wires.end(), 
        [this](const art::Ptr<recob::Wire>& wire) { 
            return _bad_channel_mask[wire->Channel()]; 
        }), 
        wires.end());
}

void ImageAnalyser::buildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
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

void ImageAnalyser::addDaughters(const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v)
{
    auto daughters = pfp_pxy->Daughters();
    slice_v.push_back(pfp_pxy);

    for (auto const &daughterid : daughters)
    {
        if (_pfpmap.find(daughterid) == _pfpmap.end())
            continue;

        auto pfp_pxy2 = pfp_pxy_col.begin();
        for (size_t j = 0; j < _pfpmap.at(daughterid); ++j)
            ++pfp_pxy2;

        this->addDaughters(*pfp_pxy2, pfp_pxy_col, slice_v);
    } 
    return;
} 

std::vector<common::ProxyPfpElem_t> ImageAnalyser::collectNeutrinoSlice(const ProxyPfpColl_t& pfp_proxy) {
    std::vector<ProxyPfpElem_t> neutrino_slice;
    
    for (const ProxyPfpElem_t& pfp_pxy : pfp_proxy) {
        if (pfp_pxy->IsPrimary() && (fabs(pfp_pxy->PdgCode()) == 12 || fabs(pfp_pxy->PdgCode()) == 14)) {
            this->addDaughters(pfp_pxy, pfp_proxy, neutrino_slice);
        }
    }

    return neutrino_slice;
}

std::vector<art::Ptr<recob::Hit>> ImageAnalyser::collectNeutrinoHits(const art::Event& e, const std::vector<ProxyPfpElem_t>& neutrino_slice) 
{
    std::vector<art::Ptr<recob::Hit>> neutrino_hits;
    
    auto clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
        e, _CLSproducer, proxy::withAssociated<recob::Hit>(_CLSproducer)
    );

    for (const auto& pfp : neutrino_slice) {
        if (pfp->IsPrimary()) continue;

        for (auto ass_clus : pfp.get<recob::Cluster>()) {
            auto clus_hit_v = clus_proxy[ass_clus.key()].get<recob::Hit>();
            neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
        }
    }

    return neutrino_hits;
}

std::pair<double, double> ImageAnalyser::calculateCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) 
{
    double sum_charge = 0.0, sum_wire = 0.0, sum_drift = 0.0;

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

DEFINE_ART_MODULE(ImageAnalyser)