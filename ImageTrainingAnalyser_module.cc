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

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"

#include "TDatabasePDG.h"

#ifdef ClassDef
#undef ClassDef
#endif
#include <torch/torch.h>
#include <torch/script.h>

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cmath>

#include "ImageProcessor.h"

class ImageTrainingAnalyser : public art::EDAnalyser
{
public:
    explicit ImageTrainingAnalyser(fhicl::ParameterSet const& pset);
    
    ImageTrainingAnalyser(ImageTrainingAnalyser const&) = delete;
    ImageTrainingAnalyser(ImageTrainingAnalyser&&) = delete;
    ImageTrainingAnalyser& operator=(ImageTrainingAnalyser const&) = delete;
    ImageTrainingAnalyser& operator=(ImageTrainingAnalyser&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;

private:
    bool _training_mode;
    bool _veto_bad_channels;

    std::string _training_output_file;
    TFile* _root_file;
    TTree* _image_tree;
    std::unique_ptr<image::ImageManager> _image_manager;

    int _image_width, _image_height;
    float _drift_step;
    float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;

    std::shared_ptr<torch::jit::script::Module> _model_u, _model_v, _model_w;

    std::vector<std::unique_ptr<signature::SignatureToolBase>> _signatureToolsVec; 

    std::string _bad_channel_file; 
    std::vector<bool> _bad_channel_mask;

    const geo::GeometryCore* _geo;
    const detinfo::DetectorProperties* _detp;

    void initialiseBadChannelMask();

    void filterBadChannels(std::vector<art::Ptr<recob::Wire>>& wires);
    void filterBadSimChannels(std::vector<art::Ptr<sim::SimChannel>>& sim_channels);
}

ConvolutionNetworkAlog::ConvolutionNetworkAlog(fhicl::ParamaterSet const& pset)
    : EDAnalyzer{pset}
    , _training_mode{pset.get<bool>("TrainingMode", true)}
    , _veto_bad_channels{pset.get<bool>("VetoBadChannels", true)} 
    , _training_output_file{pset.get<std::string>("TrainingOutputFile", "training_output")}
    , _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")} 
    , _width{pset.get<int>("ImageWidth", 256)}
    , _height{pset.get<int>("ImageHeight", 256)}
    , _WREproducer{pset.get<art::InputTag>("WireProducer", "butcher")}
    , _SCHproducer{pset.get<art::InputTag>("SimChannelProducer", "largeant")}
    , _HITproducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BKTproducer{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _PFPproducer{pset.get<art::InputTag>("PFPproducer", "pandora")}
    , _CLSproducer{pset.get<art::InputTag>("CLSproducer", "pandora")}
    , _SHRproducer{pset.get<art::InputTag>("SHRproducer", "pandora")}
    , _SLCproducer{pset.get<art::InputTag>("SLCproducer", "pandora")}
    , _VTXproducer{pset.get<art::InputTag>("VTXproducer", "pandora")}
    , _PCAproducer{pset.get<art::InputTag>("PCAproducer", "pandora")}
    , _TRKproducer{pset.get<art::InputTag>("TRKproducer", "pandora")}
{
    try {
        if (!_training_mode) 
        {
            _model_u = torch::jit::load(pset.get<std::string>("ModelFileU"));
            _model_v = torch::jit::load(pset.get<std::string>("ModelFileV"));
            _model_w = torch::jit::load(pset.get<std::string>("ModelFileW"));
        }
    } catch (const c10::Error& e) {
        throw cet::exception("ImageTrainingAnalyser") << "**error loading torch models: " << e.what() << "\n";
    }

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }

    _geo = art::ServiceHandle<geo::Geometry>()->provider();

    _drift_step = (_detp->SamplingRate()/1000.) * _detp->DriftVelocity(_detp->Efield(), _detp->Temperature());
    _wire_pitch_u = _geo->WirePitch(geo::kU);                 // U plane
    _wire_pitch_v = _geo->WirePitch(geo::kV);                 // V plane
    _wire_pitch_w = _geo->WirePitch(geo::kW);                 // W plane

    size_t n_channels = _geo->Nchannels();
    _bad_channel_mask.resize(n_channels, false);

    if(_veto_bad_channels)
        this->initialiseBadChannelMask();
}

void ImageTrainingAnalyser::initialiseBadChannelMask()
{
    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) {
            throw cet::exception("ImageTrainingAnalyser") << "**bad channel file not found: " << _bad_channel_file;
        }

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

void ImageTrainingAnalyser::beginJob() 
{
    _root_file = new TFile(_training_output_file.c_str(), "RECREATE");
    _image_tree = new TTree("ImageTree", "Tree containing training images");
    _image_manager = std::make_unique<image::ImageManager>(_image_tree, *_geo);
}

void ImageTrainingAnalyser::analyze(const art::Event* evt) 
{   
    _image_manager->reset();

    if (_training_mode) {
        this->produceTrainingSample(evt);
    } else {
        this->infer(evt);
    }
}

void ImageTrainingAnalyser::filterBadChannels(std::vector<art::Ptr<recob::Wire>>& wires)
{
    wires.erase(std::remove_if(wires.begin(), wires.end(), 
        [this](const art::Ptr<recob::Wire>& wire) { 
            return _bad_channel_mask[wire->Channel()]; 
        }), 
        wires.end());
}

void ImageTrainingAnalyser::filterBadSimChannels(std::vector<art::Ptr<sim::SimChannel>>& sim_channels)
{
    sim_channels.erase(std::remove_if(sim_channels.begin(), sim_channels.end(), 
        [this](const art::Ptr<sim::SimChannel>& sim_channel) { 
            return _bad_channel_mask[sim_channel->Channel()]; 
        }), 
        sim_channels.end());
}

void ImageTrainingAnalyser::produceTrainingSample(const art::Event* evt) 
{
    signature::Pattern pattern;
    bool pattern_found = true;

    for (auto& signatureTool : _signatureToolsVec) {
        signature::Signature signature; 
        if (!signatureTool->constructSignature(evt, signature)) {
            pattern_found = false;
            break;
        }
        pattern.push_back(signature);
    }

    if (!pattern_found && !pattern.empty())
        pattern.clear();

    int run = evt->run();
    int subrun = evt->subRun();
    int event = evt->event();

    auto const& pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(
        *evt, _PFPproducer,
        proxy::withAssociated<recob::Cluster>(_CLSproducer),
        proxy::withAssociated<recob::Slice>(_SLCproducer),
        proxy::withAssociated<recob::Hit>(_HITproducer)
    );

    std::vector<art::Ptr<recob::Hit>> neutrino_hits;
    common::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
        e, _CLSproducer, proxy::withAssociated<recob::Hit>(_CLSproducer)
    );

    for (const auto& pfp : slice_pfp_v) {
        if (pfp->IsPrimary()) continue; 

        auto clus_pxy_v = pfp.get<recob::Cluster>();

        for (auto ass_clus : clus_pxy_v) {
            const auto& clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();

            neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
        }
    }

    double sum_charge_u = 0.0, sum_wire_u = 0.0, sum_drift_u = 0.0;
    double sum_charge_v = 0.0, sum_wire_v = 0.0, sum_drift_v = 0.0;
    double sum_charge_w = 0.0, sum_wire_w = 0.0, sum_drift_w = 0.0;

    for (const auto& hit : neutrino_hits) {
        double charge = hit->Integral();
        common::PandoraView pandora_view = common::GetPandoraView(hit);
        TVector3 hit_pos = common::GetPandoraHitPosition(*evt, hit, pandora_view);

        if (pandora_view == common::TPC_VIEW_U) {
            sum_charge_u += charge;
            sum_wire_u += hit_pos.Z() * charge;
            sum_drift_u += hit_pos.X() * charge;
        } 
        else if (pandora_view == common::TPC_VIEW_V) {
            sum_charge_v += charge;
            sum_wire_v += hit_pos.Z() * charge;
            sum_drift_v += hit_pos.X() * charge;
        } 
        else if (pandora_view == common::TPC_VIEW_W) {
            sum_charge_w += charge;
            sum_wire_w += hit_pos.Z() * charge;
            sum_drift_w += hit_pos.X() * charge;
        }
    }

    double centroid_wire_u = (sum_charge_u > 0) ? sum_wire_u / sum_charge_u : 0.0;
    double centroid_drift_u = (sum_charge_u > 0) ? sum_drift_u / sum_charge_u : 0.0;

    double centroid_wire_v = (sum_charge_v > 0) ? sum_wire_v / sum_charge_v : 0.0;
    double centroid_drift_v = (sum_charge_v > 0) ? sum_drift_v / sum_charge_v : 0.0;

    double centroid_wire_w = (sum_charge_w > 0) ? sum_wire_w / sum_charge_w : 0.0;
    double centroid_drift_w = (sum_charge_w > 0) ? sum_drift_w / sum_charge_w : 0.0;

    std::vector<image::ImageMeta> metas;
    metas.emplace_back(
        centroid_wire_u - (_image_width / 2.0) * _wire_pitch_u,
        centroid_drift_u - (_image_height / 2.0) * _drift_step,
        _image_height, _image_width, _wire_pitch_u, _drift_step, geo::kU
    );
    metas.emplace_back(
        centroid_wire_v - (_image_width / 2.0) * _wire_pitch_v,
        centroid_drift_v - (_image_height / 2.0) * _drift_step,
        _image_height, _image_width, _wire_pitch_v, _drift_step, geo::kV
    );
    metas.emplace_back(
        centroid_wire_w - (_image_width / 2.0) * _wire_pitch_w,
        centroid_drift_w - (_image_height / 2.0) * _drift_step,
        _image_height, _image_width, _wire_pitch_w, _drift_step, geo::kW
    );

    _image_manager->add(run, subrun, event, pattern_found, neutrino_hits, metas);
}
