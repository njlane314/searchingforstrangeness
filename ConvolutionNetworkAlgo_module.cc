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

class ConvolutionNetworkAlgo : public art::EDAnalyzer 
{
public:
    explicit ConvolutionNetworkAlgo(fhicl::ParameterSet const& pset);

    ConvolutionNetworkAlgo(ConvolutionNetworkAlgo const&) = delete;
    ConvolutionNetworkAlgo(ConvolutionNetworkAlgo&&) = delete;
    ConvolutionNetworkAlgo& operator=(ConvolutionNetworkAlgo const&) = delete;
    ConvolutionNetworkAlgo& operator=(ConvolutionNetworkAlgo&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;

    void infer(const art::Event* evt, std::map<int, std::vector<art::Ptr<recob::Hit>>>& classified_hits);

private:
    bool _training_mode;
    
    std::string _training_output_file;
    TFile* _root_file;
    TTree* _image_tree;
    std::unique_ptr<image::ImageManager> _image_manager;

    std::shared_ptr<torch::jit::script::Module> _model_u, _model_v, _model_w;

    int _width, _height;

    float _drift_step;
    float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;

    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer;

    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;

    std::string _bad_channel_file;
    bool _veto_bad_channels;
    std::vector<bool> _bad_channel_mask;
    const geo::GeometryCore* _geo;
    
    std::vector<art::Ptr<recob::Wire>> _wires;

    void initialiseEvent(const art::Event* evt);
    void initialiseBadChannelMask();
    void prepareTrainingSample(const art::Event& evt);
    void produceTrainingSample(const art::Event& evt, const std::vector<recob::Wire>& wires, const std::vector<sim::SimChannel>& channels);
    void makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, const common::PandoraView view, torch::Tensor& network_input, std::map<art::Ptr<recob::Hit>,std::pair<int, int>>& calohit_pixel_map);
    void getNuVertex(const art::Event* evt, std::array<float, 3>& nu_vtx, bool& found_vertex);
    void calculateChargeCentroid(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hits, std::map<common::PandoraView, std::array<float, 2>>& q_cent_map, std::map<common::PandoraView, float>& tot_q_map);
};

ConvolutionNetworkAlgo::ConvolutionNetworkAlgo(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset}
    , _training_mode{pset.get<bool>("TrainingMode", true)}
    , _training_output_file{pset.get<std::string>("TrainingOutputFile", "training_output")}
    , _width{pset.get<int>("ImageWidth", 256)}
    , _height{pset.get<int>("ImageHeight", 256)}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _PFPproducer{pset.get<art::InputTag>("PFPproducer", "pandora")}
    , _CLSproducer{pset.get<art::InputTag>("CLSproducer", "pandora")}
    , _SHRproducer{pset.get<art::InputTag>("SHRproducer", "pandora")}
    , _SLCproducer{pset.get<art::InputTag>("SLCproducer", "pandora")}
    , _VTXproducer{pset.get<art::InputTag>("VTXproducer", "pandora")}
    , _PCAproducer{pset.get<art::InputTag>("PCAproducer", "pandora")}
    , _TRKproducer{pset.get<art::InputTag>("TRKproducer", "pandora")}
    , _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")}
    , _veto_bad_channels{pset.get<bool>("VetoBadChannels", true)}
{
    try {
        if (!_training_mode) 
        {
            _model_u = torch::jit::load(pset.get<std::string>("ModelFileU"));
            _model_v = torch::jit::load(pset.get<std::string>("ModelFileV"));
            _model_w = torch::jit::load(pset.get<std::string>("ModelFileW"));
        }
    } catch (const c10::Error& e) {
        throw cet::exception("ConvolutionNetworkAlgo") << "Error loading Torch models: " << e.what() << "\n";
    }

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }

    _geo = art::ServiceHandle<geo::Geometry>()->provider();
    
    _drift_step = 0.5; 
    _wire_pitch_u = _geo->WirePitch(geo::kU);                 // U plane
    _wire_pitch_v = _geo->WirePitch(geo::kV);                 // V plane
    _wire_pitch_w = _geo->WirePitch(geo::kW);                 // W plane

    size_t num_channels = _geo->Nchannels();
    _bad_channel_mask.resize(num_channels, false);

    if (_veto_bad_channels)
        this->initialiseBadChannelMask();
}

void ConvolutionNetworkAlgo::beginJob() 
{
    _root_file = new TFile(_training_output_file.c_str(), "RECREATE");
    _image_tree = new TTree("ImageTree", "Tree containing training images");
    _image_manager = std::make_unique<image::ImageManager>(_image_tree, *_geo);
}

void ConvolutionNetworkAlgo::analyze(const art::Event* evt) 
{   
    this->initialiseEvent(evt); 
    if (_training_mode) {
        this->prepareTrainingSample(evt);
    } else {
        this->infer(evt);
    }
}

void ConvolutionNetworkAlgo::initialiseEvent(const art::Event* evt)
{
    _image_manager->reset();

    art::Handle<std::vector<recob::Wire>> wire_handle;
    if(!evt.getByLabel(_WREproducer))
        return;

    
    art::fill_ptr_vector(_wires, wire_handle);
}

void ConvolutionNetworkAlgo::initialiseBadChannelMask()
{
    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) {
            throw cet::exception("ConvolutionNetworkAlgo") << "Bad channel file not found: " << _bad_channel_file;
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
        std::cout << "Loaded bad channels from: " << fullname << std::endl;
    }
}

void ConvolutionNetworkAlgo::calculateChargeCentroid(const std::vector<recob::Wire>& wires, double& centroid_x, double& centroid_y)
{
    double total_charge = 0.0;
    double weighted_x = 0.0, weighted_y = 0.0;

    for (const auto& wire : wires) {
        const geo::WireID& wire_id = wire.WireID();
        TVector3 wire_center = _geo->WireIDToWireGeo(wire_id).GetCenter();
        double charge = std::accumulate(wire.SignalROI().waveform().begin(), wire.SignalROI().waveform().end(), 0.0);

        weighted_x += wire_center.X() * charge;
        weighted_y += wire_center.Y() * charge;
        total_charge += charge;
    }

    if (total_charge > 0) {
        centroid_x = weighted_x / total_charge;
        centroid_y = weighted_y / total_charge;
    }
}

void ConvolutionNetworkAlgo::filterBadChannels(std::vector<recob::Wire>& wires)
{
    wires.erase(std::remove_if(wires.begin(), wires.end(), 
        [this](const recob::Wire& wire) { return _bad_channel_mask[wire.Channel()]; }), 
        wires.end());
}

void ConvolutionNetworkAlgo::prepareTrainingSample(const art::Event* evt) 
{
    signature::Pattern patt;
    bool patt_found = true;
    for (auto& signatureTool : _signatureToolsVec) {
        signature::Signature signature; 
        if (!signatureTool->constructSignature(evt, signature)) {
            patt_found = false;
            break;
        }

        patt.push_back(signature);
    }

    if (!patt_found && !patt.empty())
        patt.clear();

    int run = evt.run();
    int subrun = evt.subRun();
    int event = evt.event();

    
         
}

void ConvolutionNetworkAlgo::produceTrainingSample(const std::string& filename, const std::vector<float>& feat_vec, bool result)
{
}

void ConvolutionNetworkAlgo::infer(const art::Event* evt, std::map<int, std::vector<art::Ptr<recob::Hit>>>& classified_hits) 
{
}

void ConvolutionNetworkAlgo::makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, const common::PandoraView view, torch::Tensor& network_input, std::map<art::Ptr<recob::Hit>,std::pair<int, int>>& calohit_pixel_map)
{
}

void ConvolutionNetworkAlgo::endJob() 
{}

DEFINE_ART_MODULE(ConvolutionNetworkAlgo)
