#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
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

    void infer(art::Event const& evt, std::map<int, std::vector<art::Ptr<recob::Hit>>>& classified_hits);

private:
    bool _training_mode;
    int _pass;
    
    std::string _training_output_file;
    std::shared_ptr<torch::jit::script::Module> _model_u, _model_v, _model_w;

    int _width, _height;

    float _drift_step;
    float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;
    std::map<common::PandoraView, float> _wire_pitch;

    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer;

    std::map<common::PandoraView, std::array<float, 4>> _region_bounds;
    std::vector<art::Ptr<recob::Hit>> _region_hits;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

    calo::CalorimetryAlg* _calo_alg;

    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;

    std::string _bad_channel_file;
    bool _veto_bad_channels;
    std::vector<bool> _bad_channel_mask;
    const geo::GeometryCore* _geo;

    void initialiseEvent(art::Event const& evt);
    void initialiseBadChannelMask();
    void prepareTrainingSample(art::Event const& evt);
    void produceTrainingSample(const std::string& filename, const std::vector<float>& feat_vec, bool result);
    void makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, const common::PandoraView view, torch::Tensor& network_input, std::map<art::Ptr<recob::Hit>,std::pair<int, int>>& calohit_pixel_map);
    void findRegionBounds(art::Event const& evt, const std::vector<art::Ptr<recob::Hit>>& hits);
    void getNuVertex(art::Event const& evt, std::array<float, 3>& nu_vtx, bool& found_vertex);
    void calculateChargeCentroid(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hits, std::map<common::PandoraView, std::array<float, 2>>& q_cent_map, std::map<common::PandoraView, float>& tot_q_map);
    std::tuple<float, float, float, float> getBoundsForView(common::PandoraView view) const;
};

ConvolutionNetworkAlgo::ConvolutionNetworkAlgo(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset}
    , _training_mode{pset.get<bool>("TrainingMode", true)}
    , _pass{pset.get<int>("Pass", 1)}
    , _training_output_file{pset.get<std::string>("TrainingOutputFile", "training_output")}
    , _width{pset.get<int>("ImageWidth", 256)}
    , _height{pset.get<int>("ImageHeight", 256)}
    , _drift_step{pset.get<float>("DriftStep", 0.5)}
    , _wire_pitch_u{pset.get<float>("WirePitchU", 0.3)}
    , _wire_pitch_v{pset.get<float>("WirePitchU", 0.3)}
    , _wire_pitch_w{pset.get<float>("WirePitchU", 0.3)}
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
            std::cout << "In testing mode!" << std::endl;
            std::cout << pset.get<std::string>("ModelFileU") << std::endl;
            _model_u = torch::jit::load(pset.get<std::string>("ModelFileU"));
            _model_v = torch::jit::load(pset.get<std::string>("ModelFileV"));
            _model_w = torch::jit::load(pset.get<std::string>("ModelFileW"));
            std::cout << "Loaded models" << std::endl;
        }
    } catch (const c10::Error& e) {
        throw cet::exception("ConvolutionNetworkAlgo") << "Error loading Torch models: " << e.what() << "\n";
    }

    _calo_alg = new calo::CalorimetryAlg(pset.get<fhicl::ParameterSet>("CaloAlg"));

    _wire_pitch = {
        {common::TPC_VIEW_U, _wire_pitch_u},
        {common::TPC_VIEW_V, _wire_pitch_v},
        {common::TPC_VIEW_W, _wire_pitch_w}
    };

    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }

    _geo = art::ServiceHandle<geo::Geometry>()->provider();
    size_t num_channels = _geo->Nchannels();
    _bad_channel_mask.resize(num_channels, false);

    if (_veto_bad_channels)
        this->initialiseBadChannelMask();
}

void ConvolutionNetworkAlgo::analyze(art::Event const& evt) 
{   
    this->initialiseEvent(evt); 
    if (_region_hits.empty())
        return;

    try {
        if (_training_mode)
            this->prepareTrainingSample(evt);
    } catch (const c10::Error& e) {
        throw cet::exception("ConvolutionNetworkAlgo") << "Error running algorithm: " << e.what() << "\n";
    }
}

void ConvolutionNetworkAlgo::initialiseEvent(art::Event const& evt)
{
    _region_bounds.clear();
    _region_hits.clear(); 
    _mcp_bkth_assoc.reset();

    std::vector<art::Ptr<recob::Hit>> evt_hits, all_hits, sim_hits;
    art::Handle<std::vector<recob::Hit>> hit_handle;
    
    if (evt.getByLabel(_HitProducer, hit_handle))
    {
        art::fill_ptr_vector(evt_hits, hit_handle);
        _mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, evt, _BacktrackTag);

        for (const auto& hit : evt_hits) 
        {
            if (_veto_bad_channels && _bad_channel_mask[hit->Channel()]) 
                continue;
            
            all_hits.push_back(hit);

            auto assmcp = _mcp_bkth_assoc->at(hit.key());
            auto assmdt = _mcp_bkth_assoc->data(hit.key());
            for (unsigned int ia = 0; ia < assmcp.size(); ++ia)
            {
                auto amd = assmdt[ia];
                if (amd->isMaxIDE != 1)
                    continue;
                
                sim_hits.push_back(hit);
            }
        }
    }

    if (sim_hits.empty() || all_hits.empty()) 
        return;

    mf::LogInfo("ConvolutionNetworkAlgo") << "Input Hit size: " << sim_hits.size();

    this->findRegionBounds(evt, sim_hits);
    if (_region_bounds.empty())
        return;

    for (const auto& hit : sim_hits)
    {
        common::PandoraView view = common::GetPandoraView(hit);
        auto [drift_min, drift_max, wire_min, wire_max] = this->getBoundsForView(view);

        const auto pos = common::GetPandoraHitPosition(evt, hit, view);
        float x = pos.X();
        float z = pos.Z();

        if (x >= drift_min && x <= drift_max && z >= wire_min && z <= wire_max)
            _region_hits.push_back(hit);
    }

    mf::LogInfo("ConvolutionNetworkAlgo") << "Region Hit size: " << _region_hits.size();
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

void ConvolutionNetworkAlgo::findRegionBounds(art::Event const& evt, const std::vector<art::Ptr<recob::Hit>>& hits)
{
    std::map<common::PandoraView, std::array<float, 2>> q_cent_map;
    std::map<common::PandoraView, float> tot_q_map;
    common::initialiseChargeMap(q_cent_map, tot_q_map);
    this->calculateChargeCentroid(evt, hits, q_cent_map, tot_q_map);

    for (const auto& view : {common::TPC_VIEW_U, common::TPC_VIEW_V, common::TPC_VIEW_W}) 
    {
        const auto [x_centroid, z_centroid] = q_cent_map[view];

        float x_min = x_centroid - (_height / 2) * _drift_step;
        float x_max = x_centroid + (_height / 2) * _drift_step;
        float z_min = z_centroid - (_width / 2) * _wire_pitch[view];
        float z_max = z_centroid + (_width / 2) * _wire_pitch[view];

        _region_bounds[view] = {x_min, x_max, z_min, z_max};

        std::cout << "View: " 
            << (view == common::TPC_VIEW_U ? "U" : (view == common::TPC_VIEW_V ? "V" : "W")) 
            << ", X bounds: [" << x_min << ", " << x_max << "]"
            << ", Z bounds: [" << z_min << ", " << z_max << "]" << std::endl;
    }
}

std::tuple<float, float, float, float> ConvolutionNetworkAlgo::getBoundsForView(common::PandoraView view) const
{
    const auto& bounds = _region_bounds.at(view);  
    float drift_min = bounds[0]; 
    float drift_max = bounds[1]; 
    float wire_min = bounds[2];   
    float wire_max = bounds[3];   
    return std::make_tuple(drift_min, drift_max, wire_min, wire_max);
}

void ConvolutionNetworkAlgo::calculateChargeCentroid(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hits, std::map<common::PandoraView, std::array<float, 2>>& q_cent_map, std::map<common::PandoraView, float>& tot_q_map)
{
    for (const auto& hit : hits)
    {
        common::PandoraView view = common::GetPandoraView(hit);
        const TVector3 pos = common::GetPandoraHitPosition(evt, hit, view);
        float charge = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);

        q_cent_map[view][0] += pos.X() * charge;  
        q_cent_map[view][1] += pos.Z() * charge;  
        tot_q_map[view] += charge;
    }

    for (auto& [view, charge_center] : q_cent_map)
    {
        if (tot_q_map[view] > 0) 
        {
            charge_center[0] /= tot_q_map[view];
            charge_center[1] /= tot_q_map[view];
        }
    }
}

void ConvolutionNetworkAlgo::prepareTrainingSample(art::Event const& evt) 
{
    std::array<float, 3> nu_vtx = {0.0f, 0.0f, 0.0f};
    bool found_vertex = false;

    this->getNuVertex(evt, nu_vtx, found_vertex);
    if (!found_vertex) 
        return; 

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

    unsigned int n_flags = _signatureToolsVec.size(); 
    int run = evt.run();
    int subrun = evt.subRun();
    int event = evt.event();

    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> region_hits;
    for (const art::Ptr<recob::Hit>& hit : _region_hits) 
    {
        common::PandoraView view = common::GetPandoraView(hit);
        region_hits[view].push_back(hit);
    }

    for (const auto& [view, evt_view_hits] : region_hits)
    {
        float x_vtx = nu_vtx[0];
        float z_vtx = (common::ProjectToWireView(nu_vtx[0], nu_vtx[1], nu_vtx[2], view)).Z();

        auto [drift_min, drift_max, wire_min, wire_max] = this->getBoundsForView(view);
        if (x_vtx > (drift_min - 1.f) && x_vtx < (drift_max + 1.f) && z_vtx > (wire_min - 1.f) && z_vtx < (wire_max + 1.f))
        {           
            unsigned int n_hits = 0;
            unsigned int n_meta = 0;
            std::vector<float> feat_vec = { static_cast<float>(n_hits), 
                                            static_cast<float>(n_flags),
                                            static_cast<float>(n_meta),
                                            static_cast<float>(run), 
                                            static_cast<float>(subrun), 
                                            static_cast<float>(event),
                                            static_cast<float>(_height),
                                            static_cast<float>(_width),
                                            x_vtx, z_vtx, 
                                            drift_min, drift_max, 
                                            wire_min, wire_max };

            n_meta = feat_vec.size();
            feat_vec[2] = static_cast<float>(n_meta);

            for (const auto& hit : evt_view_hits)
            {
                const geo::WireID hit_wire(hit->WireID());
                if (hit_wire.Wire >= art::ServiceHandle<geo::Geometry>()->Nwires(hit_wire)) 
                    continue;

                const auto pos = common::GetPandoraHitPosition(evt, hit, static_cast<common::PandoraView>(view));
                float x = pos.X();
                float z = pos.Z();
                float q = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);

                std::vector<float> signature_flags(n_flags, 0.f);
                if (_mcp_bkth_assoc != nullptr) 
                {
                    const auto& assmcp = _mcp_bkth_assoc->at(hit.key());
                    const auto& assmdt = _mcp_bkth_assoc->data(hit.key());

                    for (unsigned int ia = 0; ia < assmcp.size(); ++ia) 
                    {
                        bool found_flag = false;
                        if (assmdt[ia]->isMaxIDE == 1) 
                        {
                            size_t sig_ctr = 0;
                            for (const auto& sig : patt) 
                            {
                                for (size_t it = 0; it < sig.size(); ++it)
                                {
                                    if (sig[it]->TrackId() == assmcp[ia]->TrackId()) 
                                    {
                                        signature_flags.at(sig_ctr) = 1.f;
                                        found_flag = true;

                                        break;
                                    }
                                }
                                sig_ctr++;
                            }
                        }

                        if (found_flag)
                            break;
                    }
                }

                feat_vec.insert(feat_vec.end(), {x, z, q});
                feat_vec.insert(feat_vec.end(), signature_flags.begin(), signature_flags.end());
                ++n_hits;
            }

            feat_vec[0] = static_cast<float>(n_hits);

            std::string view_string = (view == common::TPC_VIEW_U) ? "U" : (view == common::TPC_VIEW_V) ? "V" : "W";
            std::string training_filename = _training_output_file + "_" + view_string + ".csv";
            this->produceTrainingSample(training_filename, feat_vec, true);
        }
    }
}

void ConvolutionNetworkAlgo::getNuVertex(art::Event const& evt, std::array<float, 3>& nu_vtx, bool& found_vertex)
{
    found_vertex = false;

    auto const &mct_h = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    const simb::MCTruth& mct = mct_h->at(0);
    if (!mct.NeutrinoSet()) return;

    auto const &neutrino = mct.GetNeutrino();
    auto const &nu = neutrino.Nu();

    common::True2RecoMappingXYZ(nu.T(), nu.Vx(), nu.Vy(), nu.Vz(), nu_vtx.data());
    found_vertex = true;
}

void ConvolutionNetworkAlgo::produceTrainingSample(const std::string& filename, const std::vector<float>& feat_vec, bool result)
{
    std::ofstream out_file(filename, std::ios_base::app);
    if (!out_file.is_open()) {
        mf::LogError("ConvolutionNetworkAlgo") << "Error: Could not open file " << filename;
        return;
    }

    std::string delimiter = ",";

    for (const float &feature : feat_vec)
        out_file << feature << delimiter;

    out_file << static_cast<int>(result) << '\n';

    out_file.close();
}

void ConvolutionNetworkAlgo::infer(art::Event const& evt, std::map<int, std::vector<art::Ptr<recob::Hit>>>& classified_hits) 
{
    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> region_hits;
    for (const auto& hit : _region_hits)
        region_hits[common::GetPandoraView(hit)].push_back(hit);

    for (const auto& [view, evt_view_hits] : region_hits)
    {
        torch::Tensor network_input;
        std::map<art::Ptr<recob::Hit>,std::pair<int, int>> calohit_pixel;

        this->makeNetworkInput(evt, evt_view_hits, view, network_input, calohit_pixel);

        torch::Tensor output;
        if (view == common::TPC_VIEW_U)
            output = _model_u->forward({network_input}).toTensor();
        else if (view == common::TPC_VIEW_V)
            output = _model_v->forward({network_input}).toTensor();
        else if (view == common::TPC_VIEW_W)
            output = _model_w->forward({network_input}).toTensor();

        torch::Tensor predicted_classes = torch::argmax(output, 1); 
        auto classes_accessor{predicted_classes.accessor<int64_t, 4>()};
        
        for (size_t i = 0; i < evt_view_hits.size(); ++i)
        {
            const auto& hit = evt_view_hits[i];

            const auto pixel = calohit_pixel.at(hit);
            int predicted_class = predicted_classes[0][pixel.first][pixel.second].item<int>();

            classified_hits[predicted_class].push_back(hit);
        }
    }

    for (const auto& [class_id, hits] : classified_hits)
        std::cout << "Class " << class_id << " has " << hits.size() << " hits.";
}

void ConvolutionNetworkAlgo::makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, const common::PandoraView view, torch::Tensor& network_input, std::map<art::Ptr<recob::Hit>,std::pair<int, int>>& calohit_pixel_map)
{
    const auto [x_min, x_max, z_min, z_max] = this->getBoundsForView(view);
    std::vector<double> x_bin_edges(_width + 1);
    std::vector<double> z_bin_edges(_height + 1);

    x_bin_edges[0] = x_min;
    const double dx = (x_max - x_min) / _width;

    for (int i = 1; i <= _width; ++i)
        x_bin_edges[i] = x_bin_edges[i - 1] + dx;

    z_bin_edges[0] = z_min;
    const double dz = (z_max - z_min) / _height;

    for (int i = 1; i <= _height; ++i)
        z_bin_edges[i] = z_bin_edges[i - 1] + dz;

    network_input = torch::zeros({1, 1, _height, _width});
    auto accessor = network_input.accessor<float, 4>();
    for (const auto& hit : hit_list)
    {
        const auto pos = common::GetPandoraHitPosition(evt, hit, static_cast<common::PandoraView>(view));
        float x = pos.X();
        float z = pos.Z();

        const int pixel_x{static_cast<int>(std::floor((x - x_bin_edges[0]) / dx))};
        const int pixel_z{static_cast<int>(std::floor((z - z_bin_edges[0]) / dz))};

        if (pixel_x >= 0 && pixel_x < _width && pixel_z >= 0 && pixel_z < _height)
        {
            float q = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);
            accessor[0][0][pixel_z][pixel_x] += q;
            calohit_pixel_map.insert({hit, {pixel_z, pixel_x}});
        }
    }
}

void ConvolutionNetworkAlgo::beginJob() 
{}

void ConvolutionNetworkAlgo::endJob() 
{}

DEFINE_ART_MODULE(ConvolutionNetworkAlgo)
