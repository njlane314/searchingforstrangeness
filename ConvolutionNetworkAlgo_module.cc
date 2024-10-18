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

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"

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

    art::InputTag _HitProducer, _MCPproducer, _MCTproducer, _BacktrackTag;
    std::vector<art::Ptr<recob::Hit>> _event_hits;

    calo::CalorimetryAlg* _calo_alg;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

    std::vector<std::unique_ptr<::signature::SignatureToolBase>> _signatureToolsVec;

    void initialiseEvent(art::Event const& evt);
    void prepareTrainingSample(art::Event const& evt);
    void produceTrainingSample(const std::string& filename, const std::vector<float>& feat_vec, bool result);
    void makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, common::PandoraView view, float x_min, float x_max, float z_min, float z_max, torch::Tensor& network_input, std::map<art::Ptr<recob::Hit>,std::pair<int, int>>& m_calohit_pixel);
    void findRegionExtent(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, float& x_min, float& x_max, float& z_min, float& z_max) const;
    void getNuVertex(art::Event const& evt, std::array<float, 3>& nu_vtx, bool& found_vertex);

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
{
    try {
        if (!_training_mode) {
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
}

void ConvolutionNetworkAlgo::analyze(art::Event const& evt) 
{   
    this->prepareTrainingSample(evt);
}

void ConvolutionNetworkAlgo::initialiseEvent(art::Event const& evt)
{
    _event_hits.clear(); 
    _mcp_bkth_assoc.reset();

    art::Handle<std::vector<recob::Hit>> hit_handle;
    if (evt.getByLabel(_HitProducer, hit_handle))
    {
        art::fill_ptr_vector(_event_hits, hit_handle);
        _mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, evt, _BacktrackTag);

        mf::LogInfo("ConvNetworkAlgorithm") << "Input Hit size: " << _event_hits.size();
    }
}

void ConvolutionNetworkAlgo::prepareTrainingSample(art::Event const& evt) 
{
    std::cout << "Preparing training sample..." << std::endl;
    this->initialiseEvent(evt);

    std::array<float, 3> nu_vtx = {0.0f, 0.0f, 0.0f};
    bool found_vertex = false;

    this->getNuVertex(evt, nu_vtx, found_vertex);
    if (!found_vertex) 
        return; 

    std::vector<signature::Trace> trace_coll;
    bool found_all_signatures = true;
    for (auto& signatureTool : _signatureToolsVec)
    {
        bool found_signature = signatureTool->identifySignalParticles(evt, trace_coll);

        if (!found_signature) 
            found_all_signatures = false;
    }

    if (!found_all_signatures) 
    {
        std::cout << "Didn't find all the signatures!" << std::endl;
        return;
    }

    std::cout << "size of trace collection " << trace_coll.size() << std::endl;

    for (auto& trace : trace_coll)
    {
        std::cout << "Trace: " << trace.pdg << ", " << trace.trckid << std::endl;
    }

    unsigned int n_flags = trace_coll.size();

    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> intr_hits;
    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> evt_hits;

    for (const art::Ptr<recob::Hit>& hit : _event_hits) 
    {
        common::PandoraView view = common::GetPandoraView(hit);
        evt_hits[view].push_back(hit);

        const std::vector<art::Ptr<simb::MCParticle>>& mc_particle_vector = _mcp_bkth_assoc->at(hit.key());
        const auto& matching_data_vector = _mcp_bkth_assoc->data(hit.key());

        for (size_t i_p = 0; i_p < mc_particle_vector.size(); ++i_p) 
        {
            if (matching_data_vector[i_p]->isMaxIDE == 1) 
            {
                intr_hits[view].push_back(hit);
                break;  
            }
        }
    }

    for (const auto& [view, evt_view_hits] : evt_hits)
    {
        float intr_drft_min, intr_drft_max, intr_wire_min, intr_wire_max;
        this->findRegionExtent(evt, intr_hits[view], intr_drft_min, intr_drft_max, intr_wire_min, intr_wire_max);

        float evt_drft_min, evt_drft_max, evt_wire_min, evt_wire_max;
        this->findRegionExtent(evt, evt_view_hits, evt_drft_min, evt_drft_max, evt_wire_min, evt_wire_max);

        float x_vtx = nu_vtx[0];
        float z_vtx = (common::ProjectToWireView(nu_vtx[0], nu_vtx[1], nu_vtx[2], view)).Z();

        if (x_vtx > (evt_drft_min - 1.f) && x_vtx < (evt_drft_max + 1.f) && z_vtx > (evt_wire_min - 1.f) && z_vtx < (evt_wire_max + 1.f))
        {            
            unsigned int n_hits = 0;
            // also store the event, subrun and run numbers to the vector
            std::vector<float> feat_vec = { static_cast<float>(n_hits), 
                                            static_cast<float>(n_flags),
                                            x_vtx, z_vtx, 
                                            evt_drft_min, evt_drft_max, 
                                            evt_wire_min, evt_wire_max, 
                                            intr_drft_min, intr_drft_max,
                                            intr_wire_min, intr_wire_max };

            for (unsigned int ih = 0; ih < evt_view_hits.size(); ih++)
            {   
                const auto hit = evt_view_hits.at(ih);
                const geo::WireID hit_wire(hit->WireID());

                if (hit_wire.Wire >= art::ServiceHandle<geo::Geometry>()->Nwires(hit_wire)) 
                    continue;

                const auto hit_pos = common::GetPandoraHitPosition(evt, hit, static_cast<common::PandoraView>(view));
                float x = hit_pos.X();
                float z = hit_pos.Z();
                float q = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);

                feat_vec.insert(feat_vec.end(), {x, z, q});

                std::vector<float> sig_flags(trace_coll.size(), 0.0f);
                if (_mcp_bkth_assoc != nullptr) 
                {
                    const auto& assmcp = _mcp_bkth_assoc->at(hit.key());
                    const auto& assmdt = _mcp_bkth_assoc->data(hit.key());

                    if (!assmcp.empty() && !assmdt.empty()) 
                    {
                        for (unsigned int ia = 0; ia < assmcp.size(); ++ia) 
                        {   
                            auto mcp = assmcp[ia]; 
                            auto amd = assmdt[ia];

                            if (amd->isMaxIDE == 1) 
                            {
                                for (size_t it = 0; it < trace_coll.size(); ++it) 
                                {
                                    if (mcp->TrackId() == trace_coll[it].trckid) 
                                        sig_flags[it] = 1.0f;
                                }
                            }
                        }
                    }
                }

                ++n_hits;
                feat_vec.insert(feat_vec.end(), sig_flags.begin(), sig_flags.end());
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

void ConvolutionNetworkAlgo::findRegionExtent(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>> &hit_v, float &x_min, float &x_max, float &z_min, float &z_max) const
{   
    x_min = std::numeric_limits<float>::max();
    x_max = -std::numeric_limits<float>::max();
    z_min = std::numeric_limits<float>::max();
    z_max = -std::numeric_limits<float>::max();

    for (const art::Ptr<recob::Hit> &hit : hit_v)
    {
        common::PandoraView view = common::GetPandoraView(hit);
        TVector3 pos = common::GetPandoraHitPosition(evt, hit, view);

        x_min = std::min<float>(x_min, static_cast<float>(pos.X()));
        x_max = std::max<float>(x_max, static_cast<float>(pos.X()));
        z_min = std::min<float>(z_min, static_cast<float>(pos.Z()));
        z_max = std::max<float>(z_max, static_cast<float>(pos.Z()));
    }
}

void ConvolutionNetworkAlgo::produceTrainingSample(const std::string& filename, const std::vector<float>& feat_vec, bool result)
{
    std::cout << "Attempting to open file: " << filename << std::endl;
    std::ofstream out_file(filename, std::ios_base::app);
    if (!out_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string delimiter = ",";

    for (const float &feature : feat_vec)
        out_file << feature << delimiter;

    out_file << static_cast<int>(result) << '\n';

    std::cout << "Saving the training output!" << std::endl;
    out_file.close();
}

void ConvolutionNetworkAlgo::infer(art::Event const& evt, std::map<int, std::vector<art::Ptr<recob::Hit>>>& classified_hits) 
{
    std::cout << "Startig inference..." << std::endl;
    this->initialiseEvent(evt);

    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> evt_hits;
    for (const auto& hit : _event_hits)
        evt_hits[common::GetPandoraView(hit)].push_back(hit);

    for (const auto& [view, evt_view_hits] : evt_hits)
    {
        float drift_min = std::numeric_limits<float>::max();
        float drift_max = -std::numeric_limits<float>::max();
        float wire_min, wire_max;
        this->findRegionExtent(evt, evt_view_hits, drift_min, drift_max, wire_min, wire_max);

        torch::Tensor network_input;
        std::map<art::Ptr<recob::Hit>,std::pair<int, int>> calohit_pixel;

        this->makeNetworkInput(evt, evt_view_hits, view, drift_min, drift_max, wire_min, wire_max, network_input, calohit_pixel);
       
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
    {
        std::cout << "Class " << class_id << " has " << hits.size() << " hits." << std::endl;
    }

    std::cout << "Ending inference!" << std::endl;
}

void ConvolutionNetworkAlgo::makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, const common::PandoraView view, const float x_min, const float x_max, const float z_min, const float z_max, torch::Tensor& network_input, std::map<art::Ptr<recob::Hit>,std::pair<int, int>>& m_calohit_pixel)
{
    const float pitch = _wire_pitch[view];
    std::vector<double> x_bin_edges(_width + 1);
    std::vector<double> z_bin_edges(_height + 1);
    x_bin_edges[0] = x_min - 0.5f * _drift_step; 
    const double dx = ((x_max + 0.5f * _drift_step) - x_bin_edges[0]) / _width; 

    for (int i = 1; i <= _width; ++i) 
        x_bin_edges[i] = x_bin_edges[i - 1] + dx;

    z_bin_edges[0] = z_min - 0.5f * pitch; 
    const double dz = ((z_max + 0.5f * pitch) - z_bin_edges[0]) / _height;

    for (int i = 1; i <= _height; ++i) 
        z_bin_edges[i] = z_bin_edges[i - 1] + dz;

    network_input = torch::zeros({1, 1, _height, _width});
    auto accessor = network_input.accessor<float, 4>();

    for (const auto& hit : hit_list)
    {
        const auto hit_pos = common::GetPandoraHitPosition(evt, hit, static_cast<common::PandoraView>(view));
        float x = hit_pos.X();
        float z = hit_pos.Z();

        if (_pass > 1 && (x < x_min || x > x_max || z > z_max))
            continue;

        float q = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);
        const int pixel_x{static_cast<int>(std::floor((x - x_bin_edges[0]) / dx))};
        const int pixel_z{static_cast<int>(std::floor((z - z_bin_edges[0]) / dz))};
        accessor[0][0][pixel_z][pixel_x] += q;
        m_calohit_pixel.insert({hit,{pixel_z,pixel_x}});
    }
}

void ConvolutionNetworkAlgo::beginJob() 
{}

void ConvolutionNetworkAlgo::endJob() 
{}

DEFINE_ART_MODULE(ConvolutionNetworkAlgo)
