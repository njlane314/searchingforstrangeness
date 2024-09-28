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

class ConvNetworkAlgorithm : public art::EDAnalyzer {
public:
    explicit ConvNetworkAlgorithm(fhicl::ParameterSet const& pset);

    ConvNetworkAlgorithm(ConvNetworkAlgorithm const&) = delete;
    ConvNetworkAlgorithm(ConvNetworkAlgorithm&&) = delete;
    ConvNetworkAlgorithm& operator=(ConvNetworkAlgorithm const&) = delete;
    ConvNetworkAlgorithm& operator=(ConvNetworkAlgorithm&&) = delete;

    void analyze(art::Event const& e) override;
    void beginJob() override;
    void endJob() override;

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
    std::vector<art::Ptr<recob::Hit>> _hit_list;

    calo::CalorimetryAlg* _calo_alg;
    std::unique_ptr<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>> _mcp_bkth_assoc;

    float _PionThreshold, _MuonThreshold;

    TParticlePDG *neutral_kaon = TDatabasePDG::Instance()->GetParticle(311);
    TParticlePDG *kaon_short = TDatabasePDG::Instance()->GetParticle(310);
    TParticlePDG *kaon_long = TDatabasePDG::Instance()->GetParticle(130);
    TParticlePDG *lambda = TDatabasePDG::Instance()->GetParticle(3122);
    TParticlePDG *sigma_plus = TDatabasePDG::Instance()->GetParticle(3222); 
    TParticlePDG *sigma_minus = TDatabasePDG::Instance()->GetParticle(3112);
    TParticlePDG *sigma_zero = TDatabasePDG::Instance()->GetParticle(3212);
    TParticlePDG *muon = TDatabasePDG::Instance()->GetParticle(13);
    TParticlePDG *pion = TDatabasePDG::Instance()->GetParticle(211);

    void prepareTrainingSample(art::Event const& evt);
    void infer(art::Event const& evt);
    void produceTrainingSample(const std::string& filename, const std::vector<float>& feature_vector);
    void makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, common::PandoraView view, float x_min, float x_max, float z_min, float z_max, torch::Tensor& network_input, std::vector<std::pair<int, int>>& pixel_vector);
    void getHitRegion(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, float& x_min, float& x_max, float& z_min, float& z_max) const;
    void identifySignalParticles(art::Event const& evt, int& muon_tid, int& piplus_tid, int& piminus_tid, bool& found_signature);
};

ConvNetworkAlgorithm::ConvNetworkAlgorithm(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset}
    , _training_mode{pset.get<bool>("TrainingMode", true)}
    , _pass{pset.get<int>("Pass", 1)}
    , _training_output_file{pset.get<std::string>("TrainingOutputFile", "training_output")}
    , _width{pset.get<int>("ImageWidth", 256)}
    , _height{pset.get<int>("ImageHeight", 256)}
    , _drift_step{pset.get<float>("DriftStep", 0.5)}
    , _wire_pitch_u{pset.get<float>("WirePitchU", 0.46669998765)}
    , _wire_pitch_v{pset.get<float>("WirePitchU", 0.46669998765)}
    , _wire_pitch_w{pset.get<float>("WirePitchU", 0.46669998765)}
    , _HitProducer{pset.get<art::InputTag>("HitProducer", "gaushit")}
    , _MCPproducer{pset.get<art::InputTag>("MCPproducer", "largeant")}
    , _MCTproducer{pset.get<art::InputTag>("MCTproducer", "generator")}
    , _BacktrackTag{pset.get<art::InputTag>("BacktrackTag", "gaushitTruthMatch")}
    , _PionThreshold{pset.get<float>("PionThreshold", 0.1)}
    , _MuonThreshold{pset.get<float>("MuonThreshold", 0.1)}
{
    if (!_training_mode)
    {
        _model_u = torch::jit::load(pset.get<std::string>("ModelFileU"));
        _model_v = torch::jit::load(pset.get<std::string>("ModelFileV"));
        _model_w = torch::jit::load(pset.get<std::string>("ModelFileW"));
    }
   
    _calo_alg = new calo::CalorimetryAlg(pset.get<fhicl::ParameterSet>("CaloAlg"));

    _wire_pitch = {
        {common::TPC_VIEW_U, _wire_pitch_u},
        {common::TPC_VIEW_V, _wire_pitch_v},
        {common::TPC_VIEW_W, _wire_pitch_w}
    };
}

void ConvNetworkAlgorithm::analyze(art::Event const& evt) 
{   
    art::Handle<std::vector<recob::Hit>> hit_handle;
    if (evt.getByLabel(_HitProducer, hit_handle))
    {
        art::fill_ptr_vector(_hit_list, hit_handle);
        _mcp_bkth_assoc = std::make_unique<art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData>>(hit_handle, evt, _BacktrackTag);
    }
  
    if (_training_mode)
        this->prepareTrainingSample(evt);
    else
        this->infer(evt);
}

void ConvNetworkAlgorithm::prepareTrainingSample(art::Event const& evt) 
{
    int muon_tid, piplus_tid, piminus_tid;
    bool found_signature;

    this->identifySignalParticles(evt, muon_tid, piplus_tid, piminus_tid, found_signature);
    if (!found_signature) 
        return;

    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> hit_view_map;
    for (const auto& hit : _hit_list)
        hit_view_map[common::GetPandoraView(hit)].push_back(hit);

    for (const auto& [view, hits] : hit_view_map)
    {
        float drift_min = std::numeric_limits<float>::max();
        float drift_max = -std::numeric_limits<float>::max();
        float wire_min, wire_max;

        this->getHitRegion(evt, hits, drift_min, drift_max, wire_min, wire_max);

        std::vector<float> feature_vector = {
            drift_min, drift_max, wire_min, wire_max
        };

        unsigned int n_hits = 0;
        std::unordered_map<int, int> particle_to_flag = {
            {muon_tid, 0},    
            {piplus_tid, 1},  
            {piminus_tid, 2}  
        };

        for (unsigned int ih = 0; ih < hits.size(); ih++)
        {
            const auto hit = hits.at(ih);
            const auto hit_pos = common::GetPandoraHitPosition(evt, hit, static_cast<common::PandoraView>(view));
            float x = hit_pos.X();
            float z = hit_pos.Z();
            float q = _calo_alg->ElectronsFromADCArea(hit->Integral(), hit->WireID().Plane);

            feature_vector.insert(feature_vector.end(), {x, z, q});

            std::array<float, 3> particle_flags = {0.0f, 0.0f, 0.0f};
            auto const& assmcp = _mcp_bkth_assoc->at(ih);
            auto const& assmdt = _mcp_bkth_assoc->data(ih);

            for (unsigned int ia = 0; ia < assmcp.size(); ++ia) 
            {
                auto mcp = assmcp[ia]; 
                auto amd = assmdt[ia];

                if (amd->isMaxIDE != 1) 
                    continue;

                if (particle_to_flag.count(mcp->TrackId())) 
                    particle_flags[particle_to_flag[mcp->TrackId()]] = 1.0f;
            }

            feature_vector.insert(feature_vector.end(), particle_flags.begin(), particle_flags.end());
            ++n_hits;
        }

        if (n_hits < 10) 
            continue;

        std::string training_filename = _training_output_file + "_" + std::to_string(view) + ".csv";
        this->produceTrainingSample(training_filename, feature_vector);
    }
}


void ConvNetworkAlgorithm::identifySignalParticles(art::Event const& evt, int &muon_tid, int &piplus_tid, int &piminus_tid, bool &found_signature)
{
    muon_tid = piplus_tid = piminus_tid = -1;
    found_signature = false;

    auto const &mct_h = evt.getValidHandle<std::vector<simb::MCTruth>>(_MCTproducer);
    auto const &mcp_h = evt.getValidHandle<std::vector<simb::MCParticle>>(_MCPproducer);

    std::map<int, art::Ptr<simb::MCParticle>> mcp_map;
    for (size_t d = 0; d < mcp_h->size(); d++)
    {
        const art::Ptr<simb::MCParticle> mcp(mcp_h, d);
        mcp_map[mcp->TrackId()] = mcp;
    }

    const simb::MCTruth& mct = mct_h->at(0);
    if (!mct.NeutrinoSet()) return;

    auto const &neutrino = mct.GetNeutrino();
    if (neutrino.CCNC() != simb::kCC) return;
    const simb::MCParticle& lepton = neutrino.Lepton();

    if (abs(lepton.PdgCode()) != muon->PdgCode() || lepton.Momentum().Vect().Mag() < _MuonThreshold) return;
    
    muon_tid = lepton.TrackId();

    for (const auto &t_part : *mcp_h) 
    {
        if (abs(t_part.PdgCode()) == neutral_kaon->PdgCode() && t_part.Process() == "primary" && t_part.EndProcess() == "Decay" && t_part.NumberDaughters() == 1 && !found_signature) 
        {
            std::vector<art::Ptr<simb::MCParticle>> dtrs = common::GetDaughters(mcp_map.at(t_part.TrackId()), mcp_map);
            if (dtrs.size() != 1) continue; 

            auto g_part = dtrs.at(0);
            if (g_part->PdgCode() == kaon_short->PdgCode() && g_part->Process() == "Decay" && g_part->EndProcess() == "Decay" && g_part->NumberDaughters() == 2 && !found_signature)
            {
                auto daughters = common::GetDaughters(mcp_map.at(g_part->TrackId()), mcp_map);
                if (daughters.size() != 2) continue;

                std::vector<int> exp_dtrs = {-211, 211};
                std::vector<int> fnd_dtrs;

                for (const auto &dtr : daughters) 
                    fnd_dtrs.push_back(dtr->PdgCode());
                    
                std::sort(exp_dtrs.begin(), exp_dtrs.end());
                std::sort(fnd_dtrs.begin(), fnd_dtrs.end());
                if (fnd_dtrs == exp_dtrs) 
                {
                    for (const auto &dtr : daughters) 
                    {
                        if (dtr->PdgCode() == 211) // pion-plus
                            piplus_tid = dtr->TrackId();
                        
                        else if (dtr->PdgCode() == -211) // pion-minus
                            piminus_tid = dtr->TrackId();
                    }

                    found_signature = std::all_of(daughters.begin(), daughters.end(), [&](const auto& dtr) {
                        return dtr->Momentum().Vect().Mag() >= _PionThreshold;
                    });

                    if (found_signature) break;
                }
            }
        }
    }
}

void ConvNetworkAlgorithm::getHitRegion(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>> &hit_v, float &x_min, float &x_max, float &z_min, float &z_max) const
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

void ConvNetworkAlgorithm::produceTrainingSample(const std::string& filename, const std::vector<float>& feature_vector)
{
    std::ofstream out_file(_training_output_file, std::ios_base::app);

    if (!out_file.is_open())
        return;

    for (const float &feature : feature_vector)
        out_file << feature << ",";
    
    //out_file << static_cast<int>(result) << "\n";

    out_file.close();
}

void ConvNetworkAlgorithm::infer(art::Event const& evt) 
{
    std::vector<art::Ptr<recob::Hit>> muon_hits;
    std::vector<art::Ptr<recob::Hit>> pion_plus_hits;
    std::vector<art::Ptr<recob::Hit>> pion_minus_hits;

    std::map<common::PandoraView, std::vector<art::Ptr<recob::Hit>>> hit_view_map;
    for (const auto& hit : _hit_list)
        hit_view_map[common::GetPandoraView(hit)].push_back(hit);

    for (const auto& [view, hits] : hit_view_map)
    {
        float drift_min = std::numeric_limits<float>::max();
        float drift_max = -std::numeric_limits<float>::max();
        float wire_min, wire_max;
        this->getHitRegion(evt, hits, drift_min, drift_max, wire_min, wire_max);

        torch::Tensor network_input;
        std::vector<std::pair<int, int>> pixel_vector;
        this->makeNetworkInput(evt, hits, view, drift_min, drift_max, wire_min, wire_max, network_input, pixel_vector);

        torch::Tensor output;
        if (view == common::TPC_VIEW_U)
            output = _model_u->forward({network_input}).toTensor();
        else if (view == common::TPC_VIEW_V)
            output = _model_v->forward({network_input}).toTensor();
        else if (view == common::TPC_VIEW_W)
            output = _model_w->forward({network_input}).toTensor();

        torch::Tensor predicted_classes = torch::argmax(output, 1); 
        for (size_t i = 0; i < hits.size(); ++i)
        {
            const auto& hit = hits[i];
            int predicted_class = predicted_classes[0][pixel_vector[i].first][pixel_vector[i].second].item<int>();

            if (predicted_class == 1) 
                muon_hits.push_back(hit); 
            else if (predicted_class == 2) 
                pion_plus_hits.push_back(hit);  
            else if (predicted_class == 3) 
                pion_minus_hits.push_back(hit);  
        }
    }
}

void ConvNetworkAlgorithm::makeNetworkInput(const art::Event& evt, const std::vector<art::Ptr<recob::Hit>>& hit_list, const common::PandoraView view, const float x_min, const float x_max, const float z_min, const float z_max, torch::Tensor& network_input, std::vector<std::pair<int, int>>& pixel_vector)
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
    }

    for (int row = 0; row < _height; ++row)
    {
        for (int col = 0; col < _width; ++col)
        {
            const float value{accessor[0][0][row][col]};
            if (value > 0)
                pixel_vector.emplace_back(std::make_pair(row, col));
        }
    }
}

void ConvNetworkAlgorithm::beginJob() 
{}

void ConvNetworkAlgorithm::endJob() 
{}

DEFINE_ART_MODULE(ConvNetworkAlgorithm)