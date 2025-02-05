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

#include "CommonFunctions/Pandora.h"
#include "CommonFunctions/Scatters.h"
#include "CommonFunctions/Corrections.h"
#include "CommonFunctions/Region.h"
#include "CommonFunctions/Types.h"
#include "CommonFunctions/Geometry.h"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"

#include "SignatureTools/SignatureToolBase.h"

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

class ImageTrainingAnalyser : public art::EDAnalyzer
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

    using ProxyPfpColl_t = common::ProxyPfpColl_t;
    using ProxyPfpElem_t = common::ProxyPfpElem_t;

private:
    std::string _training_output_file;
    TFile* _root_file;
    TTree* _image_tree;
    std::unique_ptr<image::ImageTrainingHandler> _image_handler;

    int _image_width, _image_height;
    float _drift_step;
    float _wire_pitch_u, _wire_pitch_v, _wire_pitch_w;

    std::vector<std::unique_ptr<signature::SignatureToolBase>> _signatureToolsVec; 

    std::string _bad_channel_file; 
    std::vector<bool> _bad_channel_mask;

    const geo::GeometryCore* _geo;
    const detinfo::DetectorProperties* _detp;

    std::map<int, unsigned int> _pfpmap;

    art::InputTag _WREproducer, /*_SCHproducer,*/ _HITproducer, _MCPproducer, _MCTproducer, _BKTproducer, _PFPproducer, _CLSproducer, _SHRproducer, _SLCproducer, _VTXproducer, _PCAproducer, _TRKproducer;

    void initialiseBadChannelMask();

    void filterBadChannels(std::vector<art::Ptr<recob::Wire>>& wires);
    void filterBadSimChannels(std::vector<art::Ptr<sim::SimChannel>>& sim_channels);
    void produceTrainingSample(art::Event const& e);
    void addDaughters(const ProxyPfpElem_t &pfp_pxy, const ProxyPfpColl_t &pfp_pxy_col, std::vector<ProxyPfpElem_t> &slice_v);
    void buildPFPMap(const ProxyPfpColl_t &pfp_pxy_col);
};

ImageTrainingAnalyser::ImageTrainingAnalyser(fhicl::ParameterSet const& pset)
    : EDAnalyzer{pset}
    , _training_output_file{pset.get<std::string>("TrainingOutputFile", "training_output.root")}
    , _bad_channel_file{pset.get<std::string>("BadChannelFile", "badchannels.txt")} 
    , _image_width{pset.get<int>("ImageWidth", 512)}
    , _image_height{pset.get<int>("ImageHeight", 512)}
    , _WREproducer{pset.get<art::InputTag>("WREproducer", "butcher")}
    //, _SCHproducer{pset.get<art::InputTag>("SCHproducer", "simpleSC")}
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
{
    std::cout << "initialising...." << std::endl;
    const fhicl::ParameterSet &tool_psets = pset.get<fhicl::ParameterSet>("SignatureTools");
    for (auto const &tool_pset_label : tool_psets.get_pset_names())
    {
        auto const tool_pset = tool_psets.get<fhicl::ParameterSet>(tool_pset_label);
        _signatureToolsVec.push_back(art::make_tool<::signature::SignatureToolBase>(tool_pset));
    }

    _geo = art::ServiceHandle<geo::Geometry>()->provider();

    //_drift_step = (_detp->SamplingRate()/1000.) * _detp->DriftVelocity(_detp->Efield(), _detp->Temperature());
    _drift_step = 0.5;
    std::cout << "Drift step: " << _drift_step << std::endl;
    _wire_pitch_u = _geo->WirePitch(geo::kU);                 // U plane
    _wire_pitch_v = _geo->WirePitch(geo::kV);                 // V plane
    _wire_pitch_w = _geo->WirePitch(geo::kW);                 // W plane

    size_t n_channels = _geo->Nchannels();
    _bad_channel_mask.resize(n_channels, false);
    this->initialiseBadChannelMask();
    std::cout << "finished initialising.." << std::endl;
}

void ImageTrainingAnalyser::initialiseBadChannelMask()
{
    if (!_bad_channel_file.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fullname;
        sp.find_file(_bad_channel_file, fullname);
        if (fullname.empty()) {
            throw cet::exception("ImageTrainingAnalyser") << "Bad channel file not found: " << _bad_channel_file;
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
    _image_handler = std::make_unique<image::ImageTrainingHandler>(_image_tree, *_geo);
}

void ImageTrainingAnalyser::analyze(const art::Event& e) 
{   
    _image_handler->reset();
    this->produceTrainingSample(e);
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

void ImageTrainingAnalyser::buildPFPMap(const ProxyPfpColl_t &pfp_pxy_col)
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

void ImageTrainingAnalyser::addDaughters(const ProxyPfpElem_t &pfp_pxy,
                                           const ProxyPfpColl_t &pfp_pxy_col,
                                           std::vector<ProxyPfpElem_t> &slice_v)
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

void ImageTrainingAnalyser::produceTrainingSample(const art::Event& e) 
{
    signature::Pattern pattern;
    bool pattern_found = true;

    for (auto& signatureTool : _signatureToolsVec) {
        signature::Signature signature; 
        if (!signatureTool->constructSignature(e, signature)) {
            pattern_found = false;
            break;
        }
        pattern.emplace_back(signatureTool->getSignatureType(), signature);
    }

    if (!pattern_found && !pattern.empty())
    {
        pattern.clear();
        return;
    }

    int run = e.run();
    int subrun = e.subRun();
    int event = e.event();

    std::vector<art::Ptr<recob::Wire>> wire_vec;
    auto wireHandle = e.getValidHandle<std::vector<recob::Wire>>(_WREproducer);
    if (wireHandle) 
        art::fill_ptr_vector(wire_vec, wireHandle);
    else
        return;

    std::cout << "Wire size: " << wire_vec.size() << std::endl;

    this->filterBadChannels(wire_vec);

    /*std::vector<art::Ptr<sim::SimChannel>> sim_channel_vec;
    auto simChannelHandle = e.getValidHandle<std::vector<sim::SimChannel>>(_SCHproducer);
    if (simChannelHandle) 
        art::fill_ptr_vector(sim_channel_vec, simChannelHandle);
    else
        return;

    this->filterBadSimChannels(sim_channel_vec);*/

    common::ProxyPfpColl_t const &pfp_proxy = proxy::getCollection<std::vector<recob::PFParticle>>(e, _PFPproducer,
                                                        proxy::withAssociated<larpandoraobj::PFParticleMetadata>(_PFPproducer),
                                                        proxy::withAssociated<recob::Cluster>(_CLSproducer),
                                                        proxy::withAssociated<recob::Slice>(_SLCproducer),
                                                        proxy::withAssociated<recob::Track>(_TRKproducer),
                                                        proxy::withAssociated<recob::Vertex>(_VTXproducer),
                                                        proxy::withAssociated<recob::PCAxis>(_PCAproducer),
                                                        proxy::withAssociated<recob::Shower>(_SHRproducer),
                                                        proxy::withAssociated<recob::SpacePoint>(_PFPproducer));

    this->buildPFPMap(pfp_proxy);
    std::vector<ProxyPfpElem_t> neutrino_slice;
    for (const ProxyPfpElem_t &pfp_pxy : pfp_proxy)
    {
        const auto &pfParticleMetadataList = pfp_pxy.get<larpandoraobj::PFParticleMetadata>();

        if (pfp_pxy->IsPrimary() == false)
            continue;

        auto primary_pdg = fabs(pfp_pxy->PdgCode());
        if ((primary_pdg == 12) || (primary_pdg == 14))
            this->addDaughters(pfp_pxy, pfp_proxy, neutrino_slice);
    } 

    std::vector<art::Ptr<recob::Hit>> neutrino_hits;
    common::ProxyClusColl_t const& clus_proxy = proxy::getCollection<std::vector<recob::Cluster>>(
        e, _CLSproducer, proxy::withAssociated<recob::Hit>(_CLSproducer)
    );

    for (const auto& pfp : neutrino_slice) {
        if (pfp->IsPrimary()) continue; 

        auto clus_pxy_v = pfp.get<recob::Cluster>();

        for (auto ass_clus : clus_pxy_v) {
            const auto& clus = clus_proxy[ass_clus.key()];
            auto clus_hit_v = clus.get<recob::Hit>();

            neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
        }
    }

    std::cout << "Neutrino hits size: " << neutrino_hits.size() << std::endl;

    double sum_charge_u = 0.0, sum_wire_u = 0.0, sum_drift_u = 0.0;
    double sum_charge_v = 0.0, sum_wire_v = 0.0, sum_drift_v = 0.0;
    double sum_charge_w = 0.0, sum_wire_w = 0.0, sum_drift_w = 0.0;

    for (const auto& hit : neutrino_hits) {
        double charge = hit->Integral();
        common::PandoraView pandora_view = common::GetPandoraView(hit);
        TVector3 hit_pos = common::GetPandoraHitPosition(e, hit, pandora_view);

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

    std::cout << "Centroids - U: (" << sum_wire_u / sum_charge_u << ", " << sum_drift_u / sum_charge_u << ")"
              << ", V: (" << sum_wire_v / sum_charge_v << ", " << sum_drift_v / sum_charge_v << ")"
              << ", W: (" << sum_wire_w / sum_charge_w << ", " << sum_drift_w / sum_charge_w << ")\n";

    std::vector<image::ImageProperties> properties;
    properties.emplace_back(centroid_wire_u, centroid_drift_u,
        _image_height, _image_width, _wire_pitch_u, _drift_step, geo::kU
    );
    properties.emplace_back(centroid_wire_v, centroid_drift_v,
        _image_height, _image_width, _wire_pitch_v, _drift_step, geo::kV
    );
    properties.emplace_back(centroid_wire_w, centroid_drift_w,
        _image_height, _image_width, _wire_pitch_w, _drift_step, geo::kW
    );

    _image_handler->add(run, subrun, event, pattern_found, wire_vec, /*sim_channel_vec,*/ properties);
}

void ImageTrainingAnalyser::endJob() 
{
    if (_root_file) {
        _root_file->cd();         
        _image_tree->Write();     
        _root_file->Close();      

        delete _root_file;        
        _root_file = nullptr;     
    }
}

DEFINE_ART_MODULE(ImageTrainingAnalyser)