#ifndef ANALYSIS_IMAGE_CXX
#define ANALYSIS_IMAGE_CXX

#include <vector>
#include <string>
#include <map>
#include <array>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <limits>
#include <optional>
#include <numeric>
#include <set>
#include <fstream>
#include <experimental/filesystem>
#include <sstream>
#include <H5Cpp.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <lardataobj/AnalysisBase/BackTrackerMatchingData.h>
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/RecoBaseProxy/ProxyBase.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "lardata/Utilities/FindManyInChainP.h"
#include "AnalysisToolBase.h"
#include "../CommonDefs/Pandora.h"
#include "../CommonDefs/Types.h"
#include "../Image.h"
#include "../ImageAlgorithm.h"

namespace fs = std::experimental::filesystem;

namespace analysis {

    struct ModelConfig {
        std::string name;
        std::string weights_file;
    };

    class ImageAnalysis : public AnalysisToolBase {
    public:
        explicit ImageAnalysis(fhicl::ParameterSet const& p);
        virtual ~ImageAnalysis() = default;

        void configure(const fhicl::ParameterSet& p) override;
        void analyseEvent(art::Event const& e, bool _is_data) override {}
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, bool _is_data, bool selected) override;
        void setBranches(TTree* _tree) override;
        void resetTTree(TTree* _tree) override;

    private:
        art::InputTag fPFPproducer;
        art::InputTag fCLSproducer;
        art::InputTag fSLCprodcuer;
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        art::InputTag fBKTproducer;
        std::string fBadChannelFile;
        std::set<unsigned int> fBadChannels;
        std::vector<ModelConfig> fModels;
        int _image_width;
        int _image_height;
        float _adc_image_threshold;
        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;
        float _drift_step;
        float _wire_pitch_u;
        float _wire_pitch_v;
        float _wire_pitch_w;
        float _reco_neutrino_vertex_x;
        float _reco_neutrino_vertex_y;
        float _reco_neutrino_vertex_z;
        std::vector<float> _detector_image_u;
        std::vector<float> _detector_image_v;
        std::vector<float> _detector_image_w;
        std::vector<int> _semantic_image_u;
        std::vector<int> _semantic_image_v;
        std::vector<int> _semantic_image_w;
        std::vector<float> _event_detector_image_u;
        std::vector<float> _event_detector_image_v;
        std::vector<float> _event_detector_image_w;
        std::vector<int> _event_semantic_image_u;
        std::vector<int> _event_semantic_image_v;
        std::vector<int> _event_semantic_image_w;
        float _event_adc_u;
        float _event_adc_v;
        float _event_adc_w;
        std::vector<int> _slice_semantic_counts_u;
        std::vector<int> _slice_semantic_counts_v;
        std::vector<int> _slice_semantic_counts_w;
        std::vector<int> _event_semantic_counts_u;
        std::vector<int> _event_semantic_counts_v;
        std::vector<int> _event_semantic_counts_w;
        bool _is_vtx_in_image_u;
        bool _is_vtx_in_image_v;
        bool _is_vtx_in_image_w;
        std::unique_ptr<TruthLabelClassifier> _semantic_classifier;
        std::unordered_map<std::string, float> _inference_scores;

        void loadBadChannels(const std::string& filename);
        std::vector<art::Ptr<recob::Hit>> collectAllHits(const art::Event& e);
        std::vector<art::Ptr<recob::Hit>> collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
        std::pair<double, double> calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        void constructPixelImages(const art::Event& e, const std::vector<art::Ptr<recob::Hit>>& hits, const std::vector<ImageProperties>& properties, std::vector<Image<float>>& detector_images, std::vector<Image<int>>& semantic_images, bool is_data);
        float runInference(const std::vector<Image<float>>& detector_images, const std::string& absolute_scratch_dir, const std::string& work_dir, const std::string& weights_file);
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer");
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fSLCprodcuer = p.get<art::InputTag>("SLCproducer");
        fHITproducer = p.get<art::InputTag>("HITproducer");
        fWIREproducer = p.get<art::InputTag>("WIREproducer");
        fMCPproducer = p.get<art::InputTag>("MCPproducer");
        fBKTproducer = p.get<art::InputTag>("BKTproducer");
        fBadChannelFile = p.get<std::string>("BadChannelFile");
        auto model_psets = p.get<std::vector<fhicl::ParameterSet>>("Models", {});
        for (const auto& ps : model_psets) {
            fModels.push_back({ps.get<std::string>("name"), ps.get<std::string>("weights_file")});
        }
        _image_width = 512;
        _image_height = 512;
        _adc_image_threshold = 1.0;
        _geo = art::ServiceHandle<geo::Geometry>()->provider();
        _detp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->provider();
        auto clock = art::ServiceHandle<detinfo::DetectorClocksService>()->provider();
        double tick_period = clock->TPCClock().TickPeriod();
        double drift_velocity = _detp->DriftVelocity();
        _drift_step = tick_period * drift_velocity * 1e1;
        _wire_pitch_u = _geo->WirePitch(geo::kU);
        _wire_pitch_v = _geo->WirePitch(geo::kV);
        _wire_pitch_w = _geo->WirePitch(geo::kW);
        _semantic_classifier = std::make_unique<TruthLabelClassifier>(fMCPproducer);
        if (!fBadChannelFile.empty()) {
            this->loadBadChannels(fBadChannelFile);
        }
    }

    void ImageAnalysis::loadBadChannels(const std::string& filename) {
        fBadChannels.clear();
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            throw art::Exception(art::errors::Configuration) << "Could not open bad channels file: " << filename;
        }
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty() || line.find('#') != std::string::npos) {
                continue;
            }
            std::stringstream ss(line);
            unsigned int first_ch, second_ch;
            ss >> first_ch;
            if (ss >> second_ch) {
                for (unsigned int ch = first_ch; ch <= second_ch; ++ch) {
                    fBadChannels.insert(ch);
                }
            } else {
                fBadChannels.insert(first_ch);
            }
        }
    }

    void ImageAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("reco_neutrino_vertex_x", &_reco_neutrino_vertex_x, "reco_neutrino_vertex_x/F");
        _tree->Branch("reco_neutrino_vertex_y", &_reco_neutrino_vertex_y, "reco_neutrino_vertex_y/F");
        _tree->Branch("reco_neutrino_vertex_z", &_reco_neutrino_vertex_z, "reco_neutrino_vertex_z/F");
        _tree->Branch("detector_image_u", &_detector_image_u);
        _tree->Branch("detector_image_v", &_detector_image_v);
        _tree->Branch("detector_image_w", &_detector_image_w);
        _tree->Branch("semantic_image_u", &_semantic_image_u);
        _tree->Branch("semantic_image_v", &_semantic_image_v);
        _tree->Branch("semantic_image_w", &_semantic_image_w);
        _tree->Branch("event_detector_image_u", &_event_detector_image_u);
        _tree->Branch("event_detector_image_v", &_event_detector_image_v);
        _tree->Branch("event_detector_image_w", &_event_detector_image_w);
        _tree->Branch("event_semantic_image_u", &_event_semantic_image_u);
        _tree->Branch("event_semantic_image_v", &_event_semantic_image_v);
        _tree->Branch("event_semantic_image_w", &_event_semantic_image_w);
        _tree->Branch("event_adc_u", &_event_adc_u, "event_adc_u/F");
        _tree->Branch("event_adc_v", &_event_adc_v, "event_adc_v/F");
        _tree->Branch("event_adc_w", &_event_adc_w, "event_adc_w/F");
        _tree->Branch("slice_semantic_counts_u", &_slice_semantic_counts_u);
        _tree->Branch("slice_semantic_counts_v", &_slice_semantic_counts_v);
        _tree->Branch("slice_semantic_counts_w", &_slice_semantic_counts_w);
        _tree->Branch("event_semantic_counts_u", &_event_semantic_counts_u);
        _tree->Branch("event_semantic_counts_v", &_event_semantic_counts_v);
        _tree->Branch("event_semantic_counts_w", &_event_semantic_counts_w);
        _tree->Branch("is_vtx_in_image_u", &_is_vtx_in_image_u, "is_vtx_in_image_u/O");
        _tree->Branch("is_vtx_in_image_v", &_is_vtx_in_image_v, "is_vtx_in_image_v/O");
        _tree->Branch("is_vtx_in_image_w", &_is_vtx_in_image_w, "is_vtx_in_image_w/O");
        for (const auto& model : fModels) {
            _tree->Branch(("inference_score_" + model.name).c_str(), &_inference_scores[model.name], "score/F");
        }
    }

    void ImageAnalysis::resetTTree(TTree* _tree) {
        _reco_neutrino_vertex_x = std::numeric_limits<float>::lowest();
        _reco_neutrino_vertex_y = std::numeric_limits<float>::lowest();
        _reco_neutrino_vertex_z = std::numeric_limits<float>::lowest();
        _detector_image_u.clear();
        _detector_image_v.clear();
        _detector_image_w.clear();
        _semantic_image_u.clear();
        _semantic_image_v.clear();
        _semantic_image_w.clear();
        _event_detector_image_u.clear();
        _event_detector_image_v.clear();
        _event_detector_image_w.clear();
        _event_semantic_image_u.clear();
        _event_semantic_image_v.clear();
        _event_semantic_image_w.clear();
        _event_adc_u = std::numeric_limits<float>::lowest();
        _event_adc_v = std::numeric_limits<float>::lowest();
        _event_adc_w = std::numeric_limits<float>::lowest();
        _slice_semantic_counts_u.clear();
        _slice_semantic_counts_v.clear();
        _slice_semantic_counts_w.clear();
        _event_semantic_counts_u.clear();
        _event_semantic_counts_v.clear();
        _event_semantic_counts_w.clear();
        _is_vtx_in_image_u = false;
        _is_vtx_in_image_v = false;
        _is_vtx_in_image_w = false;
        _inference_scores.clear();
        for (const auto& model : fModels) {
            _inference_scores[model.name] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, bool _is_data, bool selected) {
        for (const auto& pfp : pfp_pxy_v) {
            if (pfp->IsPrimary()) {
                const auto& vtx = pfp.get<recob::Vertex>();
                if (vtx.size() > 0) {
                    const auto& vtx_pos = vtx[0]->position();
                    _reco_neutrino_vertex_x = vtx_pos.X();
                    _reco_neutrino_vertex_y = vtx_pos.Y();
                    _reco_neutrino_vertex_z = vtx_pos.Z();
                }
                break;
            }
        }

        std::vector<art::Ptr<recob::Hit>> all_hits = this->collectAllHits(e);
        std::vector<art::Ptr<recob::Hit>> neutrino_hits = this->collectSliceHits(e, pfp_pxy_v);

        auto [centroid_wire_u, centroid_drift_u] = this->calculateChargeCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateChargeCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateChargeCentroid(e, common::TPC_VIEW_W, neutrino_hits);

        std::vector<ImageProperties> properties;
        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_width, _image_height, _drift_step, _wire_pitch_u, geo::kU);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_width, _image_height, _drift_step, _wire_pitch_v, geo::kV);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_width, _image_height, _drift_step, _wire_pitch_w, geo::kW);

        std::vector<Image<float>> detector_images;
        std::vector<Image<int>> semantic_images;
        this->constructPixelImages(e, neutrino_hits, properties, detector_images, semantic_images, _is_data);

        std::vector<Image<float>> event_detector_images;
        std::vector<Image<int>> event_semantic_images;
        this->constructPixelImages(e, all_hits, properties, event_detector_images, event_semantic_images, _is_data);

        _detector_image_u = detector_images[0].data();
        _detector_image_v = detector_images[1].data();
        _detector_image_w = detector_images[2].data();
        _semantic_image_u = semantic_images[0].data();
        _semantic_image_v = semantic_images[1].data();
        _semantic_image_w = semantic_images[2].data();
        _event_detector_image_u = event_detector_images[0].data();
        _event_detector_image_v = event_detector_images[1].data();
        _event_detector_image_w = event_detector_images[2].data();
        _event_semantic_image_u = event_semantic_images[0].data();
        _event_semantic_image_v = event_semantic_images[1].data();
        _event_semantic_image_w = event_semantic_images[2].data();
        _event_adc_u = std::accumulate(_event_detector_image_u.begin(), _event_detector_image_u.end(), 0.0f);
        _event_adc_v = std::accumulate(_event_detector_image_v.begin(), _event_detector_image_v.end(), 0.0f);
        _event_adc_w = std::accumulate(_event_detector_image_w.begin(), _event_detector_image_w.end(), 0.0f);

        if (!_is_data) {
            for (size_t i = 0; i < TruthLabelClassifier::truth_primary_label_names.size(); ++i) {
                _slice_semantic_counts_u.push_back(std::count(_semantic_image_u.begin(), _semantic_image_u.end(), i));
                _slice_semantic_counts_v.push_back(std::count(_semantic_image_v.begin(), _semantic_image_v.end(), i));
                _slice_semantic_counts_w.push_back(std::count(_semantic_image_w.begin(), _semantic_image_w.end(), i));
                _event_semantic_counts_u.push_back(std::count(_event_semantic_image_u.begin(), _event_semantic_image_u.end(), i));
                _event_semantic_counts_v.push_back(std::count(_event_semantic_image_v.begin(), _event_semantic_image_v.end(), i));
                _event_semantic_counts_w.push_back(std::count(_event_semantic_image_w.begin(), _event_semantic_image_w.end(), i));
            }
        }

        std::string scratch_dir = ".";
        const char* scratch_dir_cstr = std::getenv("_CONDOR_SCRATCH_DIR");
        if (scratch_dir_cstr) {
            scratch_dir = std::string(scratch_dir_cstr);
        }
        char scratch_path_buffer[PATH_MAX];
        realpath(scratch_dir.c_str(), scratch_path_buffer);
        std::string absolute_scratch_dir(scratch_path_buffer);

        char cwd_buffer[PATH_MAX];
        getcwd(cwd_buffer, sizeof(cwd_buffer));
        std::string work_dir(cwd_buffer);

        for (const auto& model : fModels) {
            float score = this->runInference(detector_images, absolute_scratch_dir, work_dir, model.weights_file);
            _inference_scores[model.name] = score;
        }

        if (_reco_neutrino_vertex_x != std::numeric_limits<float>::lowest()) {
            TVector3 vtx_pos_3d(_reco_neutrino_vertex_x, _reco_neutrino_vertex_y, _reco_neutrino_vertex_z);
            TVector3 vtx_proj_u = common::ProjectToWireView(vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_U);
            _is_vtx_in_image_u = (properties[0].row(vtx_proj_u.X()) != static_cast<size_t>(-1) && properties[0].col(vtx_proj_u.Z()) != static_cast<size_t>(-1));
            TVector3 vtx_proj_v = common::ProjectToWireView(vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_V);
            _is_vtx_in_image_v = (properties[1].row(vtx_proj_v.X()) != static_cast<size_t>(-1) && properties[1].col(vtx_proj_v.Z()) != static_cast<size_t>(-1));
            TVector3 vtx_proj_w = common::ProjectToWireView(vtx_pos_3d.X(), vtx_pos_3d.Y(), vtx_pos_3d.Z(), common::TPC_VIEW_W);
            _is_vtx_in_image_w = (properties[2].row(vtx_proj_w.X()) != static_cast<size_t>(-1) && properties[2].col(vtx_proj_w.Z()) != static_cast<size_t>(-1));
        }
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectAllHits(const art::Event& e) {
        std::vector<art::Ptr<recob::Hit>> all_hits;
        auto hit_handle = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
        for (size_t i = 0; i < hit_handle->size(); ++i) {
            all_hits.emplace_back(hit_handle, i);
        }
        return all_hits;
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        auto pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);
        art::FindManyP<recob::Slice> sliceAssoc(pfpHandle, e, fPFPproducer);
        size_t pfpIndex = pfp_pxy_v[0].index();
        const auto& slices = sliceAssoc.at(pfpIndex);
        if (slices.empty()) return {};
        const art::Ptr<recob::Slice>& slice = slices[0];
        auto sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSLCprodcuer);
        art::FindManyP<recob::Hit> hitAssoc(sliceHandle, e, fSLCprodcuer);
        const std::vector<art::Ptr<recob::Hit>>& sliceHits = hitAssoc.at(slice.key());
        neutrino_hits.reserve(sliceHits.size());
        for (const auto& hit : sliceHits) {
            neutrino_hits.push_back(hit);
        }
        return neutrino_hits;
    }

    std::pair<double, double> ImageAnalysis::calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits) {
        double sum_charge = 0.0;
        double sum_wire = 0.0;
        double sum_drift = 0.0;
        for (const auto& hit : hits) {
            if (fBadChannels.count(hit->Channel())) {
                continue;
            }
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

    void ImageAnalysis::constructPixelImages(const art::Event& e,
                                             const std::vector<art::Ptr<recob::Hit>>& hits,
                                             const std::vector<ImageProperties>& properties,
                                             std::vector<Image<float>>& detector_images,
                                             std::vector<Image<int>>& semantic_images,
                                             bool is_data) {
        detector_images.clear();
        semantic_images.clear();
        for (const auto& prop : properties) {
            Image<float> detector_image(prop);
            detector_image.clear(0.0);
            detector_images.push_back(std::move(detector_image));
            Image<int> semantic_image(prop);
            semantic_image.clear(static_cast<int>(TruthLabelClassifier::TruthPrimaryLabel::Empty));
            semantic_images.push_back(std::move(semantic_image));
        }

        auto wire_vector = e.getValidHandle<std::vector<recob::Wire>>(fWIREproducer);
        auto hit_vector = e.getValidHandle<std::vector<recob::Hit>>(fHITproducer);
        art::Handle<std::vector<simb::MCParticle>> mcp_vector;
        bool has_mcps = e.getByLabel(fMCPproducer, mcp_vector);
        art::FindManyP<recob::Hit> wire_hit_assoc(wire_vector, e, fHITproducer);
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hit_vector, e, fBKTproducer);

        std::vector<TruthLabelClassifier::TruthPrimaryLabel> semantic_label_vector;
        if (!is_data && has_mcps && mcp_vector.isValid() && _semantic_classifier) {
            semantic_label_vector = _semantic_classifier->classifyParticles(e);
        }

        std::map<int, size_t> trackid_to_index;
        if (!is_data && has_mcps && mcp_vector.isValid()) {
            for (size_t i = 0; i < mcp_vector->size(); ++i) {
                trackid_to_index[mcp_vector->at(i).TrackId()] = i;
            }
        }

        std::set<art::Ptr<recob::Hit>> hit_set(hits.begin(), hits.end());
        for (size_t wire_idx = 0; wire_idx < wire_vector->size(); ++wire_idx) {
            const auto& wire = wire_vector->at(wire_idx);
            auto ch_id = wire.Channel();
            if (fBadChannels.count(ch_id)) {
                continue;
            }
            std::vector<geo::WireID> wire_ids = _geo->ChannelToWire(ch_id);
            if (wire_ids.empty()) continue;
            geo::View_t view = _geo->View(wire_ids.front().planeID());
            size_t view_idx = static_cast<size_t>(view);
            const geo::WireGeo* wire_geo = _geo->WirePtr(wire_ids.front());
            TVector3 center = wire_geo->GetCenter();
            TVector3 wire_center(center.X(), center.Y(), center.Z());
            double wire_coord = (view == geo::kW) ? wire_center.Z() :
                                (view == geo::kU) ? (wire_center.Z() * std::cos(1.04719758034) - wire_center.Y() * std::sin(1.04719758034)) :
                                (wire_center.Z() * std::cos(-1.04719758034) - wire_center.Y() * std::sin(-1.04719758034));
            auto hits_for_wire = wire_hit_assoc.at(wire_idx);
            std::vector<art::Ptr<recob::Hit>> filtered_hits;
            for (const auto& hit : hits_for_wire) {
                if (hit_set.count(hit)) {
                    filtered_hits.push_back(hit);
                }
            }
            std::sort(filtered_hits.begin(), filtered_hits.end(), [](const art::Ptr<recob::Hit>& a, const art::Ptr<recob::Hit>& b) {
                return a->StartTick() < b->StartTick();
            });
            size_t hit_index = 0;
            for (const auto& range : wire.SignalROI().get_ranges()) {
                const auto& adcs = range.data();
                int start_tick = range.begin_index();
                for (size_t adc_index = 0; adc_index < adcs.size(); ++adc_index) {
                    int tick = start_tick + adc_index;
                    double x_drift_pos = _detp->ConvertTicksToX(static_cast<double>(tick), wire_ids.front().planeID());
                    size_t row = properties[view_idx].row(x_drift_pos);
                    size_t col = properties[view_idx].col(wire_coord);
                    if (row == static_cast<size_t>(-1) || col == static_cast<size_t>(-1)) continue;
                    while (hit_index < filtered_hits.size() && filtered_hits[hit_index]->EndTick() <= tick) {
                        ++hit_index;
                    }
                    if (hit_index < filtered_hits.size() &&
                        filtered_hits[hit_index]->StartTick() <= tick &&
                        tick < filtered_hits[hit_index]->EndTick()) {
                        if (adcs[adc_index] > _adc_image_threshold) {
                            detector_images[view_idx].set(row, col, adcs[adc_index]);
                            if (!is_data) {
                                const art::Ptr<recob::Hit>& matched_hit = filtered_hits[hit_index];
                                TruthLabelClassifier::TruthPrimaryLabel semantic_pixel_label = TruthLabelClassifier::TruthPrimaryLabel::Cosmic;
                                if (has_mcps && mcp_vector.isValid() && _semantic_classifier) {
                                    std::vector<art::Ptr<simb::MCParticle>> mcp_particles_ass_to_hit;
                                    std::vector<anab::BackTrackerHitMatchingData const*> bkth_data_ass_to_hit;
                                    mcp_bkth_assoc.get(matched_hit.key(), mcp_particles_ass_to_hit, bkth_data_ass_to_hit);
                                    if (!bkth_data_ass_to_hit.empty()) {
                                        float max_ide_fraction = -1.0;
                                        int best_match_track_id = -1;
                                        for (size_t i_bkth = 0; i_bkth < bkth_data_ass_to_hit.size(); ++i_bkth) {
                                            if (bkth_data_ass_to_hit[i_bkth] && bkth_data_ass_to_hit[i_bkth]->ideFraction > max_ide_fraction) {
                                                max_ide_fraction = bkth_data_ass_to_hit[i_bkth]->ideFraction;
                                                best_match_track_id = mcp_particles_ass_to_hit[i_bkth]->TrackId();
                                            }
                                        }
                                        if (best_match_track_id != -1) {
                                            auto it_trackid = trackid_to_index.find(best_match_track_id);
                                            if (it_trackid != trackid_to_index.end()) {
                                                size_t particle_mcp_idx = it_trackid->second;
                                                if (particle_mcp_idx < semantic_label_vector.size()) {
                                                    semantic_pixel_label = semantic_label_vector[particle_mcp_idx];
                                                }
                                            }
                                        }
                                    }
                                }
                                semantic_images[view_idx].set(row, col, static_cast<int>(semantic_pixel_label), false);
                            }
                        }
                    }
                }
            }
        }
    }

    float ImageAnalysis::runInference(const std::vector<Image<float>>& detector_images,
        const std::string& absolute_scratch_dir,
        const std::string& work_dir,
        const std::string& weights_file) {
        std::string temp_in = absolute_scratch_dir + "/temp_test_in.h5";
        std::string temp_out = absolute_scratch_dir + "/temp_test_out.txt";
        std::string script_stdout = absolute_scratch_dir + "/py_script.out";
        std::string script_stderr = absolute_scratch_dir + "/py_script.err";

        std::string tree_name = "imagetree";
        std::vector<float> image_u = detector_images[0].data();
        std::vector<float> image_v = detector_images[1].data();
        std::vector<float> image_w = detector_images[2].data();
        const hsize_t dims[2] = {1, static_cast<hsize_t>(image_u.size())};
        H5::H5File h5file(temp_in.c_str(), H5F_ACC_TRUNC);
        H5::Group group = h5file.createGroup("/" + tree_name);
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dset_u = group.createDataSet("image_u", H5::PredType::NATIVE_FLOAT, dataspace);
        H5::DataSet dset_v = group.createDataSet("image_v", H5::PredType::NATIVE_FLOAT, dataspace);
        H5::DataSet dset_w = group.createDataSet("image_w", H5::PredType::NATIVE_FLOAT, dataspace);
        dset_u.write(image_u.data(), H5::PredType::NATIVE_FLOAT);
        dset_v.write(image_v.data(), H5::PredType::NATIVE_FLOAT);
        dset_w.write(image_w.data(), H5::PredType::NATIVE_FLOAT);
        h5file.close();

        std::string branch_name = "image_w";
        std::string container = "/cvmfs/uboone.opensciencegrid.org/containers/lantern_v2_me_06_03_prod";
        std::string wrapper_script = "run_strangeness_inference.sh";

        std::string bind_paths = absolute_scratch_dir;
        const char* bind_env = std::getenv("APPTAINER_BINDPATH");
        bool need_cvmfs = true;
        
        if (bind_env && std::string(bind_env).find("/cvmfs") != std::string::npos) {
            need_cvmfs = false;
        } else {
            std::ifstream mounts("/proc/mounts");
            std::string line;
        
            while (std::getline(mounts, line)) {
                if (line.find(" /cvmfs ") != std::string::npos) {
                    need_cvmfs = false;
                    break;
                }
            }
        }
        
        if (need_cvmfs) {
            bind_paths = "/cvmfs," + bind_paths;
        }

        std::string command = "apptainer exec --cleanenv --bind " + bind_paths + " " +
            container + " " +
            "/bin/bash ./" + wrapper_script + " " +
            temp_in + " " +
            temp_out + " " +
            work_dir + "/" + weights_file + " " +
            tree_name + " " +
            branch_name +
            " > " + script_stdout + " 2> " + script_stderr;

        mf::LogInfo("ImageAnalysis") << "Executing inference: " << command;

        auto start = std::chrono::steady_clock::now();
        int code = std::system(command.c_str());
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration<double>(end - start).count();

        bool success = (code == 0) && fs::exists(temp_out);
        if (!success) {
            std::ifstream error_stream(script_stderr);
            std::string error_message((std::istreambuf_iterator<char>(error_stream)), std::istreambuf_iterator<char>());
            throw art::Exception(art::errors::LogicError)
                << "Inference script failed with exit code " << code
                << "\n--- Script stderr ---\n" << error_message << "\n--- End Script stderr ---";
        }   

        std::ifstream result_stream(temp_out.c_str());
        if (!result_stream) {
            throw art::Exception(art::errors::LogicError) << "Could not open temporary result file: " << temp_out;
        }
        
        float score;
        result_stream >> score;

        mf::LogInfo("ImageAnalysis") << "Inference time: " << duration << " seconds";
        mf::LogInfo("ImageAnalysis") << "Predicted score: " << score;

        remove(temp_in.c_str());
        remove(temp_out.c_str());
        remove(script_stdout.c_str());
        remove(script_stderr.c_str());

        return score;
    }


    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif
