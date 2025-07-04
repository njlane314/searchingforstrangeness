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

#ifdef ClassDef
#undef ClassDef
#endif
#include "torch/torch.h"
#include "torch/script.h"

namespace analysis 
{
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
        art::InputTag fSLCEproducer;
        art::InputTag fHITproducer;
        art::InputTag fWIREproducer;
        art::InputTag fMCPproducer;
        std::string fBackTrackerLabel;
        bool fRunInference;

        int _image_width;
        int _image_height;
        float _adc_image_threshold;

        const geo::GeometryCore* _geo;
        const detinfo::DetectorProperties* _detp;

        float _drift_step;
        float _wire_pitch_u;
        float _wire_pitch_v;
        float _wire_pitch_w;

        TTree* _tree = nullptr;

        std::vector<float> _detector_image_u;
        std::vector<float> _detector_image_v;
        std::vector<float> _detector_image_w;
        std::vector<int> _semantic_image_u;
        std::vector<int> _semantic_image_v;
        std::vector<int> _semantic_image_w;

        std::unique_ptr<TruthLabelClassifier> _semantic_classifier;

        std::shared_ptr<torch::jit::script::Module> _classifier_model;
        std::string _classifier_model_path;

        std::vector<float> _resnet_probabilities;

        std::vector<art::Ptr<recob::Hit>> collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v);
        std::pair<double, double> calculateChargeCentroid(const art::Event& e, common::PandoraView view, const std::vector<art::Ptr<recob::Hit>>& hits);
        void constructPixelImages(const art::Event& e, const std::vector<art::Ptr<recob::Hit>>& neutrino_hits, const std::vector<ImageProperties>& properties, std::vector<Image<float>>& detector_images, std::vector<Image<int>>& semantic_images);
        std::vector<float> runInference(const std::vector<float>& image_data);
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fPFPproducer = p.get<art::InputTag>("PFPproducer", "pandora");
        fCLSproducer = p.get<art::InputTag>("CLSproducer", "pandora");
        fSLCEproducer = p.get<art::InputTag>("SLICEproducer", "pandora");
        fHITproducer = p.get<art::InputTag>("HITproducer");
        fWIREproducer = p.get<art::InputTag>("WIREproducer");
        fMCPproducer = p.get<art::InputTag>("MCPproducer");
        fBackTrackerLabel = p.get<std::string>("BackTrackerLabel");
        fRunInference = p.get<bool>("RunInference", false);

        _image_width = 512;
        _image_height = 512;
        _adc_image_threshold = 4.0;

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

        _classifier_model_path = p.get<std::string>("ClassifierModelPath");
        if (!fRunInference) 
            return;

        try {
            _classifier_model = torch::jit::load(_classifier_model_path);
            _classifier_model->to(torch::kCPU);
            std::cout << "Successfully loaded classifier model from: " << _classifier_model_path << std::endl;
        } catch (const c10::Error& e) {
            throw std::runtime_error("Error loading the TorchScript model: " + std::string(e.what()));
        }
    }

    void ImageAnalysis::setBranches(TTree* _tree_ptr) {
        _tree = _tree_ptr;
        _tree->Branch("detector_image_u", &_detector_image_u);
        _tree->Branch("detector_image_v", &_detector_image_v);
        _tree->Branch("detector_image_w", &_detector_image_w);
        _tree->Branch("semantic_image_u", &_semantic_image_u);
        _tree->Branch("semantic_image_v", &_semantic_image_v);
        _tree->Branch("semantic_image_w", &_semantic_image_w);

        _tree->Branch("resnet_probabilities", &_resnet_probabilities);
    }

    void ImageAnalysis::resetTTree(TTree* _tree_ptr) {
        _detector_image_u.clear();
        _detector_image_v.clear();
        _detector_image_w.clear();
        _semantic_image_u.clear();
        _semantic_image_v.clear();
        _semantic_image_w.clear();

        _resnet_probabilities.clear();
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& pfp_pxy_v, bool _is_data, bool selected) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits = this->collectSliceHits(e, pfp_pxy_v);
        auto [centroid_wire_u, centroid_drift_u] = this->calculateChargeCentroid(e, common::TPC_VIEW_U, neutrino_hits);
        auto [centroid_wire_v, centroid_drift_v] = this->calculateChargeCentroid(e, common::TPC_VIEW_V, neutrino_hits);
        auto [centroid_wire_w, centroid_drift_w] = this->calculateChargeCentroid(e, common::TPC_VIEW_W, neutrino_hits);
        std::vector<ImageProperties> properties;
        properties.emplace_back(centroid_wire_u, centroid_drift_u, _image_height, _image_width, _drift_step, _wire_pitch_u, geo::kU);
        properties.emplace_back(centroid_wire_v, centroid_drift_v, _image_height, _image_width, _drift_step, _wire_pitch_v, geo::kV);
        properties.emplace_back(centroid_wire_w, centroid_drift_w, _image_height, _image_width, _drift_step, _wire_pitch_w, geo::kW);

        std::vector<Image<float>> detector_images;
        std::vector<Image<int>> semantic_images;
        this->constructPixelImages(e, neutrino_hits, properties, detector_images, semantic_images);

        _detector_image_u.clear();
        _detector_image_v.clear();
        _detector_image_w.clear();
        _semantic_image_u.clear();
        _semantic_image_v.clear();
        _semantic_image_w.clear();

        _detector_image_u = detector_images[0].data();  
        _detector_image_v = detector_images[1].data();  
        _detector_image_w = detector_images[2].data();  
        _semantic_image_u = semantic_images[0].data();
        _semantic_image_v = semantic_images[1].data();
        _semantic_image_w = semantic_images[2].data();

        if (fRunInference)
            _resnet_probabilities = this->runInference(_detector_image_w);

        if (_tree) _tree->Fill();
    }

    std::vector<float> ImageAnalysis::runInference(const std::vector<float>& image_data) {
        if (!_classifier_model || image_data.empty()) {
            return {}; 
        }

        torch::NoGradGuard no_grad;

        torch::Tensor image_tensor = torch::from_blob(
            const_cast<float*>(image_data.data()),
            {_image_height, _image_width}, 
            torch::kFloat32
        ).clone();

        image_tensor.masked_fill_(image_tensor < _adc_image_threshold, 0.0);
        image_tensor = torch::log1p(image_tensor);
        image_tensor = image_tensor.reshape({1, 1, _image_height, _image_width});

        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(image_tensor.to(torch::kCPU));
        
        at::Tensor output_logits = _classifier_model->forward(inputs).toTensor();

        at::Tensor probabilities_tensor = torch::softmax(output_logits, /*dim=*/1);
        const float* prob_data_ptr = static_cast<const float*>(probabilities_tensor.data_ptr());
        std::vector<float> probabilities(prob_data_ptr, prob_data_ptr + probabilities_tensor.numel());

        std::cout << "Probabilities: ";
        for (const auto& prob : probabilities) {
            std::cout << prob << " ";
        }
        std::cout << std::endl;

        return probabilities;
    }

    std::vector<art::Ptr<recob::Hit>> ImageAnalysis::collectSliceHits(const art::Event& e, const std::vector<common::ProxyPfpElem_t>& pfp_pxy_v) {
        std::vector<art::Ptr<recob::Hit>> neutrino_hits;
        
        auto pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPproducer);

        art::FindManyP<recob::Slice> sliceAssoc(pfpHandle, e, fPFPproducer);
        size_t pfpIndex = pfp_pxy_v[0].index();

        const auto& slices = sliceAssoc.at(pfpIndex);
        if (slices.empty()) return {};

        const art::Ptr<recob::Slice>& slice = slices[0];
        auto sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSLCEproducer);
        art::FindManyP<recob::Hit> hitAssoc(sliceHandle, e, fSLCEproducer);
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
                                         const std::vector<art::Ptr<recob::Hit>>& neutrino_hits,
                                         const std::vector<ImageProperties>& properties,
                                         std::vector<Image<float>>& detector_images,
                                         std::vector<Image<int>>& semantic_images) {
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
        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcp_bkth_assoc(hit_vector, e, fBackTrackerLabel);

        std::vector<TruthLabelClassifier::TruthPrimaryLabel> semantic_label_vector;
        if (has_mcps && mcp_vector.isValid() && _semantic_classifier) {
            semantic_label_vector = _semantic_classifier->classifyParticles(e);
        }

        std::map<int, size_t> trackid_to_index;
        if (has_mcps && mcp_vector.isValid()) {
            for (size_t i = 0; i < mcp_vector->size(); ++i) {
                trackid_to_index[mcp_vector->at(i).TrackId()] = i;
            }
        }

        std::set<art::Ptr<recob::Hit>> neutrino_hit_set(neutrino_hits.begin(), neutrino_hits.end());

        for (size_t wire_idx = 0; wire_idx < wire_vector->size(); ++wire_idx) {
            const auto& wire = wire_vector->at(wire_idx);
            auto ch_id = wire.Channel();
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
                if (neutrino_hit_set.count(hit)) {
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
                        const art::Ptr<recob::Hit>& matched_hit = filtered_hits[hit_index];
                        TruthLabelClassifier::TruthPrimaryLabel semantic_pixel_label = TruthLabelClassifier::TruthPrimaryLabel::Other;

                        if (has_mcps && mcp_vector.isValid() && _semantic_classifier) {
                            std::vector<art::Ptr<simb::MCParticle>> mcp_particles_ass_to_hit;
                            std::vector<anab::BackTrackerHitMatchingData const*> bkth_data_ass_to_hit;
                            mcp_bkth_assoc.get(matched_hit.key(), mcp_particles_ass_to_hit, bkth_data_ass_to_hit);

                            if (!bkth_data_ass_to_hit.empty()) {
                                for (size_t i_bkth = 0; i_bkth < bkth_data_ass_to_hit.size(); ++i_bkth) {
                                    if (bkth_data_ass_to_hit[i_bkth] && bkth_data_ass_to_hit[i_bkth]->isMaxIDE) {
                                        int track_id = mcp_particles_ass_to_hit[i_bkth]->TrackId();
                                        auto it_trackid = trackid_to_index.find(track_id);
                                        if (it_trackid != trackid_to_index.end()) {
                                            size_t particle_mcp_idx = it_trackid->second;
                                            if (particle_mcp_idx < semantic_label_vector.size()) {
                                                semantic_pixel_label = semantic_label_vector[particle_mcp_idx];
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }

                        if (adcs[adc_index] > _adc_image_threshold) {
                            detector_images[view_idx].set(row, col, adcs[adc_index]);
                            semantic_images[view_idx].set(row, col, static_cast<int>(semantic_pixel_label), false);
                        }
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif