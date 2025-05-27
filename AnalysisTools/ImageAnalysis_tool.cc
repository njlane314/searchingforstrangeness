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
#include <fstream>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TDirectoryFile.h>
#include <TVector3.h>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/System/TriggerNamesService.h"

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
        void analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool _is_data, bool selected) override;
        void setBranches(TTree* tree) override;
        void resetTTree(TTree* tree) override;

    private:
        torch::Tensor imageToTorchTensor(const Image<float>& img, geo::View_t view) const;
        std::string findFileInPath(const std::string& file, const std::string& path_env_var) const;
        void loadTorchModel(const std::string& model_path, std::shared_ptr<torch::jit::script::Module>& model, const std::string& view_name = "");

        art::InputTag fCLSproducer;
        std::string fModelFileU;
        std::string fModelFileV;
        std::string fModelFileW;
        bool fProcessUView;
        bool fProcessVView;
        bool fProcessWView;
        bool fEnableInference;

        std::unique_ptr<ImageAlgorithm> image_algo_;

        std::shared_ptr<torch::jit::script::Module> fTorchModelU;
        std::shared_ptr<torch::jit::script::Module> fTorchModelV;
        std::shared_ptr<torch::jit::script::Module> fTorchModelW;

        std::vector<float> _raw_image_u_data;
        std::vector<float> _raw_image_v_data;
        std::vector<float> _raw_image_w_data;
        std::vector<int> _reco_image_u_data;
        std::vector<int> _reco_image_v_data;
        std::vector<int> _reco_image_w_data;
        std::vector<int> _true_image_u_data;
        std::vector<int> _true_image_v_data;
        std::vector<int> _true_image_w_data;

        float _inference_score_u;
        float _inference_score_v;
        float _inference_score_w;
    };

    ImageAnalysis::ImageAnalysis(const fhicl::ParameterSet& pset) {
        this->configure(pset);
    }

    std::string ImageAnalysis::findFileInPath(const std::string& file, const std::string& path_env_var) const {
        const char* search_path_cstr = std::getenv(path_env_var.c_str());
        if (search_path_cstr) {
            std::string search_path(search_path_cstr);
            std::istringstream iss(search_path);
            std::string segment;
            while(std::getline(iss, segment, ':')) {
                if (!segment.empty() && segment.back() == '/') {
                    segment.pop_back();
                }
                std::string full_path = segment + "/" + file;
                std::ifstream f(full_path.c_str());
                if (f.good()) {
                    return full_path;
                }
            }
        }
        std::ifstream f(file.c_str());
        if (f.good()) {
            return file;
        }
        return "";
    }

    void ImageAnalysis::loadTorchModel(const std::string& model_path, std::shared_ptr<torch::jit::script::Module>& model, const std::string& view_name) {
        if (model_path.empty()) {
            throw art::Exception(art::errors::Configuration) << "Model path for " << view_name << " view is empty.";
        }
        std::string full_model_path = this->findFileInPath(model_path, "FW_SEARCH_PATH");
        if (full_model_path.empty()) {
            throw art::Exception(art::errors::NotFound)
                << "TorchScript model file not found for " << view_name << " view: " << model_path
                << ". Searched FW_SEARCH_PATH.";
        }
        try {
            model = torch::jit::load(full_model_path);
            mf::LogInfo("ImageAnalysis") << "Successfully loaded TorchScript model for " << view_name 
                                        << " view from: " << full_model_path;
        } catch (const c10::Error& e) {
            throw art::Exception(art::errors::Configuration)
                << "Error loading TorchScript model for " << view_name << " view from " 
                << full_model_path << ": " << e.what();
        }
    }

    void ImageAnalysis::configure(const fhicl::ParameterSet& p) {
        fCLSproducer = p.get<art::InputTag>("CLSproducer");
        fEnableInference = p.get<bool>("EnableInference", true);

        fhicl::ParameterSet imageAlgoParams = p.get<fhicl::ParameterSet>("ImageAlgo", p);
        image_algo_ = std::make_unique<ImageAlgorithm>(imageAlgoParams);

        if (fEnableInference) {
            fModelFileU = p.get<std::string>("ModelFileU", "");
            fModelFileV = p.get<std::string>("ModelFileV", "");
            fModelFileW = p.get<std::string>("ModelFileW", "");
            fProcessUView = p.get<bool>("ProcessUView", !fModelFileU.empty());
            fProcessVView = p.get<bool>("ProcessVView", !fModelFileV.empty());
            fProcessWView = p.get<bool>("ProcessWView", !fModelFileW.empty());

            if (fProcessUView) {
                this->loadTorchModel(fModelFileU, fTorchModelU, "U");
            }
            if (fProcessVView) {
                this->loadTorchModel(fModelFileV, fTorchModelV, "V");
            }
            if (fProcessWView) {
                this->loadTorchModel(fModelFileW, fTorchModelW, "W");
            }
        }
        mf::LogInfo("ImageAnalysis") << "Configuration complete. Inference Mode: " << (fEnableInference ? "Enabled" : "Disabled");
    }

    torch::Tensor ImageAnalysis::imageToTorchTensor(const Image<float>& img, geo::View_t view) const {
        long height = static_cast<long>(img.height());
        long width = static_cast<long>(img.width());

        if (height == 0 || width == 0) {
            mf::LogWarning("ImageAnalysis") << "View " << geo::PlaneGeo::ViewName(view)
                                           << ": Attempting to convert an empty image (0 height or width) to tensor.";
            return torch::empty({0});
        }

        auto tensor_options = torch::TensorOptions().dtype(torch::kFloat32);
        torch::Tensor tensor = torch::empty({height, width}, tensor_options);
        auto accessor = tensor.accessor<float, 2>();

        for (long r = 0; r < height; ++r) {
            for (long c = 0; c < width; ++c) {
                accessor[r][c] = img.get(r, c);
            }
        }
        tensor = tensor.reshape({1, 1, height, width});
        return tensor;
    }

    void ImageAnalysis::setBranches(TTree* _tree) {
        _tree->Branch("raw_image_u", &_raw_image_u_data);
        _tree->Branch("raw_image_v", &_raw_image_v_data);
        _tree->Branch("raw_image_w", &_raw_image_w_data);
        _tree->Branch("reco_image_u", &_reco_image_u_data);
        _tree->Branch("reco_image_v", &_reco_image_v_data);
        _tree->Branch("reco_image_w", &_reco_image_w_data);
        _tree->Branch("true_image_u", &_true_image_u_data);
        _tree->Branch("true_image_v", &_true_image_v_data);
        _tree->Branch("true_image_w", &_true_image_w_data);

        if (fEnableInference) {
            if (fProcessUView) _tree->Branch("inference_score_u", &_inference_score_u, "inference_score_u/F");
            if (fProcessVView) _tree->Branch("inference_score_v", &_inference_score_v, "inference_score_v/F");
            if (fProcessWView) _tree->Branch("inference_score_w", &_inference_score_w, "inference_score_w/F");
        }
    }

    void ImageAnalysis::resetTTree(TTree* ) {
        _raw_image_u_data.clear();
        _raw_image_v_data.clear();
        _raw_image_w_data.clear();
        _reco_image_u_data.clear();
        _reco_image_v_data.clear();
        _reco_image_w_data.clear();
        _true_image_u_data.clear();
        _true_image_v_data.clear();
        _true_image_w_data.clear();

        if (fEnableInference) {
            _inference_score_u = std::numeric_limits<float>::signaling_NaN();
            _inference_score_v = std::numeric_limits<float>::signaling_NaN();
            _inference_score_w = std::numeric_limits<float>::signaling_NaN();
        }
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool , bool selected) {
        this->resetTTree(nullptr);

        if (!selected && !slice_pfp_v.empty()) {
            mf::LogDebug("ImageAnalysis") << "Slice not selected and not empty, skipping.";
        }

        std::vector<art::Ptr<recob::Hit>> neutrino_hits;

        if (fCLSproducer.label().empty()) {
            mf::LogError("ImageAnalysis") << "CLSproducer InputTag is not configured. Cannot extract hits for slice.";
            return;
        }

        try {
            auto clus_handle = e.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
            art::FindManyP<recob::Hit> fmHits(clus_handle, e, fCLSproducer.label());

            for (const auto& pfp_proxy : slice_pfp_v) {
                if (pfp_proxy->IsPrimary()) continue;

                for (const art::Ptr<recob::Cluster>& clus_ptr : pfp_proxy.get<recob::Cluster>()) {
                    if (clus_ptr.key() < fmHits.size()) {
                        const auto& clus_hit_v = fmHits.at(clus_ptr.key());
                        neutrino_hits.insert(neutrino_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
                    } else {
                        mf::LogWarning("ImageAnalysis") << "Cluster key " << clus_ptr.key()
                                                        << " out of bounds for FindManyP<recob::Hit> (size " << fmHits.size() << ").";
                    }
                }
            }
            std::sort(neutrino_hits.begin(), neutrino_hits.end());
            neutrino_hits.erase(std::unique(neutrino_hits.begin(), neutrino_hits.end()), neutrino_hits.end());

        } catch (const art::Exception& ex) {
            if (ex.categoryCode() == art::errors::ProductNotFound) {
                mf::LogWarning("ImageAnalysis") << "Product not found for CLSproducer " << fCLSproducer.encode()
                                                << " or its associations. Skipping slice hit extraction.";
            } else {
                throw;
            }
            return;
        }

        if (neutrino_hits.empty()) {
            if (!slice_pfp_v.empty()) {
                mf::LogDebug("ImageAnalysis") << "No unique neutrino_hits extracted for the slice, though PFP list was not empty. Skipping image generation.";
            } else {
                mf::LogDebug("ImageAnalysis") << "Slice is empty (no PFPs, no hits). Skipping.";
            }
            return;
        }

        std::vector<Image<float>> out_raw_images;
        std::vector<Image<int>> out_reco_images;
        std::vector<Image<int>> out_true_images;
        image_algo_->generateSliceImages(e, neutrino_hits, out_raw_images, out_reco_images, out_true_images);

        if (out_raw_images.size() >= 3) {
            _raw_image_u_data = out_raw_images[geo::kU].data();
            _raw_image_v_data = out_raw_images[geo::kV].data();
            _raw_image_w_data = out_raw_images[geo::kW].data();
        } else {
            mf::LogWarning("ImageAnalysis") << "generateSliceImages produced " << out_raw_images.size() << " raw images, expected at least 3.";
        }

        if (out_reco_images.size() >= 3) {
            _reco_image_u_data = out_reco_images[geo::kU].data();
            _reco_image_v_data = out_reco_images[geo::kV].data();
            _reco_image_w_data = out_reco_images[geo::kW].data();
        } else {
            mf::LogWarning("ImageAnalysis") << "generateSliceImages produced " << out_reco_images.size() << " reco images, expected at least 3.";
        }

        if (out_true_images.size() >= 3) {
            _true_image_u_data = out_true_images[geo::kU].data();
            _true_image_v_data = out_true_images[geo::kV].data();
            _true_image_w_data = out_true_images[geo::kW].data();
        } else {
            mf::LogWarning("ImageAnalysis") << "generateSliceImages produced " << out_true_images.size() << " true images, expected at least 3.";
        }

        if (fEnableInference) {
            // U-View Inference
            if (fProcessUView && fTorchModelU && out_raw_images.size() > geo::kU) {
                torch::Tensor input_tensor_u = imageToTorchTensor(out_raw_images[geo::kU], geo::kU);
                if (input_tensor_u.numel() > 0) {
                    std::vector<torch::jit::IValue> inputs_u;
                    inputs_u.push_back(input_tensor_u.to(torch::kCPU));
                    try {
                        // Set the U model to evaluation mode
                        at::Tensor output_u = fTorchModelU->forward(inputs_u).toTensor();
                        if (output_u.numel() == 1) _inference_score_u = output_u.item<float>();
                        else mf::LogWarning("ImageAnalysis") << "U-View model output has " << output_u.numel() << " elements, expected 1 for score.";
                    } catch (const c10::Error& err) {
                        mf::LogError("ImageAnalysis") << "Torch U-View Inference Error: " << err.what();
                    }
                }
            }

            // V-View Inference
            if (fProcessVView && fTorchModelV && out_raw_images.size() > geo::kV) {
                torch::Tensor input_tensor_v = imageToTorchTensor(out_raw_images[geo::kV], geo::kV);
                if (input_tensor_v.numel() > 0) {
                    std::vector<torch::jit::IValue> inputs_v;
                    inputs_v.push_back(input_tensor_v.to(torch::kCPU));
                    try {
                        // Set the V model to evaluation mode
                        at::Tensor output_v = fTorchModelV->forward(inputs_v).toTensor();
                        if (output_v.numel() == 1) _inference_score_v = output_v.item<float>();
                        else mf::LogWarning("ImageAnalysis") << "V-View model output has " << output_v.numel() << " elements, expected 1 for score.";
                    } catch (const c10::Error& err) {
                        mf::LogError("ImageAnalysis") << "Torch V-View Inference Error: " << err.what();
                    }
                }
            }

            // W-View Inference
            if (fProcessWView && fTorchModelW && out_raw_images.size() > geo::kW) {
                torch::Tensor input_tensor_w = imageToTorchTensor(out_raw_images[geo::kW], geo::kW);
                if (input_tensor_w.numel() > 0) {
                    std::vector<torch::jit::IValue> inputs_w;
                    inputs_w.push_back(input_tensor_w.to(torch::kCPU));
                    try {
                        // Set the W model to evaluation mode
                        at::Tensor output_w = fTorchModelW->forward(inputs_w).toTensor();
                        if (output_w.numel() == 1) _inference_score_w = output_w.item<float>();
                        else mf::LogWarning("ImageAnalysis") << "W-View model output has " << output_w.numel() << " elements, expected 1 for score.";
                    } catch (const c10::Error& err) {
                        mf::LogError("ImageAnalysis") << "Torch W-View Inference Error: " << err.what();
                    }
                }
            }
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif