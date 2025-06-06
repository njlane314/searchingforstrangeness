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

        std::vector<float> _detector_plane_image_u;
        std::vector<float> _detector_plane_image_v;
        std::vector<float> _detector_plane_image_w;

        std::vector<int> _semantic_plane_image_u;
        std::vector<int> _semantic_plane_image_v;
        std::vector<int> _semantic_plane_image_w;

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
        _tree->Branch("detector_image_u", &_detector_plane_image_u);
        _tree->Branch("detector_image_v", &_detector_plane_image_v);
        _tree->Branch("detector_image_w", &_detector_plane_image_w);
        _tree->Branch("semantic_image_u", &_semantic_plane_image_u);
        _tree->Branch("semantic_image_v", &_semantic_plane_image_v);
        _tree->Branch("semantic_image_w", &_semantic_plane_image_w);

        /*if (fEnableInference) {
            if (fProcessUView) _tree->Branch("inference_score_u", &_inference_score_u, "inference_score_u/F");
            if (fProcessVView) _tree->Branch("inference_score_v", &_inference_score_v, "inference_score_v/F");
            if (fProcessWView) _tree->Branch("inference_score_w", &_inference_score_w, "inference_score_w/F");
        }*/
    }

    void ImageAnalysis::resetTTree(TTree* ) {
        _detector_plane_image_u.clear();
        _detector_plane_image_v.clear();
        _detector_plane_image_w.clear();
        _semantic_plane_image_u.clear();
        _semantic_plane_image_v.clear();
        _semantic_plane_image_w.clear();

        /*if (fEnableInference) {
            _inference_score_u = std::numeric_limits<float>::signaling_NaN();
            _inference_score_v = std::numeric_limits<float>::signaling_NaN();
            _inference_score_w = std::numeric_limits<float>::signaling_NaN();
        }*/
    }

    void ImageAnalysis::analyseSlice(art::Event const& e, std::vector<common::ProxyPfpElem_t>& slice_pfp_v, bool , bool selected) {
        mf::LogWarning("ImageAnalysis") << "Processing slice with " << slice_pfp_v.size() << " PFPs.";
        this->resetTTree(nullptr);

        if (!selected && !slice_pfp_v.empty()) {
            mf::LogWarning("ImageAnalysis") << "Slice not selected and not empty, skipping.";
        }

        std::vector<art::Ptr<recob::Hit>> slice_hits;
        if (fCLSproducer.label().empty()) {
            mf::LogError("ImageAnalysis") << "CLSproducer InputTag is not configured. Cannot extract hits for slice.";
            return;
        }

        try {
            auto clus_handle = e.getValidHandle<std::vector<recob::Cluster>>(fCLSproducer);
            art::FindManyP<recob::Hit> cluster_hits(clus_handle, e, fCLSproducer.label());

            for (const auto& pfp_proxy : slice_pfp_v) {
                if (pfp_proxy->IsPrimary()) continue;

                for (const art::Ptr<recob::Cluster>& clus_ptr : pfp_proxy.get<recob::Cluster>()) {
                    if (clus_ptr.key() < cluster_hits.size()) {
                        const auto& clus_hit_v = cluster_hits.at(clus_ptr.key());
                        slice_hits.insert(slice_hits.end(), clus_hit_v.begin(), clus_hit_v.end());
                    } else {
                        mf::LogWarning("ImageAnalysis") << "Cluster key " << clus_ptr.key()
                                                        << " out of bounds for FindManyP<recob::Hit> (size " << cluster_hits.size() << ").";
                    }
                }
            }
            std::sort(slice_hits.begin(), slice_hits.end());
            slice_hits.erase(std::unique(slice_hits.begin(), slice_hits.end()), slice_hits.end());

        } catch (const art::Exception& ex) {
            if (ex.categoryCode() == art::errors::ProductNotFound) {
                mf::LogWarning("ImageAnalysis") << "Product not found for CLSproducer " << fCLSproducer.encode()
                                                << " or its associations. Skipping slice hit extraction.";
            } else {
                throw;
            }
            return;
        }

        if (slice_hits.empty()) {
            if (!slice_pfp_v.empty()) {
                mf::LogWarning("ImageAnalysis") << "No unique slice_hits extracted for the slice, though PFP list was not empty. Skipping image generation.";
            } else {
                mf::LogWarning("ImageAnalysis") << "Slice is empty (no PFPs, no hits). Skipping.";
            }
            return;
        }

        std::vector<Image<float>> detector_plane_tensor;
        std::vector<Image<int>> semantic_plane_tensor;
        try {
            image_algo_->generateDetectorPlaneImages(e, slice_hits, detector_plane_tensor, semantic_plane_tensor);
        } catch (const std::exception& ex) {
            mf::LogError("ImageAnalysis") << "Error generating detector plane images: " << ex.what();
            return;
        }

    
        for (size_t i = 0; i < detector_plane_tensor.size(); ++i) {
            mf::LogWarning log("ImageAnalysis");
            log << "Raw Image " << i << " (" << geo::PlaneGeo::ViewName(static_cast<geo::View_t>(i)) << ") non-zero ADC values:";
            const auto& pixels = detector_plane_tensor[i].data();
            size_t count = 0;
            for (size_t j = 0; j < pixels.size(); ++j) {
                if (pixels[j] > 0.0f) { 
                    log << " " << pixels[j];
                    count++;
                }
                if (count >= 10) break; 
            }
            if (count == 0) {
                log << " (none found)";
            } else {
                log << " (total non-zero count: " << count << ")";
            }
        }   

        if (detector_plane_tensor.size() >= 3) {
            _detector_plane_image_u = detector_plane_tensor[geo::kU].data();
            _detector_plane_image_v = detector_plane_tensor[geo::kV].data();
            _detector_plane_image_w = detector_plane_tensor[geo::kW].data();
        } else {
            mf::LogWarning("ImageAnalysis") << "generateSliceImages produced " << detector_plane_tensor.size() << " raw images, expected at least 3.";
        }

        if (semantic_plane_tensor.size() >= 3) {
            _semantic_plane_image_u = semantic_plane_tensor[geo::kU].data();
            _semantic_plane_image_v = semantic_plane_tensor[geo::kV].data();
            _semantic_plane_image_w = semantic_plane_tensor[geo::kW].data();
        } else {
            mf::LogWarning("ImageAnalysis") << "generateSliceImages produced " << semantic_plane_tensor.size() << " reco images, expected at least 3.";
        }
    }

    DEFINE_ART_CLASS_TOOL(ImageAnalysis)
}

#endif